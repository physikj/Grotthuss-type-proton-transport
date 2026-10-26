%% H + H₂ → H₃ 各向异性反应截面分析 (高性能版)
%  从VASP AIMD轨迹中识别 H + H₂ 碰撞事件，
%  统计H₃形成截面随碰撞角度的变化 σ(θ)
%
%  反应判据: 游离H进入H₂的bond_threshold内并持续≥2帧 → H₃稳定形成
%
%  性能优化:
%    - 预计算分数坐标，避免重复矩阵乘法
%    - 矩阵存储+sortrows代替containers.Map动态cell数组
%    - 内联碰撞参数计算，消除函数调用开销

clear; close all; clc;

%% ==================== 1. 参数设置 ====================
% --- 文件与体系 ---
filename      = 'XDATCAR_60_2500_4';
n_H           = 432;
n_O           = 72;
n_atom        = n_H + n_O;

% --- 时间步长 ---
timestep       = 0.5;              % fs (VASP POTIM)
output_interval = 10;              % 每N步输出
frame_dt       = timestep * output_interval;

% --- 距离判据 (Å, 带滞回) ---
r_HH_min = 0.50;                   % H₂键长下限 (排除噪声)
r_HH_on  = 1.00;                   % H₂形成阈值
r_HH_off = 1.20;                   % H₂断裂阈值
r_H3_on  = 1.00;                   % H₃形成阈值
r_H3_off = 1.20;                   % H₃分离阈值
collision_threshold = 4.0;         % 碰撞检测半径
b_max              = collision_threshold;

% --- 稳定性要求 ---
min_h3_frames = 1;                 % H₃至少连续N帧才算反应

% --- 读取控制 ---
max_frames = 100000;

% --- 角度分箱 (0-90°, 同核双原子) ---
theta_edges   = 0:10:90;
theta_centers = (theta_edges(1:end-1) + theta_edges(2:end)) / 2;
n_theta_bins  = length(theta_centers);

% --- 能量截断 (eV) ---
E_min = 0.001;
E_max = 10.0;

% --- 矩阵预分配 (每帧碰撞对估计) ---
PAIRS_PER_FRAME_EST = 500;         % 粗略估计，超出自动扩展

%% ==================== 2. 读取XDATCAR (预计算分数坐标) ====================
fprintf('══════════════════════════════════════════\n');
fprintf('  H + H₂ → H₃  各向异性反应截面 (高性能版)\n');
fprintf('══════════════════════════════════════════\n\n');
fprintf('[1/5] 读取并预计算分数坐标...\n');

fp = fopen(filename, 'r');
assert(fp ~= -1, '无法打开文件: %s', filename);

fgetl(fp);
scale = sscanf(fgetl(fp), '%f', 1);

lat = zeros(3, 3);
for i = 1:3
    lat(i,:) = sscanf(fgetl(fp), '%f %f %f')';
end
lat = lat * scale;

fgetl(fp);
counts = sscanf(fgetl(fp), '%d %d');
fprintf('  原子: H=%d, O=%d\n', counts(1), counts(2));

% 跳转到坐标
tline = fgetl(fp);
safety = 0;
while ischar(tline) && ~startsWith(strtrim(tline), 'Direct') ...
                   && ~startsWith(strtrim(tline), 'Cartesian')
    tline = fgetl(fp);
    safety = safety + 1;
    assert(safety < 200, '未找到Direct/Cartesian行');
end

% 读取所有帧 → 同时存笛卡尔和分数坐标
frames_cart = {};                   % 笛卡尔 (Å)
frames_frac = {};                   % 分数坐标 (预计算)
n_read = 0;

while ischar(tline) && n_read < max_frames
    tag = strtrim(tline);
    if startsWith(tag, 'Direct') || startsWith(tag, 'Cartesian')
        n_read = n_read + 1;
        crd = zeros(n_atom, 3);
        for a = 1:n_atom
            tline = fgetl(fp);
            if ~ischar(tline), break; end
            crd(a,:) = sscanf(tline, '%f %f %f')';
        end
        cart = crd * lat;
        frames_cart{n_read} = cart;
        frames_frac{n_read} = crd;  % 直接存分数坐标
    end
    tline = fgetl(fp);
end
fclose(fp);

n_frames = n_read;
fprintf('  读取 %d 帧  (dt = %.1f fs)\n', n_frames, frame_dt);
assert(n_frames >= 10, '帧数太少');

%% ==================== 3. 逐帧识别 H₂/H₃ (带滞回) + 碰撞对 ====================
fprintf('[2/5] 逐帧识别H₂/H₃ (滞回: H₂ %g→%g, H₃ %g→%g) + 碰撞对...\n', ...
    r_HH_on, r_HH_off, r_H3_on, r_H3_off);

% 预分配碰撞对矩阵
% 列: [encoded_key, frame, h1, h2, free_H, is_H3]
max_pairs = PAIRS_PER_FRAME_EST * n_frames;
pair_mat  = zeros(max_pairs, 6);
n_pairs   = 0;

bond_lens = [];

% 跨帧持久状态 (滞回)
h2_bonds  = containers.Map('KeyType','char','ValueType','any');  % 当前H₂键
h3_active = containers.Map('KeyType','char','ValueType','any');  % 当前活跃H₃

for fi = 1:n_frames
    if mod(fi, 100) == 0
        fprintf('  帧 %4d / %d  (已收集 %d 对)\n', fi, n_frames, n_pairs);
    end

    frac = frames_frac{fi};

    % ---- 1. H₂键 (滞回) ----
    bonded_atoms = false(n_H, 1);
    new_h2 = containers.Map('KeyType','char','ValueType','any');

    % 1a. 延续上一帧H₂键 (断裂阈值 r_HH_off)
    old_keys = keys(h2_bonds);
    for k = 1:length(old_keys)
        key  = old_keys{k};
        pair = h2_bonds(key);
        i = pair(1); j = pair(2);
        if pbc_dist(frac, i, j, lat) <= r_HH_off
            bonded_atoms(i) = true;
            bonded_atoms(j) = true;
            new_h2(key) = pair;
        end
    end

    % 1b. 形成新H₂键 (形成阈值 r_HH_on, 贪心)
    for i = 1:n_H
        if bonded_atoms(i), continue; end
        best_j = 0; best_dist = inf;
        for j = i+1:n_H
            if bonded_atoms(j), continue; end
            d = pbc_dist(frac, i, j, lat);
            if d < best_dist
                best_dist = d; best_j = j;
            end
        end
        if best_j > 0 && best_dist >= r_HH_min && best_dist <= r_HH_on
            bonded_atoms(i) = true;
            bonded_atoms(best_j) = true;
            new_h2(sprintf('%d_%d', i, best_j)) = [i, best_j];
        end
    end

    h2_bonds = new_h2;

    % 收集H₂键长
    h2_keys = keys(h2_bonds);
    for k = 1:length(h2_keys)
        pair = h2_bonds(h2_keys{k});
        bond_lens(end+1) = pbc_dist(frac, pair(1), pair(2), lat);
    end

    % 成键查找表
    bond_partner = zeros(n_H, 1);
    for k = 1:length(h2_keys)
        pair = h2_bonds(h2_keys{k});
        bond_partner(pair(1)) = pair(2);
        bond_partner(pair(2)) = pair(1);
    end

    % ---- 2. H₃ (滞回) ----
    new_h3 = containers.Map('KeyType','char','ValueType','any');
    used_freeH = false(n_H, 1);

    % 2a. 延续上一帧H₃ (分离阈值 r_H3_off), 重新推导角色
    old_h3_keys = keys(h3_active);
    for k = 1:length(old_h3_keys)
        key   = old_h3_keys{k};
        atoms = h3_active(key);
        a = atoms(1); b = atoms(2); c = atoms(3);
        if bond_partner(a) == b
            h1 = a; h2 = b; fh = c;
        elseif bond_partner(a) == c
            h1 = a; h2 = c; fh = b;
        elseif bond_partner(b) == c
            h1 = b; h2 = c; fh = a;
        else
            continue;
        end
        d_min = min(pbc_dist(frac, h1, fh, lat), pbc_dist(frac, h2, fh, lat));
        if d_min <= r_H3_off
            new_h3(key) = [h1, h2, fh];
            used_freeH(fh) = true;
        end
    end

    % 2b. 形成新H₃ (形成阈值 r_H3_on)
    for k = 1:length(h2_keys)
        pair = h2_bonds(h2_keys{k});
        h1 = pair(1); h2 = pair(2);
        already = false;
        nk = keys(new_h3);
        for kk = 1:length(nk)
            rl = new_h3(nk{kk});
            if (rl(1)==h1 && rl(2)==h2) || (rl(1)==h2 && rl(2)==h1)
                already = true; break;
            end
        end
        if already, continue; end
        for at = 1:n_H
            if at == h1 || at == h2, continue; end
            if bonded_atoms(at), continue; end
            if used_freeH(at), continue; end
            d_min = min(pbc_dist(frac, h1, at, lat), pbc_dist(frac, h2, at, lat));
            if d_min <= r_H3_on
                s = sort([h1, h2, at]);
                new_h3(sprintf('%d_%d_%d', s(1), s(2), s(3))) = [h1, h2, at];
                used_freeH(at) = true;
            end
        end
    end

    h3_active = new_h3;

    % ---- 3. 碰撞对收集 (碰撞半径 collision_threshold) ----
    free_H = find(~bonded_atoms);
    for h = 1:length(h2_keys)
        pair = h2_bonds(h2_keys{h});
        h1 = pair(1); h2 = pair(2);
        for f = 1:length(free_H)
            fh = free_H(f);

            % 碰撞半径判断 (游离H到H₂中心)
            d_frac = frac(fh,:) - (frac(h1,:) + frac(h2,:))/2;
            d_frac = d_frac - round(d_frac);
            d_center = norm(d_frac * lat);
            if d_center > collision_threshold, continue; end

            % 是否H₃ (用h3_active规范key判断, 带滞回)
            s = sort([h1, h2, fh]);
            is_H3 = isKey(h3_active, sprintf('%d_%d_%d', s(1), s(2), s(3)));

            n_pairs = n_pairs + 1;
            if n_pairs > size(pair_mat, 1)
                pair_mat = [pair_mat; zeros(max_pairs, 6)];
            end
            encoded = s(1)*1000000 + s(2)*1000 + s(3);
            pair_mat(n_pairs,:) = [encoded, fi, h1, h2, fh, double(is_H3)];
        end
    end
end

% 裁剪到实际大小
pair_mat = pair_mat(1:n_pairs, :);
fprintf('  完成: %d 个碰撞对, %d 个H₂键长\n', n_pairs, length(bond_lens));

% 📊 诊断: is_H3统计
n_is_H3 = sum(pair_mat(:,6));
fprintf('  📊 其中 is_H3=1: %d 对 (%.1f%%)\n', n_is_H3, 100*n_is_H3/max(1,n_pairs));

% 📊 诊断: 多少帧出现了H₃
h3_frames = unique(pair_mat(pair_mat(:,6)==1, 2));
fprintf('  📊 H₃出现的帧数: %d (共%d帧分析)\n', length(h3_frames), n_frames);

if n_pairs == 0
    fprintf('\n❌ 未找到任何碰撞对！\n');
    return;
end

%% ==================== 4. 排序 + 找连续段 → 碰撞遭遇 ====================
fprintf('[3/5] 排序碰撞对并提取碰撞遭遇...\n');

% 按 encoded_key 排序，同key内按帧号排序
pair_mat = sortrows(pair_mat, [1, 2]);

% 单遍扫描找连续段 (same key, consecutive frames)
% events: [b, theta, E_trans, is_reaction]
events    = zeros(ceil(n_pairs/2), 5);  % 预分配
n_events  = 0;

% 📊 诊断计数器
diag_encounters   = 0;   % 总遭遇数
diag_has_h3       = 0;   % 遭遇中有is_H3=1的帧
diag_has_stable_h3 = 0;  % 遭遇中有≥2连续is_H3=1
diag_zero_vel     = 0;   % 因零速度跳过
diag_energy_skip  = 0;   % 因能量截断跳过

idx = 1;
while idx <= n_pairs
    % 当前段的key和起始
    key_start  = pair_mat(idx, 1);
    seg_start  = idx;
    seg_frames = pair_mat(idx, 2);
    seg_is_H3  = pair_mat(idx, 6);

    % 向前扫描直到key变或帧不连续
    idx = idx + 1;
    while idx <= n_pairs && ...
          pair_mat(idx, 1) == key_start && ...
          pair_mat(idx, 2) == pair_mat(idx-1, 2) + 1
        seg_frames(end+1) = pair_mat(idx, 2);
        seg_is_H3(end+1)  = pair_mat(idx, 6);
        idx = idx + 1;
    end

    seg_len = length(seg_frames);

    % 📊 诊断
    diag_encounters = diag_encounters + 1;
    if any(seg_is_H3 == 1)
        diag_has_h3 = diag_has_h3 + 1;
    end

    % ---- 检查是否形成稳定H₃ (≥2连续帧的is_H3=1) ----
    is_reaction = false;
    h3_run = 0;
    for s = 1:seg_len
        if seg_is_H3(s) == 1
            h3_run = h3_run + 1;
            if h3_run >= min_h3_frames
                is_reaction = true;
                break;
            end
        else
            h3_run = 0;
        end
    end

    if is_reaction
        diag_has_stable_h3 = diag_has_stable_h3 + 1;
    end

    % ---- 在首帧计算碰撞参数 ----
    fs = pair_mat(seg_start, 2);        % 首帧号
    h1 = pair_mat(seg_start, 3);        % H₂_H1
    h2 = pair_mat(seg_start, 4);        % H₂_H2
    fh = pair_mat(seg_start, 5);        % 游离H索引

    % 速度 (首帧, PBC修正)
    if fs < 2
        % 前向差分
        frac_t  = frames_frac{fs};
        frac_t1 = frames_frac{fs+1};
        dt_sign = 1;
    else
        % 后向差分
        frac_t  = frames_frac{fs};
        frac_t1 = frames_frac{fs-1};
        dt_sign = -1;
    end

    % 自由H速度
    d_frac = frac_t(fh,:) - frac_t1(fh,:);
    d_frac = d_frac - round(d_frac);
    v_free = dt_sign * (d_frac * lat) / frame_dt;

    % H₂两个H的速度
    d_frac = frac_t(h1,:) - frac_t1(h1,:);
    d_frac = d_frac - round(d_frac);
    v_h1 = dt_sign * (d_frac * lat) / frame_dt;

    d_frac = frac_t(h2,:) - frac_t1(h2,:);
    d_frac = d_frac - round(d_frac);
    v_h2 = dt_sign * (d_frac * lat) / frame_dt;

    v_H2   = (v_h1 + v_h2) / 2;
    v_rel  = v_free - v_H2;
    v_mag  = norm(v_rel);

    if v_mag < 1e-6
        diag_zero_vel = diag_zero_vel + 1;
        continue;  % 零速度，跳过
    end

    % 碰撞能量 E_trans (eV)
    m_H  = 1.00784;                    % amu
    mu   = (m_H * 2*m_H) / (m_H + 2*m_H);
    E_trans = 0.5 * mu * 1.660539e-27 * (v_mag*1e5)^2 / 1.602177e-19;

    if E_trans < E_min || E_trans > E_max
        diag_energy_skip = diag_energy_skip + 1;
        continue;
    end

    % 碰撞角度 θ (同核双原子, 0-90°)
    % H₂键轴 (从h2指向h1)
    d_frac = frac_t(h1,:) - frac_t(h2,:);
    d_frac = d_frac - round(d_frac);
    ax = (d_frac * lat) / norm(d_frac * lat);

    v_unit   = v_rel / v_mag;
    cos_theta = abs(dot(ax, v_unit));
    cos_theta = min(max(cos_theta, 0), 1);
    theta = acosd(cos_theta);

    % 碰撞参数 b
    cent_frac = (frac_t(h1,:) + frac_t(h2,:)) / 2;
    d_frac = frac_t(fh,:) - cent_frac;
    d_frac = d_frac - round(d_frac);
    r_rel_vec = d_frac * lat;
    b = norm(cross(r_rel_vec, v_unit));

    % 存储事件
    n_events = n_events + 1;
    if n_events > size(events, 1)
        events = [events; zeros(ceil(n_pairs/4), 5)];
    end
    events(n_events,:) = [b, theta, E_trans, double(is_reaction), seg_len];
end

events = events(1:n_events, :);
fprintf('  碰撞遭遇: %d  |  稳定H₃反应: %d  (%.1f%%)\n', ...
    n_events, sum(events(:,4)), 100*sum(events(:,4))/max(1,n_events));

% 📊 诊断报告
fprintf('\n  ╔═══════════════════════════════════╗\n');
fprintf('  ║  📊 数据流诊断                   ║\n');
fprintf('  ╠═══════════════════════════════════╣\n');
fprintf('  ║  总碰撞对(帧×对) : %8d      ║\n', n_pairs);
fprintf('  ║  其中 is_H₃=1    : %8d (%4.1f%%) ║\n', n_is_H3, 100*n_is_H3/max(1,n_pairs));
fprintf('  ║  H₃出现帧数      : %8d          ║\n', length(h3_frames));
fprintf('  ╠═══════════════════════════════════╣\n');
fprintf('  ║  总遭遇(连续段)  : %8d          ║\n', diag_encounters);
fprintf('  ║  有H₃帧的遭遇   : %8d (%4.1f%%) ║\n', diag_has_h3, 100*diag_has_h3/max(1,diag_encounters));
fprintf('  ║  稳定H₃(≥%d帧)  : %8d (%4.1f%%) ║\n', min_h3_frames, diag_has_stable_h3, 100*diag_has_stable_h3/max(1,diag_encounters));
fprintf('  ║  零速度跳过      : %8d          ║\n', diag_zero_vel);
fprintf('  ║  能量截断跳过    : %8d          ║\n', diag_energy_skip);
fprintf('  ║  最终有效事件    : %8d          ║\n', n_events);
fprintf('  ╚═══════════════════════════════════╝\n');

if n_events == 0
    fprintf('\n❌ 无有效碰撞事件。\n');
    return;
end

%% ==================== 5. 角度分辨截面 ====================
fprintf('[4/5] 计算各向异性截面 σ(θ)...\n');

sigma_by_theta = zeros(n_theta_bins, 7);

for ib = 1:n_theta_bins
    t_lo = theta_edges(ib);
    t_hi = theta_edges(ib+1);
    in_bin = events(:,2) >= t_lo & events(:,2) < t_hi;
    N_total = sum(in_bin);
    N_r = sum(events(in_bin, 4));

    if N_total > 0
        P_r = N_r / N_total;
        P_r_err = sqrt(P_r * (1-P_r) / N_total);
    else
        P_r = 0;
        P_r_err = 0;
    end
    sigma     = P_r * pi * b_max^2;
    sigma_err = P_r_err * pi * b_max^2;
    sigma_by_theta(ib,:) = [theta_centers(ib), sigma, sigma_err, ...
                            P_r, P_r_err, N_total, N_r];
end

P_r_total    = sum(events(:,4)) / n_events;
sigma_total  = P_r_total * pi * b_max^2;
sigma_total_err = sqrt(P_r_total*(1-P_r_total)/n_events) * pi * b_max^2;

fprintf('  总截面: σ = %.6f ± %.6f Å²\n', sigma_total, sigma_total_err);
fprintf('  (计算: ΣN_rxn=%d / N_coll=%d = %.6f × π×%.1f²=%.2f)\n', ...
    sum(events(:,4)), n_events, P_r_total, b_max, pi*b_max^2);
fprintf('\n  角度分辨:\n');
fprintf('    θ(°)   N_coll  N_rxn   P_r     σ(Å²)\n');
fprintf('    --------------------------------------\n');
for ib = 1:n_theta_bins
    if sigma_by_theta(ib,6) > 0
        fprintf('    %2.0f-%2.0f  %5d   %4d   %.3f   %.3f\n', ...
            theta_edges(ib), theta_edges(ib+1), ...
            sigma_by_theta(ib,6), sigma_by_theta(ib,7), ...
            sigma_by_theta(ib,4), sigma_by_theta(ib,2));
    end
end

%% ==================== 6. 可视化 ====================
fprintf('[5/5] 绘图并导出...\n');

c_rxn   = [0.85 0.20 0.15];
c_norxn = [0.25 0.45 0.75];
c_sig   = [0.90 0.40 0.10];

figure('Position', [30 30 1500 900]);

% (a) σ(θ)
subplot(2,3,1);
has = sigma_by_theta(:,6) > 0;
errorbar(sigma_by_theta(has,1), sigma_by_theta(has,2), ...
         sigma_by_theta(has,3), 'o-', 'Color', c_sig, ...
         'LineWidth', 2.5, 'MarkerSize', 10, ...
         'MarkerFaceColor', c_sig, 'CapSize', 6);
xlabel('碰撞角度 θ (°)', 'FontSize', 13);
ylabel('反应截面 σ(θ) (Å²)', 'FontSize', 13);
title(sprintf('各向异性反应截面 (σ_{tot}=%.2f Å²)', sigma_total), 'FontSize', 14);
grid on; xlim([0 90]); yl = ylim; ylim([0 yl(2)]);

% (b) P_r(θ)
subplot(2,3,2);
errorbar(sigma_by_theta(has,1), sigma_by_theta(has,4), ...
         sigma_by_theta(has,5), 's-', 'Color', [0.2 0.55 0.2], ...
         'LineWidth', 2.5, 'MarkerSize', 10, ...
         'MarkerFaceColor', [0.2 0.55 0.2], 'CapSize', 6);
xlabel('碰撞角度 θ (°)', 'FontSize', 13);
ylabel('反应概率 P_r(θ)', 'FontSize', 13);
title('反应概率 vs 碰撞角度', 'FontSize', 14);
grid on; xlim([0 90]); ylim([0 1]);

% (c) 角度分布
subplot(2,3,3);
idx_r = events(:,4) == 1; idx_n = events(:,4) == 0;
th_r = events(idx_r, 2); th_n = events(idx_n, 2);
h1 = histogram(th_n, theta_edges, 'FaceColor', c_norxn, ...
               'EdgeColor', 'none', 'FaceAlpha', 0.65); hold on;
h2 = histogram(th_r, theta_edges, 'FaceColor', c_rxn, ...
               'EdgeColor', 'none', 'FaceAlpha', 0.65);
xlabel('碰撞角度 θ (°)', 'FontSize', 13); ylabel('事件数', 'FontSize', 13);
title('碰撞角度分布', 'FontSize', 14);
legend([h1 h2], {sprintf('非反应 (N=%d)', sum(idx_n)), ...
                 sprintf('稳定H₃ (N=%d)', sum(idx_r))}, ...
       'Location', 'northwest', 'FontSize', 10);
grid on; xlim([0 90]);

% (d) (θ, E) 散点
subplot(2,3,4);
scatter(events(idx_n,2), events(idx_n,3), 12, c_norxn, ...
        'filled', 'MarkerFaceAlpha', 0.4); hold on;
scatter(events(idx_r,2), events(idx_r,3), 18, c_rxn, ...
        'filled', 'MarkerFaceAlpha', 0.7);
xlabel('碰撞角度 θ (°)', 'FontSize', 13);
ylabel('碰撞能量 E_{trans} (eV)', 'FontSize', 13);
title('角度-能量分布', 'FontSize', 14);
legend('非反应', '稳定H₃', 'Location', 'best', 'FontSize', 10);
grid on; xlim([0 90]);

% (e) b分布
subplot(2,3,5);
b_r = events(idx_r,1); b_n = events(idx_n,1);
b_ed = linspace(0, b_max, 20);
histogram(b_n, b_ed, 'FaceColor', c_norxn, 'EdgeColor', 'none', 'FaceAlpha', 0.65); hold on;
histogram(b_r, b_ed, 'FaceColor', c_rxn, 'EdgeColor', 'none', 'FaceAlpha', 0.65);
xlabel('碰撞参数 b (Å)', 'FontSize', 13); ylabel('事件数', 'FontSize', 13);
title('碰撞参数分布', 'FontSize', 14);
legend('非反应', '稳定H₃', 'Location', 'northwest', 'FontSize', 10);
grid on;

% (f) 遭遇寿命
subplot(2,3,6);
histogram(events(:,5), 'BinMethod', 'integers', ...
          'FaceColor', [0.5 0.4 0.7], 'EdgeColor', 'k');
xlabel('遭遇持续帧数', 'FontSize', 13); ylabel('事件数', 'FontSize', 13);
title(sprintf('遭遇寿命 (均值 %.1f帧)', mean(events(:,5))), 'FontSize', 14);
grid on;

sgtitle(sprintf(['H + H_2 → H_3 各向异性截面  ' ...
       '(σ_{tot}=%.2f±%.2f Å², %d碰撞, %d反应, ≥%d帧稳定H₃)'], ...
       sigma_total, sigma_total_err, n_events, sum(idx_r), min_h3_frames), ...
       'FontSize', 15, 'FontWeight', 'bold');

saveas(gcf, 'Anisotropic_Cross_Section.png');
fprintf('  图表: Anisotropic_Cross_Section.png\n');

%% ==================== 7. 导出 ====================
T = table(theta_centers', sigma_by_theta(:,2), sigma_by_theta(:,3), ...
          sigma_by_theta(:,4), sigma_by_theta(:,5), ...
          sigma_by_theta(:,6), sigma_by_theta(:,7), ...
    'VariableNames', {'Theta_deg','Sigma_A2','Sigma_err','P_r','P_r_err', ...
                      'N_coll','N_rxn'});
writetable(T, 'anisotropic_cross_section.csv');

E = array2table(events, 'VariableNames', ...
    {'b_A','theta_deg','E_eV','is_reaction','duration'});
writetable(E, 'collision_events.csv');
fprintf('  数据: anisotropic_cross_section.csv, collision_events.csv\n');

%% ==================== 汇总 ====================
fprintf('\n');
fprintf('══════════════════════════════════════════\n');
fprintf('              分析完成\n');
fprintf('══════════════════════════════════════════\n');
fprintf('  帧数            : %d  (dt=%.1f fs)\n', n_frames, frame_dt);
fprintf('  H₂平均键长      : %.4f Å\n', mean(bond_lens));
fprintf('  碰撞遭遇数      : %d\n', n_events);
fprintf('  稳定H₃反应数    : %d  (≥%d帧, %.1f%%)\n', ...
        sum(events(:,4)), min_h3_frames, 100*sum(events(:,4))/n_events);
fprintf('  总反应截面      : %.6f ± %.6f Å²\n', sigma_total, sigma_total_err);
fprintf('══════════════════════════════════════════\n');

%% ==================== 辅助函数 ====================
function d = pbc_dist(frac, i, j, lat)
% 最小镜像距离: 原子i与j在分数坐标frac下的PBC距离
    v = frac(i,:) - frac(j,:);
    v = v - round(v);
    d = norm(v * lat);
end
