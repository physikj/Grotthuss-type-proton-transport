%% H₃分子内部RDF分析 (含H交换筛选)
%  从VASP XDATCAR中识别稳定H₃复合体，
%  筛选发生H交换的事件 (H + H₂ → H₃ → H' + H₂')，
%  统计内部H-H距离分布 (按短/中/长分类) 及键角分布

clear; close all; clc;

%% ==================== 参数 ====================
fpath = 'XDATCAR_60_2500_3';

% H₂识别 (滞回: 形成<r_HH_on, 断裂>r_HH_off)
r_HH_min = 0.50;      % H₂键长下限 (排除数值噪声)
r_HH_on  = 1.0;      % 形成阈值: d < 1.0 识别为H₂
r_HH_off = 1.2;      % 断裂阈值: d > 1.2 才认为分离

% H₃识别 (滞回: 形成<r_H3_on, 分离>r_H3_off)
r_H3_on  = 1.0;      % 形成阈值: 游离H < 1.2 识别为H₃
r_H3_off = 1.2;      % 分离阈值: 游离H > 1.5 才认为分离
min_h3_frames = 10;    % H₃最少持续帧数

% 分析窗口
n_frames = 20000;            % 读取最后N帧

% 交换检测
post_window = 5;             % H₃事件结束后搜索窗口 (帧)

% PDF参数
dr   = 0.04;
r_max = 5.0;
r_edges  = 0:dr:r_max;
r_centers = (r_edges(1:end-1) + r_edges(2:end)) / 2;

% 键角参数
dtheta = 5;
ang_edges = 0:dtheta:180;
ang_centers = (ang_edges(1:end-1) + ang_edges(2:end)) / 2;

%% ==================== 读取XDATCAR ====================
fprintf('读取 %s ...\n', fpath);

fp = fopen(fpath, 'r');
assert(fp ~= -1, '无法打开文件: %s', fpath);

fgetl(fp);
scale = fscanf(fp, '%f', 1);
fgetl(fp);

lat = zeros(3, 3);
for i = 1:3
    lat(i,:) = sscanf(fgetl(fp), '%f %f %f')';
end

fgetl(fp);
counts = sscanf(fgetl(fp), '%d %d');
n_H = counts(1);
n_O = counts(2);
n_atom = n_H + n_O;

fprintf('  H: %d, O: %d\n', n_H, n_O);

% 定位坐标行
tline = fgetl(fp);
safety = 0;
while ischar(tline) && ~contains(tline, 'Direct') && ~contains(tline, 'Cartesian')
    tline = fgetl(fp);
    safety = safety + 1;
    assert(safety < 200, '未找到坐标行');
end

%% ==================== 读取所有帧 ====================
frames = {};
while ischar(tline)
    if contains(tline, 'Direct') || contains(tline, 'Cartesian')
        crd = zeros(n_atom, 3);
        for a = 1:n_atom
            tline = fgetl(fp);
            if ~ischar(tline), break; end
            crd(a,:) = sscanf(tline, '%f %f %f')';
        end
        frames{end+1} = crd * lat;
    end
    tline = fgetl(fp);
end
fclose(fp);

n_total = length(frames);
start_frame = max(1, n_total - n_frames + 1);
n_analyzed = n_total - start_frame + 1;
fprintf('  总帧: %d, 分析最后 %d 帧 (帧 %d → %d)\n', ...
    n_total, n_analyzed, start_frame, n_total);

%% ==================== 逐帧识别H₂ 和 H₃ (带滞回) ====================
fprintf('逐帧识别H₂和H₃ (滞回: H₂ %g→%g, H₃ %g→%g) ...\n', ...
    r_HH_on, r_HH_off, r_H3_on, r_H3_off);

frame_h3 = cell(n_analyzed, 1);
bond_partner_all = zeros(n_H, n_analyzed);   % 每帧成键查找表 (供回溯上一帧定H₂身份)

% 跨帧持久状态
h2_bonds  = containers.Map('KeyType','char','ValueType','any');  % 当前H₂键: key→[i,j]
h3_active = containers.Map('KeyType','char','ValueType','any');  % 当前活跃H₃: key→[h1,h2,fh]

for fi = 1:n_analyzed
    if mod(fi, 200) == 0
        fprintf('  帧 %d/%d\n', fi, n_analyzed);
    end

    frac = frames{start_frame + fi - 1} / lat;

    % ========== 1. 更新H₂键 (滞回) ==========
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

    % 建立成键查找表: bond_partner(i) = j 若 i-j 成键, 否则 0
    bond_partner = zeros(n_H, 1);
    h2_keys_all = keys(h2_bonds);
    for k = 1:length(h2_keys_all)
        pair = h2_bonds(h2_keys_all{k});
        bond_partner(pair(1)) = pair(2);
        bond_partner(pair(2)) = pair(1);
    end
    bond_partner_all(:, fi) = bond_partner;   % 存下本帧成键关系

    % ========== 2. 更新H₃ (滞回) ==========
    new_h3 = containers.Map('KeyType','char','ValueType','any');
    used_freeH = false(n_H, 1);          % 游离H不可重复归属

    % 2a. 延续上一帧H₃ (分离阈值 r_H3_off), 并按当前H₂键重新推导角色
    old_h3_keys = keys(h3_active);
    for k = 1:length(old_h3_keys)
        key   = old_h3_keys{k};
        atoms = h3_active(key);          % [a, b, c] 三个原子 (旧角色仅用于取集合)
        a = atoms(1); b = atoms(2); c = atoms(3);

        % 确定三个原子中当前哪两个成键
        if bond_partner(a) == b
            h1 = a; h2 = b; fh = c;
        elseif bond_partner(a) == c
            h1 = a; h2 = c; fh = b;
        elseif bond_partner(b) == c
            h1 = b; h2 = c; fh = a;
        else
            continue;                     % 无成键对 → H₃已解散
        end

        % 游离H仍在分离阈值内才延续
        d_min = min(pbc_dist(frac, h1, fh, lat), pbc_dist(frac, h2, fh, lat));
        if d_min <= r_H3_off
            new_h3(key) = [h1, h2, fh];   % 用重新推导的角色
            used_freeH(fh) = true;
        end
    end

    % 2b. 形成新H₃ (形成阈值 r_H3_on)
    h2_keys = keys(h2_bonds);
    for k = 1:length(h2_keys)
        pair = h2_bonds(h2_keys{k});
        h1 = pair(1); h2 = pair(2);

        % 该H₂若已在活跃H₃中则跳过
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
            if bonded_atoms(at), continue; end      % 必须是游离H
            if used_freeH(at), continue; end        % 未被其他H₃占用
            d_min = min(pbc_dist(frac, h1, at, lat), pbc_dist(frac, h2, at, lat));
            if d_min <= r_H3_on
                s = sort([h1, h2, at]);
                key = sprintf('%d_%d_%d', s(1), s(2), s(3));
                new_h3(key) = [h1, h2, at];
                used_freeH(at) = true;
            end
        end
    end

    h3_active = new_h3;

    % ========== 3. 记录本帧活跃H₃ ==========
    h3_keys = keys(h3_active);
    h3_this = zeros(length(h3_keys), 6);
    for k = 1:length(h3_keys)
        roles = h3_active(h3_keys{k});
        s = sort(roles);
        % 前3列=规范key, 后3列=实际角色 [h1,h2,freeH]
        h3_this(k, :) = [s(1), s(2), s(3), roles(1), roles(2), roles(3)];
    end
    frame_h3{fi} = h3_this;
end

%% ==================== 提取稳定H₃事件 ====================
fprintf('提取稳定H₃事件 (≥%d帧)...\n', min_h3_frames);

% 用规范key (排序三元组) 分组, 同时记录每帧的实际H₂配对
triplet2entries = containers.Map('KeyType','char','ValueType','any');
for fi = 1:n_analyzed
    h3 = frame_h3{fi};
    if isempty(h3), continue; end
    for r = 1:size(h3, 1)
        key = sprintf('%d_%d_%d', h3(r,1), h3(r,2), h3(r,3));  % 规范key
        % entry: [frame, h1_actual, h2_actual, freeH_actual]
        entry = [fi, h3(r,4), h3(r,5), h3(r,6)];
        if ~isKey(triplet2entries, key)
            triplet2entries(key) = [];
        end
        tmp = triplet2entries(key);
        tmp(end+1, :) = entry;
        triplet2entries(key) = tmp;
    end
end

events = [];
all_keys = keys(triplet2entries);
for k = 1:length(all_keys)
    key = all_keys{k};
    entries = triplet2entries(key);           % [frame, h1, h2, fh] 每行一帧

    % 按帧号排序
    [~, idx] = sort(entries(:,1));
    entries = entries(idx, :);

    % 提取连续段
    frm = entries(:,1);
    seg_start = 1;
    for e = 2:length(frm)
        if frm(e) ~= frm(e-1) + 1
            seg_len = e - seg_start;
            if seg_len >= min_h3_frames
                seg_entries = entries(seg_start:e-1, :);
                roles = get_first_roles(seg_entries, frames, start_frame, lat, bond_partner_all);
                events(end+1, :) = [roles, frm(seg_start), frm(e-1), seg_len];
            end
            seg_start = e;
        end
    end
    seg_len = length(frm) - seg_start + 1;
    if seg_len >= min_h3_frames
        seg_entries = entries(seg_start:end, :);
        roles = get_first_roles(seg_entries, frames, start_frame, lat, bond_partner_all);
        events(end+1, :) = [roles, frm(seg_start), frm(end), seg_len];
    end
end

n_events = size(events, 1);
fprintf('  稳定H₃事件: %d\n', n_events);

if n_events == 0
    fprintf('未找到稳定H₃！\n');
    return;
end

%% ==================== H交换检测 ====================
%  判据: H₃结束后向前搜索post_window帧,
%        若游离H与H₂中某H形成新H₂键 (三对距离中最小值),
%        则判定为交换反应: H(1)-H(2) + H(f) → [H₃]‡ → H(1)-H(f) + H(2)
fprintf('\n===== H交换检测 (窗口=%d帧) =====\n', post_window);

exchange_flags = false(n_events, 1);
exch_info      = zeros(n_events, 2);  % [new_pair_code, exchange_frame]
                                       % 2=H1-fH, 3=H2-fH

for ev = 1:n_events
    h1 = events(ev,1);
    h2 = events(ev,2);
    fh = events(ev,3);
    fe_local = events(ev,5);          % 事件结束帧 (相对索引)

    for fd = fe_local+1 : min(fe_local+post_window, n_analyzed)
        frac = frames{start_frame + fd - 1} / lat;

        d12 = norm((frac(h1,:) - frac(h2,:) - ...
               round(frac(h1,:) - frac(h2,:))) * lat);
        d1f = norm((frac(h1,:) - frac(fh,:) - ...
               round(frac(h1,:) - frac(fh,:))) * lat);
        d2f = norm((frac(h2,:) - frac(fh,:) - ...
               round(frac(h2,:) - frac(fh,:))) * lat);

        [~, idx] = min([d12, d1f, d2f]);

        if idx == 1 && d12 <= r_HH_on
            % H₁-H₂恢复 → 游离H离开 → 无交换
            break;
        elseif (idx == 2 && d1f <= r_HH_on) || (idx == 3 && d2f <= r_HH_on)
            % 游离H与H₁或H₂成键 → 交换!
            exchange_flags(ev) = true;
            exch_info(ev, :) = [idx+1, fd];   % idx+1: 2=H1-fH, 3=H2-fH
            break;
        end
    end
end

n_exchange = sum(exchange_flags);
fprintf('  发生H交换: %d / %d (%.1f%%)\n', ...
    n_exchange, n_events, 100 * n_exchange / n_events);

if n_exchange > 0
    fprintf('\n  %4s %6s %6s %6s %6s %8s %10s %12s\n', ...
        '#','H1','H2','freeH','Dur','新配对','交换帧','新游离H');
    fprintf('  %s\n', repmat('-',1,62));
    cnt = 0;
    for ev = 1:n_events
        if exchange_flags(ev)
            cnt = cnt + 1;
            info = exch_info(ev, :);
            if info(1) == 2
                pair_str = 'H1-fH ';
                orphan   = events(ev,2);
            else
                pair_str = 'H2-fH ';
                orphan   = events(ev,1);
            end
            fprintf('  %4d %6d %6d %6d %6d %8s %10d %12d\n', ...
                cnt, events(ev,1), events(ev,2), events(ev,3), ...
                events(ev,6), pair_str, ...
                info(2) + start_frame - 1, orphan);
        end
    end
    fprintf('  %s\n\n', repmat('-',1,62));
end

if n_exchange == 0
    fprintf('\n未找到任何H交换事件！终止分析。\n');
    return;
end

% 过滤: 只保留交换事件
events = events(exchange_flags, :);
n_events = n_exchange;
fprintf('  后续分析: %d 个交换事件\n\n', n_events);

%% ==================== 收集H₃内部距离 + 键角 ====================
dist_all   = [];
dist_short = [];  % 每帧3个距离: 最短
dist_mid   = [];  % 每帧3个距离: 中间
dist_long  = [];  % 每帧3个距离: 最长
angles     = [];

ev_table = zeros(n_events, 10);
% 列: [H1, H2, freeH, dur, d12_mean, d13_mean, d23_mean, ang_mean, start_fr, end_fr]

for ev = 1:n_events
    h1_idx = events(ev,1);
    h2_idx = events(ev,2);
    fh_idx = events(ev,3);
    fs = events(ev,4);
    fe = events(ev,5);
    dur = events(ev,6);

    ev_d12 = []; ev_d13 = []; ev_d23 = []; ev_ang = [];

    for fd = fs:fe
        frac = frames{start_frame + fd - 1} / lat;

        % d12: H₁-H₂
        d12 = norm((frac(h1_idx,:) - frac(h2_idx,:) - ...
               round(frac(h1_idx,:) - frac(h2_idx,:))) * lat);

        % d13: freeH-H₁
        d13 = norm((frac(fh_idx,:) - frac(h1_idx,:) - ...
               round(frac(fh_idx,:) - frac(h1_idx,:))) * lat);

        % d23: freeH-H₂
        d23 = norm((frac(fh_idx,:) - frac(h2_idx,:) - ...
               round(frac(fh_idx,:) - frac(h2_idx,:))) * lat);

        dist_all = [dist_all, d12, d13, d23];
        ev_d12(end+1) = d12; ev_d13(end+1) = d13; ev_d23(end+1) = d23;

        % 按距离排序: 最短/中间/最长
        d_sorted = sort([d12, d13, d23]);
        dist_short(end+1) = d_sorted(1);
        dist_mid(end+1)   = d_sorted(2);
        dist_long(end+1)  = d_sorted(3);

        % 键角: 顶点 = 最短两条键的公共原子 (即与最长键相对的原子)
        if d12 >= d13 && d12 >= d23
            % d12(H₁-H₂) 最长 → 顶点 = 游离H, 角由 d13, d23 组成
            cos_ang = (d13^2 + d23^2 - d12^2) / (2 * d13 * d23);
        elseif d13 >= d12 && d13 >= d23
            % d13(freeH-H₁) 最长 → 顶点 = H₂, 角由 d12, d23 组成
            cos_ang = (d12^2 + d23^2 - d13^2) / (2 * d12 * d23);
        else
            % d23(freeH-H₂) 最长 → 顶点 = H₁, 角由 d12, d13 组成
            cos_ang = (d12^2 + d13^2 - d23^2) / (2 * d12 * d13);
        end
        cos_ang = min(max(cos_ang, -1), 1);
        ang = acosd(cos_ang);
        angles(end+1) = ang;
        ev_ang(end+1) = ang;
    end

    ev_table(ev,:) = [h1_idx, h2_idx, fh_idx, dur, ...
        mean(ev_d12), mean(ev_d13), mean(ev_d23), mean(ev_ang), ...
        fs + start_frame - 1, fe + start_frame - 1];
end

n_pairs    = length(dist_all);
n_triplets = length(dist_short);
n_angles   = length(angles);
fprintf('H-H距离样本: %d  (每组分类: %d)\n', n_pairs, n_triplets);
fprintf('键角样本: %d\n', n_angles);

%% ==================== PDF (全部距离) ====================
count_r = histcounts(dist_all, r_edges);
P_r = count_r / (n_pairs * dr);

%% ==================== PDF 分类 (短/中/长) ====================
count_short = histcounts(dist_short, r_edges);
count_mid   = histcounts(dist_mid,   r_edges);
count_long  = histcounts(dist_long,  r_edges);

P_short = count_short / (n_triplets * dr);
P_mid   = count_mid   / (n_triplets * dr);
P_long  = count_long  / (n_triplets * dr);

fprintf('\n--- 分类统计 ---\n');
fprintf('  最短键: %.4f ± %.4f Å\n', mean(dist_short), std(dist_short));
fprintf('  中间键: %.4f ± %.4f Å\n', mean(dist_mid),   std(dist_mid));
fprintf('  最长键: %.4f ± %.4f Å\n', mean(dist_long),  std(dist_long));

%% ==================== 键角分布 ====================
count_ang = histcounts(angles, ang_edges);
P_theta_raw = count_ang / (n_angles * dtheta);

sin_theta = sind(ang_centers);
sin_theta(sin_theta < 0.15) = NaN;
P_theta_corr = P_theta_raw ./ sin_theta;
valid = ~isnan(P_theta_corr);
P_theta_corr(valid) = P_theta_corr(valid) / sum(P_theta_corr(valid) * dtheta);

%% ==================== 导出数据 ====================
% 全部距离PDF
writematrix([r_centers(:), P_r(:)], 'H3_PDF.csv');
fprintf('\n数据已保存: H3_PDF.csv (r, P(r))\n');

% 分类PDF (短/中/长)
writematrix([r_centers(:), P_short(:), P_mid(:), P_long(:)], ...
    'H3_PDF_classified.csv');
fprintf('数据已保存: H3_PDF_classified.csv (r, P_short, P_mid, P_long)\n');

% 事件详情
ev_header = {'H1_idx','H2_idx','freeH_idx','duration', ...
    'd12_mean_A','d13_mean_A','d23_mean_A','angle_mean_deg', ...
    'start_frame','end_frame'};
ev_table_out = array2table(ev_table, 'VariableNames', ev_header);
writetable(ev_table_out, 'H3_events.csv');
fprintf('数据已保存: H3_events.csv (%d事件)\n', n_events);

% 键角
writematrix([ang_centers(:), P_theta_raw(:), P_theta_corr(:)], ...
    'H3_angle.csv');
fprintf('数据已保存: H3_angle.csv (θ, P_raw, P_corr)\n');

%% ==================== 绘图 ====================
figure('Position', [50 100 1300 850]);

c_short = [0.85 0.20 0.20];
c_mid   = [0.20 0.60 0.25];
c_long  = [0.20 0.30 0.75];

% (a) 全部距离PDF
subplot(2,2,1);
plot(r_centers, P_r, 'k-', 'LineWidth', 2);
xlabel('r (Å)', 'FontSize', 14);
ylabel('P(r) (Å^{-1})', 'FontSize', 14);
title(sprintf('H₃内部H-H距离PDF (全部 %d样本)', n_pairs), 'FontSize', 13);
grid on; xlim([0 r_max]);
set(gca, 'FontSize', 12);

% (b) 分类PDF: 最短/中间/最长
subplot(2,2,2);
hold on;
plot(r_centers, P_short, '-', 'Color', c_short, 'LineWidth', 2.5);
plot(r_centers, P_mid,   '-', 'Color', c_mid,   'LineWidth', 2.5);
plot(r_centers, P_long,  '-', 'Color', c_long,  'LineWidth', 2.5);
hold off;
xlabel('r (Å)', 'FontSize', 14);
ylabel('P(r) (Å^{-1})', 'FontSize', 14);
title(sprintf('H₃距离分类PDF (每组 %d样本)', n_triplets), 'FontSize', 13);
legend({sprintf('最短 <d>=%.3fÅ', mean(dist_short)), ...
        sprintf('中间 <d>=%.3fÅ', mean(dist_mid)), ...
        sprintf('最长 <d>=%.3fÅ', mean(dist_long))}, ...
    'FontSize', 11, 'Location', 'northeast');
grid on; xlim([0 r_max]);
set(gca, 'FontSize', 12);

% (c) 键角分布
subplot(2,2,3);
yyaxis left;
bar(ang_centers, P_theta_raw, 1.0, 'FaceColor', [0.6 0.7 0.9], ...
    'EdgeColor', 'none');
ylabel('P_{raw}(θ) (°^{-1})', 'FontSize', 14);

yyaxis right;
plot(ang_centers, P_theta_corr, 'ro-', 'LineWidth', 2, ...
    'MarkerSize', 6, 'MarkerFaceColor', [0.9 0.3 0.3]);
ylabel('P_{corr}(θ) = P_{raw}/sinθ (°^{-1})', 'FontSize', 14);
xlabel('键角 θ (°)', 'FontSize', 14);
title(sprintf('H₃键角分布 (%d样本)', n_angles), 'FontSize', 13);
grid on; xlim([0 180]);
legend({'P_{raw}(θ)', 'P_{raw}/sinθ (内禀偏好)'}, ...
    'FontSize', 10, 'Location', 'northwest');
set(gca, 'FontSize', 12);

% (d) 统计汇总
subplot(2,2,4);
axis off;
hold on;

txt_lines = {
    sprintf('===== H3 分析统计 ====='), ...
    sprintf(' '), ...
    sprintf('事件数: %d (>=%d帧, 含H交换)', n_events, min_h3_frames), ...
    sprintf('交换检测窗口: %d 帧', post_window), ...
    sprintf('总距离样本: %d', n_pairs), ...
    sprintf('每组分类样本: %d', n_triplets), ...
    sprintf(' '), ...
    sprintf('--- 距离分类统计 ---'), ...
    sprintf('最短键: %.4f +/- %.4f A', mean(dist_short), std(dist_short)), ...
    sprintf('中间键: %.4f +/- %.4f A', mean(dist_mid),   std(dist_mid)), ...
    sprintf('最长键: %.4f +/- %.4f A', mean(dist_long),  std(dist_long)), ...
    sprintf(' '), ...
    sprintf('--- 峰位 (分类PDF) ---'), ...
    sprintf('最短键峰: ~%.2f A', r_centers(find(P_short == max(P_short), 1))), ...
    sprintf('中间键峰: ~%.2f A', r_centers(find(P_mid   == max(P_mid),   1))), ...
    sprintf('最长键峰: ~%.2f A', r_centers(find(P_long  == max(P_long),  1))), ...
};

for i = 1:length(txt_lines)
    text(0.05, 1 - i*0.065, txt_lines{i}, 'FontSize', 11, ...
        'FontName', 'Consolas', 'VerticalAlignment', 'top');
end
hold off;

sgtitle(sprintf('H3分子内部分析 (H交换事件, %d个, >=%d帧)', ...
    n_events, min_h3_frames), 'FontSize', 15, 'FontWeight', 'bold');

saveas(gcf, 'H3_PDF_angle.png');
fprintf('图表已保存: H3_PDF_angle.png\n');

fprintf('\n完成。\n');

%% ==================== 辅助函数 ====================
function roles = get_first_roles(seg_entries, frames, start_frame, lat, bond_partner_all)
% 确定H₂配对: 优先回溯上一帧的H₂键, 其次退回首帧键长排序
%   seg_entries: [frame, a, b, c] (三个原子, 已按帧号排序)
%   返回: [H1, H2, freeH] — H1-H2为H₂, freeH为外来H

    first_frame = seg_entries(1, 1);
    a = seg_entries(1, 2); b = seg_entries(1, 3); c = seg_entries(1, 4);

    % 回溯上一帧的成键关系
    h1 = 0; h2 = 0; fh = 0;
    if first_frame > 1
        bp = bond_partner_all(:, first_frame - 1);
        if bp(a) == b
            h1 = a; h2 = b; fh = c;
        elseif bp(a) == c
            h1 = a; h2 = c; fh = b;
        elseif bp(b) == c
            h1 = b; h2 = c; fh = a;
        end
    end
    if h1 == 0
        % 回溯失败 (首帧 / 上一帧这三原子互不成键) → 退回首帧键长排序
        frac = frames{start_frame + first_frame - 1} / lat;
        dab = pbc_dist(frac, a, b, lat);
        dac = pbc_dist(frac, a, c, lat);
        dbc = pbc_dist(frac, b, c, lat);
        if dab <= dac && dab <= dbc
            h1 = a; h2 = b; fh = c;
        elseif dac <= dab && dac <= dbc
            h1 = a; h2 = c; fh = b;
        else
            h1 = b; h2 = c; fh = a;
        end
    end
    roles = [h1, h2, fh];
end

function d = pbc_dist(frac, i, j, lat)
% 最小镜像距离: 原子i与j在分数坐标frac下的PBC距离
    v = frac(i,:) - frac(j,:);
    v = v - round(v);
    d = norm(v * lat);
end
