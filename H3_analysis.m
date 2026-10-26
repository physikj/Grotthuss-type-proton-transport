%% H2-H2O System -- Spatial (rho, zeta) Distribution of H in Stable H3 Complexes
%  Reads VASP XDATCAR -> identifies H2 molecules -> finds free H -> filters
%  H3 persisting >= N frames -> exchange reaction filter -> 2D distribution plot.
%
%  Coordinate system (H2 molecule as reference):
%    zeta : along H2 bond axis, origin at bond centre
%    rho  : perpendicular distance to the bond axis
%
%  Exchange filter: only counts H3 events where dissociation leaves with
%    one of the original H2 atoms (not the incoming free H).

clear; close all; clc;

%% ======================== User parameters ========================
fpath = 'XDATCAR_60_2500_3';          % VASP trajectory file

% --- H-H bond criteria (Angstrom, 带滞回) ---
r_HH_min = 0.50;                 % H₂键长下限 (排除噪声)
r_HH_on  = 1.00;                 % H₂形成阈值
r_HH_off = 1.20;                 % H₂断裂阈值

% --- H3 identification (带滞回) ---
r_H3_on  = 1.00;                 % H₃形成阈值: 游离H < 1.0
r_H3_off = 1.20;                 % H₃分离阈值: 游离H > 1.2

% --- Stability requirement ---
min_h3_frames = 1;               % minimum consecutive frames for "stable" H3

% --- Exchange detection window ---
post_window = 5;                 % frames after event end to search

% --- Grid ---
npoint   = 200;                  % grid points per dimension
rho_max  = 3;                  % rho axis upper bound (Angstrom)
zeta_max = 1.5;                  % |zeta| axis upper bound (Angstrom)

% --- Analysis window ---
n_frames = 20000;                 % take last N frames

% --- Display options ---
plot_lifespan_hist = true;       % whether to plot H3 lifetime histogram

%% ======================== Publication-quality figure setup ========================
set(0, 'DefaultAxesFontName',   'Arial');
set(0, 'DefaultTextFontName',   'Arial');
set(0, 'DefaultColorbarFontName','Arial');
set(0, 'DefaultAxesFontSize',   12);
set(0, 'DefaultTextFontSize',   12);
set(0, 'DefaultLineLineWidth',  1.2);

colormap_choice = 'parula';
export_dpi = 300;

%% ======================== Step 1: Read XDATCAR ========================
fprintf('========================================\n');
fprintf('  H2-H2O  H3 Spatial Distribution\n');
fprintf('  (+ exchange reaction filter)\n');
fprintf('========================================\n\n');
fprintf('[1/6] Reading XDATCAR: %s\n', fpath);

fp = fopen(fpath, 'r');
assert(fp ~= -1, 'ERROR: cannot open file: %s', fpath);

fgetl(fp);
scale = fscanf(fp, '%f', 1);
fgetl(fp);

lat = zeros(3, 3);
for i = 1:3
    lat(i, :) = sscanf(fgetl(fp), '%f %f %f')';
end

fgetl(fp);
counts = sscanf(fgetl(fp), '%d %d');
n_H = counts(1);
n_O = counts(2);
n_atom = n_H + n_O;
fprintf('    Atoms: H = %d, O = %d\n', n_H, n_O);

tline = fgetl(fp);
iter = 0;
while ischar(tline) && ~contains(tline, 'Direct') && ~contains(tline, 'Cartesian')
    tline = fgetl(fp);
    iter = iter + 1;
    assert(iter < 200, 'ERROR: coordinate header not found');
end

%% ======================== Step 2: Read all frames ========================
fprintf('[2/6] Loading trajectory coordinates ...\n');

frames = {};
while ischar(tline)
    if contains(tline, 'Direct') || contains(tline, 'Cartesian')
        crd = zeros(n_atom, 3);
        for i = 1:n_atom
            tline = fgetl(fp);
            if ~ischar(tline), break; end
            crd(i, :) = sscanf(tline, '%f %f %f')';
        end
        frames{end+1} = crd * lat;
    end
    tline = fgetl(fp);
end
fclose(fp);

n_total = length(frames);
start_frame = max(1, n_total - n_frames + 1);
n_analyzed  = n_total - start_frame + 1;
fprintf('    Analysis window: frame %d -> %d (%d frames)\n', ...
        start_frame, n_total, n_analyzed);

%% ======================== Step 3: Per-frame H2 & H3 detection (带滞回) ========================
fprintf('[3/6] Detecting H2/H3 (滞回: H2 %g->%g, H3 %g->%g) ...\n', ...
    r_HH_on, r_HH_off, r_H3_on, r_H3_off);

frame_h3  = cell(n_analyzed, 1);
bond_lens = [];
bond_partner_all = zeros(n_H, n_analyzed);   % 每帧成键查找表 (供回溯上一帧定H₂身份)

% 跨帧持久状态 (滞回)
h2_bonds  = containers.Map('KeyType','char','ValueType','any');  % 当前H₂键
h3_active = containers.Map('KeyType','char','ValueType','any');  % 当前活跃H₃

for fi = 1:n_analyzed
    if mod(fi, 100) == 0
        fprintf('    Frame %4d / %d\n', fi, n_analyzed);
    end

    frac = frames{start_frame + fi - 1} / lat;

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

    % 收集当前所有H₂键长 (供统计)
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
    bond_partner_all(:, fi) = bond_partner;   % 存下本帧成键关系

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
            continue;                     % 无成键对 → H₃已解散
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
            if bonded_atoms(at), continue; end      % 必须是游离H
            if used_freeH(at), continue; end        % 未被其他H₃占用
            d_min = min(pbc_dist(frac, h1, at, lat), pbc_dist(frac, h2, at, lat));
            if d_min <= r_H3_on
                s = sort([h1, h2, at]);
                new_h3(sprintf('%d_%d_%d', s(1), s(2), s(3))) = [h1, h2, at];
                used_freeH(at) = true;
            end
        end
    end

    h3_active = new_h3;

    % ---- 3. 记录本帧H₃ (规范key + 角色) ----
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

fprintf('    Done. %d H2 bonds identified.\n', length(bond_lens));

%% ======================== Step 4: Extract stable H3 events ========================
fprintf('[4/6] Extracting stable H3 events (>= %d frames) ...\n', min_h3_frames);

% 用规范key (排序三元组) 分组, 记录每帧角色
triplet2entries = containers.Map('KeyType', 'char', 'ValueType', 'any');
for fi = 1:n_analyzed
    h3 = frame_h3{fi};
    if isempty(h3), continue; end
    for r = 1:size(h3, 1)
        key = sprintf('%d_%d_%d', h3(r,1), h3(r,2), h3(r,3));  % 规范key
        entry = [fi, h3(r,4), h3(r,5), h3(r,6)];  % [frame, h1, h2, fh]
        if ~isKey(triplet2entries, key)
            triplet2entries(key) = [];
        end
        tmp = triplet2entries(key);
        tmp(end+1, :) = entry;
        triplet2entries(key) = tmp;
    end
end

% events: [h1, h2, fh, zeta, rho, dur, start, end, bond_len]
events = [];
all_keys = keys(triplet2entries);
for k = 1:length(all_keys)
    key = all_keys{k};
    entries = triplet2entries(key);
    [~, idx] = sort(entries(:,1));
    entries = entries(idx, :);

    frm = entries(:,1);
    seg_start = 1;
    for e = 2:length(frm)
        if frm(e) ~= frm(e-1) + 1
            seg_len = e - seg_start;
            if seg_len >= min_h3_frames
                events(end+1,:) = build_event(entries(seg_start:e-1,:), ...
                    frames, start_frame, lat, seg_len, frm(seg_start), frm(e-1), ...
                    bond_partner_all);
            end
            seg_start = e;
        end
    end
    seg_len = length(frm) - seg_start + 1;
    if seg_len >= min_h3_frames
        events(end+1,:) = build_event(entries(seg_start:end,:), ...
            frames, start_frame, lat, seg_len, frm(seg_start), frm(end), ...
            bond_partner_all);
    end
end

%% ======================== Step 5: Exchange reaction filter ========================
fprintf('[5/6] Filtering exchange reactions (窗口=%d帧) ...\n', post_window);

exchange_idx = false(size(events,1), 1);

for ev = 1:size(events,1)
    h1 = events(ev, 1);   % H₂原子1 (键长排序后)
    h2 = events(ev, 2);   % H₂原子2
    fh = events(ev, 3);   % 外来H
    fe_local = events(ev, 8);

    for fd = fe_local+1 : min(fe_local+post_window, n_analyzed)
        frac = frames{start_frame + fd - 1} / lat;
        d12 = pbc_dist(frac, h1, h2, lat);
        d1f = pbc_dist(frac, h1, fh, lat);
        d2f = pbc_dist(frac, h2, fh, lat);
        [~, idx] = min([d12, d1f, d2f]);

        if idx == 1 && d12 <= r_HH_on
            % H1-H2恢复 -> 外来H离开 -> 无交换
            break;
        elseif (idx == 2 && d1f <= r_HH_on) || (idx == 3 && d2f <= r_HH_on)
            % 外来H-H1 或 外来H-H2成键 -> 交换! (原H₂原子被置换)
            exchange_idx(ev) = true;
            break;
        end
    end
end

n_total_events = size(events, 1);
n_exchange = sum(exchange_idx);
fprintf('    Exchange reactions: %d / %d (%.1f%%)\n', ...
    n_exchange, n_total_events, 100*n_exchange/max(1,n_total_events));

events = events(exchange_idx, :);
n_events = size(events, 1);

%% ======================== Check results ========================
if n_events == 0
    fprintf('\n*** No exchange H3 events found. ***\n');
    fprintf('Try: increase r_H3_on, decrease min_h3_frames\n');
    return;
end

fprintf('    Stable H3 + exchange events: %d\n', n_events);
fprintf('    Lifetime: min=%d, max=%d, mean=%.1f frames\n', ...
    min(events(:,6)), max(events(:,6)), mean(events(:,6)));

% --- Print event table ---
fprintf('\n');
fprintf('    ===== Exchange H3 Event List (>= %d frames, %d total) =====\n', ...
    min_h3_frames, n_events);
fprintf('    %4s %8s %8s %8s %6s %10s %10s %12s\n', ...
    '#','H2_H1','H2_H2','FreeH','Dur','Start','End','Abs.Frames');
fprintf('    %s\n', repmat('-',1,72));
for e = 1:n_events
    fprintf('    %4d %8d %8d %8d %6d %10d %10d %8d-%d\n', ...
        e, events(e,1), events(e,2), events(e,3), events(e,6), ...
        events(e,7), events(e,8), ...
        events(e,7)+start_frame-1, events(e,8)+start_frame-1);
end
fprintf('    %s\n\n', repmat('-',1,72));

% --- Duration distribution ---
dur = events(:,6);
dur_vals = unique(dur);
fprintf('    Lifetime distribution:\n');
for d = 1:length(dur_vals)
    n_d = sum(dur == dur_vals(d));
    fprintf('      %2d frames: %4d events  (%5.1f%%)\n', ...
        dur_vals(d), n_d, 100*n_d/n_events);
end

mismatch = find((events(:,8) - events(:,7) + 1) ~= events(:,6));
if ~isempty(mismatch)
    fprintf('\n    WARNING: %d events have frame/duration mismatch!\n', ...
        length(mismatch));
else
    fprintf('\n    Internal consistency: PASSED\n');
end

%% ======================== Step 6: 2D distribution & figures ========================
fprintf('[6/6] Plotting (rho, zeta) distribution ...\n');

grid_rho  = linspace(0, rho_max, npoint);
grid_zeta = linspace(-zeta_max, zeta_max, npoint);

% --- 分箱: 逐帧统计游离H 与 H₂原子 (每帧取最短键 = H₂, 第三个 = 游离H) ---
n_samples = sum(events(:, 6));        % 总事件-帧数 (每帧一个游离H样本)
count_free = zeros(npoint, npoint);   % 游离H计数
zeta_h2 = [];                          % H₂原子的 ζ 坐标 (ρ=0, ζ=±bl/2)

for e = 1:n_events
    a = events(e, 1);                  % 三个原子 (任意顺序)
    b = events(e, 2);
    c = events(e, 3);
    fs = events(e, 7);                 % 事件起止帧 (相对索引)
    fe = events(e, 8);

    for fd = fs:fe
        frac = frames{start_frame + fd - 1} / lat;

        % 每帧取最短键: 三原子中距离最短的两原子 = H₂, 第三个 = 游离H
        [h1, h2, fh] = get_h3_roles(frac, a, b, c, lat);

        % 键轴 + 质心 (瞬时 H₂)
        d_frac = frac(h2,:) - frac(h1,:);
        d_frac = d_frac - round(d_frac);
        bl = norm(d_frac * lat);                 % 该帧键长
        cent_frac = frac(h1,:) + d_frac/2;
        ax = (d_frac * lat) / bl;
        if rand() < 0.5, ax = -ax; end           % H₂同核, ±ζ对称化

        % 游离H相对质心的 (ζ, ρ)
        d_c = frac(fh,:) - cent_frac;
        d_c = d_c - round(d_c);
        d_pbc = d_c * lat;
        zeta = dot(d_pbc, ax);
        rho  = sqrt(max(0, norm(d_pbc)^2 - zeta^2));

        % 游离H分箱
        [~, ir] = min(abs(grid_rho  - rho));
        [~, iz] = min(abs(grid_zeta - zeta));
        count_free(ir, iz) = count_free(ir, iz) + 1;

        % H₂两个原子标记
        zeta_h2(end+1) = -bl/2;   % H₁
        zeta_h2(end+1) = +bl/2;   % H₂
    end
end

% --- 分箱体积归一化: 均匀网格面积 × 柱坐标雅可比 2πρ ---
drho  = rho_max  / (npoint - 1);           % ρ 方向格点间距 (Å)
dzeta = 2*zeta_max / (npoint - 1);         % ζ 方向格点间距 (Å)
bin_area = drho * dzeta;                   % 单个格点面积 (Å²)

rho_col = grid_rho(:);                      % npoint×1, 行对应 ρ
jac = repmat(2*pi*rho_col, 1, npoint);      % 体积元 = 2πρ dρ dζ 中的 2πρ
jac(jac == 0) = NaN;                        % ρ=0 处奇异, 屏蔽

% --- 增强平滑: 先平滑计数, 再除体积 (避免把雅可比奇异值抹开) ---
sigma_g = 2.0;                              % 原来0.8, 加大去散点
kernel_size = ceil(6 * sigma_g);            % ~13, 核至少6σ宽
if mod(kernel_size, 2) == 0
    kernel_size = kernel_size + 1;          % 保证奇数
end
gf = fspecial('gaussian', [kernel_size kernel_size], sigma_g);
count_free_smooth = imfilter(count_free, gf, 'replicate');

% 游离H三维概率密度 (Å^-3): P(ρ,ζ) = count / (N_frames · ΔρΔζ · 2πρ), ∫P dV = 1
prob_H_smooth = count_free_smooth ./ (n_samples * bin_area .* jac);
prob_H_smooth(isnan(prob_H_smooth)) = 0;    % ρ=0 处无意义, 置 0

% ---- FIGURE 1: Main (rho, zeta) distribution ----
fig1 = figure('Position', [100, 100, 760, 640], ...
              'Color', 'w', 'Renderer', 'painters');

imagesc(grid_rho, grid_zeta, prob_H_smooth');
set(gca, 'YDir', 'normal');
shading interp;                             % 平滑渲染, 消除像素感
colormap(gca, colormap_choice);

cb = colorbar;
cb.Label.String = 'Probability density (Å^{-3})';
cb.Label.FontSize = 20;
cb.FontSize = 20;
cb.Box = 'on';
cb.LineWidth = 1.0;

avg_bond = mean(bond_lens);

% --- H₂ 原子: 黄色平滑小球 (高斯涂抹) ---
count_h2 = zeros(npoint, npoint);
[~, ir0] = min(abs(grid_rho - 0));
for k = 1:length(zeta_h2)
    [~, iz] = min(abs(grid_zeta - zeta_h2(k)));
    count_h2(ir0, iz) = count_h2(ir0, iz) + 1;
end
count_h2_smooth = imfilter(count_h2, gf, 'replicate');

hold on;
contourf(grid_rho, grid_zeta, count_h2_smooth', ...
    [0.15*max(count_h2_smooth(:)), max(count_h2_smooth(:))], ...
    'FaceColor', 'y', 'EdgeColor', 'none');
hold off;

xlabel('\rho  (Å)',  'FontSize', 23);
ylabel('\zeta  (Å)', 'FontSize', 23);
title(sprintf(['Spatial distribution of free H around H_2\n' ...
       '(exchange H3, >= %d frames, %d events, <d_{HH}> = %.3f %s)'], ...
       min_h3_frames, n_events, avg_bond, char(197)), ...
       'FontSize', 14, 'FontWeight', 'bold');

% --- 颜色范围: 压低上限让中心过曝、低强度区域可见 ---
clim([0, max(prob_H_smooth(:)) * 1]);

% --- 坐标范围: ρ轴左移给H₂标记留空间 ---
xlim([0, rho_max]);
ylim([-zeta_max, zeta_max]);
box on;
set(gca, 'FontSize', 20, 'LineWidth', 1.2);
set(gca, 'TickDir', 'out', 'TickLength', [0.005 0.005]);
% ★ 设置 x 轴刻度：从 0 到 rho_max，步长 0.5
xticks(0:0.5:rho_max);
exportgraphics(fig1, 'H3_spatial_distribution.png', ...
    'Resolution', export_dpi, 'BackgroundColor', 'white');
exportgraphics(fig1, 'H3_spatial_distribution.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'white');
fprintf('    Figures saved: H3_spatial_distribution.png / .pdf\n');

% ---- FIGURE 2: Lifetime histogram ----
if plot_lifespan_hist
    fig2 = figure('Position', [880, 100, 420, 360], ...
                  'Color', 'w', 'Renderer', 'painters');
    histogram(events(:,6), 'BinMethod', 'integers', ...
              'FaceColor', [0.25 0.55 0.85], 'EdgeColor', 'k', ...
              'LineWidth', 1.0, 'FaceAlpha', 0.85);
    xlabel('H_3 lifetime (frames)', 'FontSize', 15, 'FontWeight', 'bold');
    ylabel('Number of events', 'FontSize', 15, 'FontWeight', 'bold');
    title(sprintf('H_3 lifetime (%d exchange events)', n_events), ...
          'FontSize', 13, 'FontWeight', 'bold');
    grid on; box on;
    set(gca, 'FontSize', 13, 'LineWidth', 1.1);
    set(gca, 'TickDir', 'out', 'TickLength', [0.005 0.005]);

    exportgraphics(fig2, 'H3_lifespan_histogram.png', ...
        'Resolution', export_dpi, 'BackgroundColor', 'white');
    exportgraphics(fig2, 'H3_lifespan_histogram.pdf', ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    fprintf('    Figures saved: H3_lifespan_histogram.png / .pdf\n');
end

%% ======================== Summary ========================
fprintf('\n');
fprintf('========================================\n');
fprintf('           Analysis complete\n');
fprintf('========================================\n');
fprintf('  Frames analysed        : %d\n', n_analyzed);
fprintf('  Mean H2 bond length    : %.4f %s\n', avg_bond, char(197));
fprintf('  H3 distance threshold  : <= %.2f %s (form) / %.2f (break)\n', ...
    r_H3_on, char(197), r_H3_off);
fprintf('  Stability requirement  : >= %d consecutive frames\n', min_h3_frames);
fprintf('  Total stable H3 events : %d\n', n_total_events);
fprintf('  Exchange H3 events     : %d (%.1f%%)\n', n_exchange, ...
    100*n_exchange/max(1,n_total_events));
fprintf('  Mean event lifetime    : %.1f frames\n', mean(events(:,6)));
fprintf('  Colormap               : %s\n', colormap_choice);
fprintf('  Output formats         : PNG (%d dpi) + PDF (vector)\n', export_dpi);
fprintf('========================================\n');

%% ======================== 辅助函数 ========================
function ev = build_event(seg_entries, frames, start_frame, lat, dur, fstart, fend, bond_partner_all)
% 确定事件角色: 优先回溯上一帧的H₂键, 其次退回首帧键长排序
%   seg_entries: [frame, a, b, c] (已按帧号排序)
%   返回 ev: [h1, h2, fh, zeta, rho, dur, fstart, fend, bond_len]
%   h1-h2 为原本的H₂(回溯上一帧定), fh 为外来H; 首帧定死, 整个事件不变

    fs = seg_entries(1, 1);
    a = seg_entries(1, 2); b = seg_entries(1, 3); c = seg_entries(1, 4);

    % --- 确定 H₂ 身份: 回溯上一帧的成键关系 ---
    h1 = 0; h2 = 0; fh = 0;
    if fs > 1
        bp = bond_partner_all(:, fs - 1);   % 上一帧各H的成键对象
        if bp(a) == b
            h1 = a; h2 = b; fh = c;         % 上一帧 a-b 成键 → 原本H₂ = a-b
        elseif bp(a) == c
            h1 = a; h2 = c; fh = b;
        elseif bp(b) == c
            h1 = b; h2 = c; fh = a;
        end
    end
    if h1 == 0
        % 回溯失败 (首帧 / 上一帧这三原子互不成键) → 退回首帧键长排序
        frac0 = frames{start_frame + fs - 1} / lat;
        dab = pbc_dist(frac0, a, b, lat);
        dac = pbc_dist(frac0, a, c, lat);
        dbc = pbc_dist(frac0, b, c, lat);
        if dab <= dac && dab <= dbc
            h1 = a; h2 = b; fh = c;
        elseif dac <= dab && dac <= dbc
            h1 = a; h2 = c; fh = b;
        else
            h1 = b; h2 = c; fh = a;
        end
    end

    % H₂质心 + 键轴 (分数坐标, 用首帧坐标)
    frac0 = frames{start_frame + fs - 1} / lat;
    d_frac = frac0(h2,:) - frac0(h1,:);
    d_frac = d_frac - round(d_frac);
    bond_len = norm(d_frac * lat);
    cent_frac = frac0(h1,:) + d_frac/2;
    ax = (d_frac * lat) / bond_len;

    % 随机翻转键轴 (H₂同核, ±ζ等价)
    if rand() < 0.5, ax = -ax; end

    % 游离H相对质心
    d_frac_c = frac0(fh,:) - cent_frac;
    d_frac_c = d_frac_c - round(d_frac_c);
    d_pbc = d_frac_c * lat;
    zeta = dot(d_pbc, ax);
    dist2c = norm(d_pbc);
    rho = sqrt(max(0, dist2c^2 - zeta^2));

    ev = [h1, h2, fh, zeta, rho, dur, fstart, fend, bond_len];
end

function [h1, h2, fh] = get_h3_roles(frac, a, b, c, lat)
% 三原子中键长最短的两原子 = H₂, 第三个 = 游离H
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

function d = pbc_dist(frac, i, j, lat)
% 最小镜像距离: 原子i与j在分数坐标frac下的PBC距离
    v = frac(i,:) - frac(j,:);
    v = v - round(v);
    d = norm(v * lat);
end
