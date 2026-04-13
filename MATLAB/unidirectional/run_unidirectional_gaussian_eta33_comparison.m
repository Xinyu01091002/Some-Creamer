% Plot finite-depth Gaussian narrow-band eta33 comparison:
%   MF12, Creamer(H3), Creamer(H3+H4).

clearvars -except Akp bandwidth_ratio Nx x_vec h kp max_active_modes max_triad_active_modes mf12_max_active_modes triad_backend n_lambda_steps;
clc;

this_dir = fileparts(mfilename('fullpath'));
matlab_dir = fileparts(this_dir);
repo_root = fileparts(matlab_dir);
addpath(this_dir);
addpath(fullfile(matlab_dir, 'core'));
addpath(fullfile(matlab_dir, 'validation'));

if ~exist('Akp', 'var') || isempty(Akp)
    Akp = 0.05;
end
if ~exist('bandwidth_ratio', 'var') || isempty(bandwidth_ratio)
    bandwidth_ratio = 0.05;
end
if ~exist('kp', 'var') || isempty(kp)
    kp = 0.0279;
end
if ~exist('Nx', 'var') || isempty(Nx)
    Nx = 2048;
end
if ~exist('h', 'var') || isempty(h)
    h = 1 / kp;
end
if ~exist('max_triad_active_modes', 'var') || isempty(max_triad_active_modes)
    max_triad_active_modes = 32;
end
if ~exist('max_active_modes', 'var') || isempty(max_active_modes)
    max_active_modes = 256;
end
if ~exist('mf12_max_active_modes', 'var') || isempty(mf12_max_active_modes)
    mf12_max_active_modes = 512;
end
if ~exist('n_lambda_steps', 'var') || isempty(n_lambda_steps)
    n_lambda_steps = 12;
end
if ~exist('triad_backend', 'var') || isempty(triad_backend)
    triad_backend = 'cpp';
end

g = 9.81;
lambda_p = 2*pi/kp;
if ~exist('x_vec', 'var') || isempty(x_vec)
    x_vec = linspace(-64 * lambda_p, 64 * lambda_p, Nx);
else
    x_vec = x_vec(:).';
    Nx = numel(x_vec);
end
dx = mean(diff(x_vec));
points_per_wavelength = lambda_p / dx;
domain_wavelengths = dx * Nx / lambda_p;

eta_lin = linear_focus_gaussian_spectrum_1d(Akp, bandwidth_ratio, x_vec, struct('kp', kp));

mf12_cfg = struct('g', g, 'h', h, 't', 0, ...
    'energy_fraction', 1, 'max_active_modes', mf12_max_active_modes, 'max_order', 3);
t_mf12 = tic;
mf12 = mf12_from_linear_focus_1d(eta_lin, x_vec, mf12_cfg);
mf12_runtime_s = toc(t_mf12);

sep_cfg = struct();
sep_cfg.g = g;
sep_cfg.depth_h = h;
sep_cfg.energy_fraction = 1;
sep_cfg.min_active_modes = 0;
sep_cfg.max_active_modes = max_active_modes;
sep_cfg.n_lambda_steps = n_lambda_steps;
sep_cfg.lambda_flow_model = 'canonical_pair';
sep_cfg.creamer_backend = 'matlab';
sep_cfg.propagation_direction_deg = 0;
sep_cfg.preserve_mean = true;
sep_cfg.verbose = false;
t_sep = tic;
creamer_sep = creamer_four_phase_separation_1d(eta_lin, x_vec, sep_cfg);
sep_runtime_s = toc(t_sep);

triad_cfg = struct('g', g, 'max_triad_active_modes', max_triad_active_modes);
t_creamer = tic;
if strcmpi(triad_backend, 'cpp')
    triad_cfg.cpp_exe = fullfile(repo_root, 'cpp', 'creamer_flow', 'build', 'creamer_eta33_triad_1d.exe');
    triad_cfg.cpp_job_dir = fullfile(tempdir, 'creamer_cpp_gaussian_eta33_triad');
    creamer = creamer_eta33_h3h4_triad_1d_cpp(eta_lin, x_vec, h, triad_cfg);
else
    creamer = creamer_eta33_h3h4_triad_1d(eta_lin, x_vec, h, triad_cfg);
end
creamer_runtime_s = toc(t_creamer);

mx = [0:(Nx/2), (-Nx/2 + 1):-1];
k = (2*pi/(dx*Nx)) * mx;
k_over_kp = k / kp;
pos_mask = mx >= 0;
spec_mask = (k_over_kp(pos_mask) >= 0) & (k_over_kp(pos_mask) <= 6);
x_over_lambda = x_vec / lambda_p;
xlim_wave = [-2, 2];

spec = struct();
spec.mf12_eta22 = abs(fft(mf12.eta22) / Nx);
spec.creamer_eta22 = abs(fft(creamer_sep.eta22) / Nx);
spec.mf12_eta33 = abs(fft(mf12.eta33) / Nx);
spec.creamer_h3_eta33 = abs(creamer.eta33_h3_hat);
spec.creamer_h3h4_eta33 = abs(creamer.eta33_h3h4_hat);
spec.eta_lin = abs(fft(eta_lin) / Nx);
env = struct();
env.mf12_eta22 = abs(hilbert(mf12.eta22));
env.creamer_eta22 = abs(hilbert(creamer_sep.eta22));
env.mf12_eta33 = abs(hilbert(mf12.eta33));
env.creamer_h3_eta33 = abs(hilbert(creamer.eta33_h3));
env.creamer_h3h4_eta33 = abs(hilbert(creamer.eta33_h3h4));

[~, idx_3kp_full] = min(abs(k_over_kp - 3));
mf12_3kp = spec.mf12_eta33(idx_3kp_full);
creamer_h3_3kp = spec.creamer_h3_eta33(idx_3kp_full);
creamer_h3h4_3kp = spec.creamer_h3h4_eta33(idx_3kp_full);
kpos = k_over_kp(pos_mask);
triad_active_count = local_triad_active_count(creamer);

out_dir = fullfile(repo_root, 'MATLAB', 'output', 'unidirectional');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

fig = figure('Color', 'w', 'Position', [120 120 1480 920]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(x_over_lambda, mf12.eta22, '-', 'LineWidth', 2.2, 'Color', [0.00 0.45 0.74]);
hold on;
plot(x_over_lambda, creamer_sep.eta22, '--', 'LineWidth', 1.45, 'Color', [0.85 0.33 0.10]);
plot(x_over_lambda, creamer_sep.eta22, ':', 'LineWidth', 2.0, 'Color', [0.47 0.67 0.19]);
plot(x_over_lambda, env.mf12_eta22, '-', 'LineWidth', 1.9, 'Color', [0.45 0.70 0.90]);
plot(x_over_lambda, env.creamer_eta22, '--', 'LineWidth', 1.7, 'Color', [0.95 0.62 0.40]);
plot(x_over_lambda, env.creamer_eta22, ':', 'LineWidth', 2.0, 'Color', [0.70 0.84 0.42]);
grid on;
xlim(xlim_wave);
xlabel('x / \lambda_p');
ylabel('\eta_{22}');
title('\eta_{22}');
legend({'MF12', 'Creamer(H3)', 'Creamer(H3+H4)', ...
    'MF12 envelope', 'Creamer(H3) envelope', 'Creamer(H3+H4) envelope'}, ...
    'Location', 'best');

nexttile;
plot(x_over_lambda, creamer_sep.eta22 - mf12.eta22, '--', 'LineWidth', 1.45, 'Color', [0.85 0.33 0.10]);
hold on;
plot(x_over_lambda, creamer_sep.eta22 - mf12.eta22, ':', 'LineWidth', 2.0, 'Color', [0.47 0.67 0.19]);
yline(0, ':', 'Color', [0.4 0.4 0.4]);
grid on;
xlim(xlim_wave);
xlabel('x / \lambda_p');
ylabel('\Delta \eta_{22}');
title('\eta_{22} residual');
legend({'Creamer(H3)-MF12', 'Creamer(H3+H4)-MF12'}, 'Location', 'best');

nexttile;
eta22_mf12_pos = spec.mf12_eta22(pos_mask);
eta22_creamer_pos = spec.creamer_eta22(pos_mask);
semilogy(kpos(spec_mask), eta22_mf12_pos(spec_mask) + eps, '-', 'LineWidth', 2.2, 'Color', [0.00 0.45 0.74]);
hold on;
semilogy(kpos(spec_mask), eta22_creamer_pos(spec_mask) + eps, '--', 'LineWidth', 1.45, 'Color', [0.85 0.33 0.10]);
semilogy(kpos(spec_mask), eta22_creamer_pos(spec_mask) + eps, ':', 'LineWidth', 2.0, 'Color', [0.47 0.67 0.19]);
for kv = 1:6
    xline(kv, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
end
grid on;
xlim([0, 6]);
ylim(local_semilogy_limits([
    eta22_mf12_pos(spec_mask), ...
    eta22_creamer_pos(spec_mask), ...
    eta22_creamer_pos(spec_mask)]));
xlabel('k / k_p');
ylabel('Amplitude');
title('\eta_{22} spectrum');
legend({'MF12', 'Creamer(H3)', 'Creamer(H3+H4)'}, 'Location', 'best');

nexttile;
plot(x_over_lambda, mf12.eta33, '-', 'LineWidth', 2.3, 'Color', [0.00 0.45 0.74]);
hold on;
plot(x_over_lambda, creamer.eta33_h3, '--', 'LineWidth', 1.45, 'Color', [0.85 0.33 0.10]);
plot(x_over_lambda, creamer.eta33_h3h4, ':', 'LineWidth', 2.0, 'Color', [0.47 0.67 0.19]);
plot(x_over_lambda, env.mf12_eta33, '-', 'LineWidth', 1.9, 'Color', [0.45 0.70 0.90]);
plot(x_over_lambda, env.creamer_h3_eta33, '--', 'LineWidth', 1.7, 'Color', [0.95 0.62 0.40]);
plot(x_over_lambda, env.creamer_h3h4_eta33, ':', 'LineWidth', 2.0, 'Color', [0.70 0.84 0.42]);
grid on;
xlim(xlim_wave);
xlabel('x / \lambda_p');
ylabel('\eta_{33}');
title('\eta_{33}');
legend({'MF12', 'Creamer(H3)', 'Creamer(H3+H4)', ...
    'MF12 envelope', 'Creamer(H3) envelope', 'Creamer(H3+H4) envelope'}, ...
    'Location', 'best');

nexttile;
plot(x_over_lambda, creamer.eta33_h3 - mf12.eta33, '--', 'LineWidth', 1.25, 'Color', [0.85 0.33 0.10]);
hold on;
plot(x_over_lambda, creamer.eta33_h3h4 - mf12.eta33, '-', 'LineWidth', 1.25, 'Color', [0.47 0.67 0.19]);
yline(0, ':', 'Color', [0.4 0.4 0.4]);
grid on;
xlim(xlim_wave);
xlabel('x / \lambda_p');
ylabel('\Delta \eta_{33}');
title('\eta_{33} residual');
legend({'Creamer(H3)-MF12', 'Creamer(H3+H4)-MF12'}, 'Location', 'best');

nexttile;
mf12_pos = spec.mf12_eta33(pos_mask);
h3_pos = spec.creamer_h3_eta33(pos_mask);
h3h4_pos = spec.creamer_h3h4_eta33(pos_mask);
semilogy(kpos(spec_mask), mf12_pos(spec_mask) + eps, '-', 'LineWidth', 2.3, 'Color', [0.00 0.45 0.74]);
hold on;
semilogy(kpos(spec_mask), h3_pos(spec_mask) + eps, '--', 'LineWidth', 1.45, 'Color', [0.85 0.33 0.10]);
semilogy(kpos(spec_mask), h3h4_pos(spec_mask) + eps, ':', 'LineWidth', 2.0, 'Color', [0.47 0.67 0.19]);
for kv = 1:6
    xline(kv, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
end
grid on;
xlim([0, 6]);
ylim(local_semilogy_limits([
    mf12_pos(spec_mask), ...
    h3_pos(spec_mask), ...
    h3h4_pos(spec_mask)]));
xlabel('k / k_p');
ylabel('Amplitude');
title('\eta_{33} spectrum');
legend({'MF12', 'Creamer(H3)', 'Creamer(H3+H4)'}, 'Location', 'best');

sgtitle(sprintf(['Gaussian finite-depth eta22/eta33 | k_p=%.4g, Ak_p=%.3f, \\sigma_k/k_p=%.3f, ', ...
    'k_p h=%.3g, N_x=%d, ppw=%.2f, L/\\lambda_p=%.1f, keep=%d, backend=%s'], ...
    kp, Akp, bandwidth_ratio, kp*h, Nx, points_per_wavelength, ...
    domain_wavelengths, max_triad_active_modes, triad_backend));
annotation(fig, 'textbox', [0.56 0.905 0.40 0.032], ...
    'String', sprintf('MF12 active modes = %d | Creamer H3/H3+H4 parent modes = %d (same set)', ...
    mf12.active_mode_count, triad_active_count), ...
    'FitBoxToText', 'off', 'EdgeColor', 'none', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'middle', 'FontSize', 9, 'Color', [0.2 0.2 0.2]);

out_png = fullfile(out_dir, sprintf('unidirectional_gaussian_eta22_eta33_kp%04d_bw%04d_Akp%04d_kph%03d_Nx%d_keep%d_%s.png', ...
    round(1000 * kp), round(10000 * bandwidth_ratio), round(10000 * Akp), round(100 * kp * h), ...
    Nx, max_triad_active_modes, triad_backend));
exportgraphics(fig, out_png, 'Resolution', 170);

fprintf('\nSaved Gaussian eta33 comparison figure to:\n  %s\n', out_png);
fprintf('MF12 runtime                            = %.3f s\n', mf12_runtime_s);
fprintf('Creamer eta22 separation runtime         = %.3f s\n', sep_runtime_s);
fprintf('Creamer triad runtime                   = %.3f s\n', creamer_runtime_s);
fprintf('Bandwidth sigma_k/kp                    = %.12g\n', bandwidth_ratio);
fprintf('Points per wavelength                   = %.6g\n', points_per_wavelength);
fprintf('Domain length / lambda_p                = %.6g\n', domain_wavelengths);
fprintf('Triad backend                           = %s\n', triad_backend);
fprintf('Triad active parent modes               = %d\n', max_triad_active_modes);
fprintf('Peak-bin eta33 amplitude near 3kp:\n');
fprintf('  MF12 eta33                            = %.12g\n', mf12_3kp);
fprintf('  Creamer(H3) eta33                     = %.12g\n', creamer_h3_3kp);
fprintf('  Creamer(H3+H4) eta33                  = %.12g\n', creamer_h3h4_3kp);
fprintf('  Creamer(H3+H4)/MF12                   = %.12g\n', creamer_h3h4_3kp / mf12_3kp);

function lim = local_semilogy_limits(vals)
vals = vals(isfinite(vals) & vals > 0);
if isempty(vals)
    lim = [1e-12, 1];
    return;
end
peak_val = max(vals);
floor_val = max(min(vals), peak_val * 1e-4);
if floor_val >= peak_val
    floor_val = peak_val / 10;
end
lim = [floor_val, 1.2 * peak_val];
end

function n = local_triad_active_count(creamer)
n = NaN;
if isfield(creamer, 'diagnostics')
    if isfield(creamer.diagnostics, 'positive_modes')
        n = numel(creamer.diagnostics.positive_modes);
    elseif isfield(creamer.diagnostics, 'n_positive')
        n = creamer.diagnostics.n_positive;
    end
end
end
