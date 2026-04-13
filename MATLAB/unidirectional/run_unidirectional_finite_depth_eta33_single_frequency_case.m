% Compare finite-depth single-frequency eta33:
%   MF12 eta33, Creamer(H3) eta33, Creamer(H3+H4) eta33.

clearvars -except Akp Nx k0 h k0h n_periods min_active_modes max_active_modes mf12_max_active_modes n_lambda_steps lambda_flow_model t phishift;
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
if ~exist('Nx', 'var') || isempty(Nx)
    Nx = 2048;
end
if ~exist('k0', 'var') || isempty(k0)
    k0 = 0.0279;
end
if ~exist('k0h', 'var') || isempty(k0h)
    k0h = 2;
end
if ~exist('h', 'var') || isempty(h)
    h = k0h / k0;
else
    k0h = k0 * h;
end
if ~exist('n_periods', 'var') || isempty(n_periods)
    n_periods = 32;
end
if ~exist('max_active_modes', 'var') || isempty(max_active_modes)
    max_active_modes = 256;
end
if ~exist('mf12_max_active_modes', 'var') || isempty(mf12_max_active_modes)
    mf12_max_active_modes = 1;
end
if ~exist('min_active_modes', 'var') || isempty(min_active_modes)
    min_active_modes = 16;
end
if ~exist('n_lambda_steps', 'var') || isempty(n_lambda_steps)
    n_lambda_steps = 24;
end
if ~exist('lambda_flow_model', 'var') || isempty(lambda_flow_model)
    lambda_flow_model = 'canonical_pair';
end
if ~exist('t', 'var') || isempty(t)
    t = 0;
end
if ~exist('phishift', 'var') || isempty(phishift)
    phishift = 0;
end

g = 9.81;
lambda0 = 2*pi/k0;
Lx = n_periods * lambda0;
x_vec = (0:Nx-1) * (Lx / Nx);
A = Akp / k0;
eta_lin = A * cos(k0 * x_vec + phishift);

cfg = struct();
cfg.g = g;
cfg.depth_h = h;
cfg.energy_fraction = 1;
cfg.min_active_modes = min_active_modes;
cfg.max_active_modes = max_active_modes;
cfg.n_lambda_steps = n_lambda_steps;
cfg.lambda_flow_model = lambda_flow_model;
cfg.propagation_direction_deg = 0;
cfg.preserve_mean = true;
cfg.verbose = true;

mf12_cfg = struct();
mf12_cfg.g = g;
mf12_cfg.h = h;
mf12_cfg.t = t;
mf12_cfg.energy_fraction = 1;
mf12_cfg.max_active_modes = mf12_max_active_modes;
mf12_cfg.max_order = 3;

t_mf12 = tic;
mf12 = mf12_from_linear_focus_1d(eta_lin, x_vec, mf12_cfg);
mf12_runtime_s = toc(t_mf12);

t_creamer = tic;
creamer_sep = creamer_four_phase_separation_1d(eta_lin, x_vec, cfg);
creamer_runtime_s = toc(t_creamer);

creamer_h3h4 = creamer_eta33_h3h4_single_frequency_1d(eta_lin, x_vec, h);

dx = mean(diff(x_vec));
mx = [0:(Nx/2), (-Nx/2 + 1):-1];
k = (2*pi/Lx) * mx;
k_over_k0 = k / k0;
pos_mask = mx >= 0;

spec_mf12_eta33 = abs(fft(mf12.eta33) / Nx);
spec_creamer_h3_eta33 = abs(fft(creamer_sep.eta33) / Nx);
spec_creamer_h3_self_eta33 = abs(fft(creamer_h3h4.eta33_h3) / Nx);
spec_creamer_h3h4_eta33 = abs(fft(creamer_h3h4.eta33_h3h4) / Nx);

[~, idx_3k0_full] = min(abs(k_over_k0 - 3));

eta33_mf12_3k0 = spec_mf12_eta33(idx_3k0_full);
eta33_creamer_h3_3k0 = spec_creamer_h3_eta33(idx_3k0_full);
eta33_creamer_h3_self_3k0 = spec_creamer_h3_self_eta33(idx_3k0_full);
eta33_creamer_h3h4_3k0 = spec_creamer_h3h4_eta33(idx_3k0_full);

sigma = tanh(k0h);
c33_h3 = local_creamer_h3_c33(sigma);
c33_h3h4 = local_stokes_c33(sigma);
eta33_expected_h3_amp = c33_h3 * Akp^3 / k0;
eta33_expected_h3h4_amp = c33_h3h4 * Akp^3 / k0;

result = struct();
result.Akp = Akp;
result.k0 = k0;
result.k0h = k0h;
result.h = h;
result.x_vec = x_vec;
result.lambda_p = lambda0;
result.xlim_wave = [0, 4];
result.k_over_kp = k_over_k0(pos_mask);
result.points_per_wavelength = lambda0 / dx;
result.max_active_modes = max_active_modes;
result.min_active_modes = min_active_modes;
result.n_lambda_steps = n_lambda_steps;
result.lambda_flow_model = lambda_flow_model;
result.creamer_backend = 'matlab';
result.mf12 = struct('eta33', mf12.eta33);
result.creamer = struct( ...
    'eta33_flow_h3', creamer_sep.eta33, ...
    'eta33', creamer_h3h4.eta33_h3, ...
    'eta33_h3_self', creamer_h3h4.eta33_h3, ...
    'eta33_h3h4', creamer_h3h4.eta33_h3h4);
result.spec = struct( ...
    'mf12_eta33', spec_mf12_eta33(pos_mask), ...
    'creamer_eta33_flow_h3', spec_creamer_h3_eta33(pos_mask), ...
    'creamer_eta33', spec_creamer_h3_self_eta33(pos_mask), ...
    'creamer_eta33_h3_self', spec_creamer_h3_self_eta33(pos_mask), ...
    'creamer_eta33_h3h4', spec_creamer_h3h4_eta33(pos_mask));

out_dir = fullfile(repo_root, 'MATLAB', 'output', 'unidirectional');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

fig = figure('Color', 'w', 'Position', [120 120 1300 820]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

x_over_lambda = x_vec / lambda0;
spec_mask = result.k_over_kp <= 6;

nexttile;
plot(x_over_lambda, mf12.eta33, '-.', 'LineWidth', 1.35);
hold on;
plot(x_over_lambda, creamer_h3h4.eta33_h3, '--', 'LineWidth', 1.25);
plot(x_over_lambda, creamer_h3h4.eta33_h3h4, '-', 'LineWidth', 1.25);
grid on;
xlim(result.xlim_wave);
xlabel('x / \lambda_0');
ylabel('\eta_{33}');
title('Single-frequency \eta_{33}');
legend({'MF12', 'Creamer(H3)', 'Creamer(H3+H4)'}, 'Location', 'best');

nexttile;
plot(x_over_lambda, creamer_h3h4.eta33_h3 - mf12.eta33, '--', 'LineWidth', 1.25);
hold on;
plot(x_over_lambda, creamer_h3h4.eta33_h3h4 - mf12.eta33, '-', 'LineWidth', 1.25);
yline(0, ':');
grid on;
xlim(result.xlim_wave);
xlabel('x / \lambda_0');
ylabel('\Delta \eta_{33}');
title('Residuals');
legend({'Creamer(H3)-MF12', 'Creamer(H3+H4)-MF12'}, 'Location', 'best');

nexttile;
semilogy(result.k_over_kp(spec_mask), result.spec.mf12_eta33(spec_mask) + eps, '-.', 'LineWidth', 1.25);
hold on;
semilogy(result.k_over_kp(spec_mask), result.spec.creamer_eta33_h3_self(spec_mask) + eps, '--', 'LineWidth', 1.25);
semilogy(result.k_over_kp(spec_mask), result.spec.creamer_eta33_h3h4(spec_mask) + eps, '-', 'LineWidth', 1.25);
xline(3, ':');
grid on;
xlim([0, 6]);
xlabel('k / k_0');
ylabel('Amplitude');
title('Spectrum');
legend({'MF12', 'Creamer(H3)', 'Creamer(H3+H4)'}, 'Location', 'best');

nexttile;
bar(categorical({'MF12','Creamer(H3)','Creamer(H3+H4)'}), ...
    [eta33_mf12_3k0, eta33_creamer_h3_self_3k0, eta33_creamer_h3h4_3k0]);
grid on;
ylabel('3k_0 amplitude');
title('Peak-bin \eta_{33}');

sgtitle(sprintf('Finite-depth single-frequency eta33 | Ak=%.4g, k h=%.4g, Nx=%d, Nlambda=%d', ...
    Akp, k0h, Nx, n_lambda_steps));

out_png = fullfile(out_dir, sprintf('unidirectional_finite_depth_eta33_single_frequency_Akp%04d_kh%03d_Nx%d_nl%d.png', ...
    round(10000*Akp), round(100*k0h), Nx, n_lambda_steps));
exportgraphics(fig, out_png, 'Resolution', 170);

fprintf('\nSaved single-frequency eta33 comparison figure to:\n  %s\n', out_png);
fprintf('MF12 runtime                            = %.3f s\n', mf12_runtime_s);
fprintf('Creamer four-phase runtime              = %.3f s\n', creamer_runtime_s);
fprintf('A k0                                    = %.12g\n', Akp);
fprintf('k0 h                                    = %.12g\n', k0h);
fprintf('sigma=tanh(k0 h)                        = %.12g\n', sigma);
fprintf('C33 Creamer(H3)                         = %.12g\n', c33_h3);
fprintf('C33 Creamer(H3+H4)/Stokes               = %.12g\n', c33_h3h4);
fprintf('Expected physical eta33 amplitudes:\n');
fprintf('  Creamer(H3)                           = %.12g\n', eta33_expected_h3_amp);
fprintf('  Creamer(H3+H4)/Stokes                 = %.12g\n', eta33_expected_h3h4_amp);
fprintf('Peak-bin eta33 amplitude near 3k0:\n');
fprintf('  MF12 eta33                            = %.12g\n', eta33_mf12_3k0);
fprintf('  Current flow-separated Creamer(H3)     = %.12g\n', eta33_creamer_h3_3k0);
fprintf('  Creamer(H3) eta33                     = %.12g\n', eta33_creamer_h3_self_3k0);
fprintf('  Creamer(H3+H4) eta33                  = %.12g\n', eta33_creamer_h3h4_3k0);

function c = local_creamer_h3_c33(sigma)
c = (3 * (1 + sigma.^2) .* ...
    (54 + 39*sigma.^2 + 70*sigma.^4 - 28*sigma.^6 - 4*sigma.^8 - 3*sigma.^10)) ./ ...
    (64 * sigma.^6 .* (9 + 14*sigma.^2 + 9*sigma.^4));
end

function c = local_stokes_c33(sigma)
c = (27 - 9*sigma.^2 + 9*sigma.^4 - 3*sigma.^6) ./ (64 * sigma.^6);
end
