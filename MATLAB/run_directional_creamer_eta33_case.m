% Compare four-phase-separated Creamer eta33 against MF12 eta33.

clearvars -except akp_index alpha_index chi_index spread_index max_active_modes n_lambda_steps lambda_flow_model n_picard_iters;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
data_file = fullfile(repo_root, ...
    'directional_waves_multi_Akp0.020-0.180_Alpha1-8_Chi0-0_Spread5-30_adaptive100pc_20260409_154705.mat');

S = load(data_file, 'results', 'x_vec', 'y_vec', 'parameters');

if ~exist('akp_index', 'var') || isempty(akp_index)
    akp_index = 2;
end
if ~exist('alpha_index', 'var') || isempty(alpha_index)
    alpha_index = 2;
end
if ~exist('chi_index', 'var') || isempty(chi_index)
    chi_index = 1;
end
if ~exist('spread_index', 'var') || isempty(spread_index)
    spread_index = 6;
end
if ~exist('max_active_modes', 'var') || isempty(max_active_modes)
    max_active_modes = 3000;
end
if ~exist('n_lambda_steps', 'var') || isempty(n_lambda_steps)
    n_lambda_steps = 12;
end
if ~exist('lambda_flow_model', 'var') || isempty(lambda_flow_model)
    lambda_flow_model = 'canonical_pair';
end
if ~exist('n_picard_iters', 'var') || isempty(n_picard_iters)
    n_picard_iters = 4;
end

eta_lin = S.results.eta_results{akp_index, alpha_index, chi_index, spread_index};
x_vec = S.x_vec;
y_vec = S.y_vec;
akp_val = S.parameters.Akp_values(akp_index);
alpha_val = S.parameters.Alpha_values(alpha_index);
chi_deg = S.parameters.propagation_directions(chi_index);
spread_deg = S.parameters.spread_angles(spread_index);

cfg = struct();
cfg.g = 9.81;
cfg.energy_fraction = 0.9999999999;
cfg.max_active_modes = max_active_modes;
cfg.n_lambda_steps = n_lambda_steps;
cfg.n_picard_iters = n_picard_iters;
cfg.lambda_flow_model = lambda_flow_model;
cfg.propagation_direction_deg = chi_deg;
cfg.preserve_mean = true;
cfg.verbose = true;

validation_file = fullfile(repo_root, ...
    'directional_validation_results_mf12_linear_groups_linear_groups_kd50_mc1500_fixakp_20260409_190137.mat');
validation_case = load_validation_third_order_case( ...
    validation_file, alpha_val, chi_deg, spread_deg, 50, akp_val);

t_sep = tic;
creamer_sep = creamer_four_phase_separation(eta_lin, x_vec, y_vec, cfg);
creamer_runtime_s = toc(t_sep);

[X_target, Y_target] = meshgrid(x_vec, y_vec);
[X_val, Y_val] = meshgrid(validation_case.x(:).', validation_case.y(:));
eta33_val = interp2(X_val, Y_val, validation_case.eta33, X_target, Y_target, 'linear', 0);
eta31_val = interp2(X_val, Y_val, validation_case.eta31, X_target, Y_target, 'linear', 0);

[~, iy0] = min(abs(y_vec - S.parameters.y_focus));
eta_lin_analytic_x = hilbert(eta_lin.').';
row_envelope_peak = max(abs(eta_lin_analytic_x), [], 2);
target_off_centerline = 0.5 * max(row_envelope_peak);
candidate_mask = (y_vec(:) >= S.parameters.y_focus) & ((1:numel(y_vec)).' ~= iy0);
candidate_idx = find(candidate_mask);
[~, best_off_idx_local] = min(abs(row_envelope_peak(candidate_idx) - target_off_centerline));
iy_off = candidate_idx(best_off_idx_local);
y_off = y_vec(iy_off);

center_eta33_val = eta33_val(iy0, :);
center_eta33_creamer = creamer_sep.eta33(iy0, :);
off_eta33_val = eta33_val(iy_off, :);
off_eta33_creamer = creamer_sep.eta33(iy_off, :);

dx = mean(diff(x_vec));
nx = numel(x_vec);
Lx = dx * nx;
kp = S.parameters.kp;
lambda_p = 2*pi/kp;
x_over_lambda = x_vec / lambda_p;
x_focus_n = S.parameters.x_focus / lambda_p;
xlim_wave = [x_focus_n - 4.0, x_focus_n + 4.0];

mx = [0:(nx/2), (-nx/2 + 1):-1];
kx = (2*pi/Lx) * mx;
pos_mask = (mx >= 0);
kx_over_kp = kx(pos_mask) / kp;

spec_center_eta33_val = abs(fft(center_eta33_val) / nx);
spec_center_eta33_creamer = abs(fft(center_eta33_creamer) / nx);
spec_off_eta33_val = abs(fft(off_eta33_val) / nx);
spec_off_eta33_creamer = abs(fft(off_eta33_creamer) / nx);

out_dir = fullfile(repo_root, 'MATLAB', 'output');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

tag = sprintf('alpha%d_chi%d_spread%d_Akp%04d_mc%d_nl%d_%s_pi%d', ...
    round(alpha_val), round(chi_deg), round(spread_deg), round(10000 * akp_val), ...
    round(max_active_modes), round(n_lambda_steps), lambda_flow_model, round(n_picard_iters));

fig = figure('Color', 'w', 'Position', [120 120 1450 980]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(x_over_lambda, center_eta33_val, '-.', 'LineWidth', 1.3);
hold on;
plot(x_over_lambda, center_eta33_creamer, '--', 'LineWidth', 1.3);
grid on;
xlim(xlim_wave);
legend({'MF12 \eta_{33}', 'Creamer \eta_{33}'}, 'Location', 'best');
title(sprintf('\\eta_{33} centerline at y = %.6g', y_vec(iy0)));
xlabel('x / \lambda_p');
ylabel('\eta');

nexttile;
plot(x_over_lambda, center_eta33_val, '-.', 'LineWidth', 1.3);
hold on;
plot(x_over_lambda, center_eta33_creamer, '--', 'LineWidth', 1.3);
grid on;
xlim(xlim_wave);
legend({'MF12 \eta_{33}', 'Creamer \eta_{33}'}, 'Location', 'best');
title(sprintf('\\eta_{33} centerline zoom at y = %.6g', y_vec(iy0)));
xlabel('x / \lambda_p');
ylabel('\eta');

nexttile;
semilogy(kx_over_kp, spec_center_eta33_val(pos_mask) + eps, '-.', 'LineWidth', 1.2);
hold on;
semilogy(kx_over_kp, spec_center_eta33_creamer(pos_mask) + eps, '-.', 'LineWidth', 1.2, 'Color', [0.6 0.2 0.6]);
grid on;
xlim([0, 5]);
legend({'MF12 \eta_{33}', 'Creamer \eta_{33}'}, 'Location', 'best');
title('Centerline spatial spectrum');
xlabel('k_x / k_p');
ylabel('Amplitude');

nexttile;
plot(x_over_lambda, off_eta33_val, '-.', 'LineWidth', 1.3);
hold on;
plot(x_over_lambda, off_eta33_creamer, '--', 'LineWidth', 1.3);
grid on;
xlim(xlim_wave);
legend({'MF12 \eta_{33}', 'Creamer \eta_{33}'}, 'Location', 'best');
title(sprintf('\\eta_{33} off-center at y = %.6g', y_off));
xlabel('x / \lambda_p');
ylabel('\eta');

nexttile;
plot(x_over_lambda, off_eta33_val, '-.', 'LineWidth', 1.3);
hold on;
plot(x_over_lambda, off_eta33_creamer, '--', 'LineWidth', 1.3);
grid on;
xlim(xlim_wave);
legend({'MF12 \eta_{33}', 'Creamer \eta_{33}'}, 'Location', 'best');
title(sprintf('\\eta_{33} off-center zoom at y = %.6g', y_off));
xlabel('x / \lambda_p');
ylabel('\eta');

nexttile;
semilogy(kx_over_kp, spec_off_eta33_val(pos_mask) + eps, '-.', 'LineWidth', 1.2);
hold on;
semilogy(kx_over_kp, spec_off_eta33_creamer(pos_mask) + eps, '-.', 'LineWidth', 1.2, 'Color', [0.6 0.2 0.6]);
grid on;
xlim([0, 5]);
legend({'MF12 \eta_{33}', 'Creamer \eta_{33}'}, 'Location', 'best');
title('Off-center spatial spectrum');
xlabel('k_x / k_p');
ylabel('Amplitude');

sgtitle(sprintf(['Third-order comparison: \\eta_{33} and spectra | ', ...
    'A k_p = %.3f, \\chi = %g^\\circ, \\alpha = %g, spread = %g^\\circ, ', ...
    'modes = %d, N_\\lambda = %d, model = %s'], ...
    akp_val, chi_deg, alpha_val, spread_deg, max_active_modes, n_lambda_steps, lambda_flow_model));

out_png = fullfile(out_dir, ['directional_creamer_eta33_' tag '.png']);
exportgraphics(fig, out_png, 'Resolution', 160);

fprintf('\nSaved third-order comparison figure to:\n  %s\n', out_png);
fprintf('Creamer four-phase runtime             = %.3f s\n', creamer_runtime_s);
fprintf('Validation source                      = %s\n', validation_file);
fprintf('Configured max active modes            = %d\n', max_active_modes);
fprintf('Configured lambda steps                = %d\n', n_lambda_steps);
fprintf('Configured lambda model                = %s\n', lambda_flow_model);
fprintf('Configured Picard iterations           = %d\n', n_picard_iters);
fprintf('Centerline max |MF12 eta33|            = %.6g\n', max(abs(center_eta33_val)));
fprintf('Centerline max |Creamer eta33|         = %.6g\n', max(abs(center_eta33_creamer)));
fprintf('Off-center max |MF12 eta33|            = %.6g\n', max(abs(off_eta33_val)));
fprintf('Off-center max |Creamer eta33|         = %.6g\n', max(abs(off_eta33_creamer)));
