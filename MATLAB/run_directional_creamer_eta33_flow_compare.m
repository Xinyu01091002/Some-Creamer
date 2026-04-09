% Compare eta33 from two lambda-flow closures against MF12 on one case.

clearvars -except akp_index alpha_index chi_index spread_index max_active_modes n_lambda_steps;
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
    spread_index = 1;
end
if ~exist('max_active_modes', 'var') || isempty(max_active_modes)
    max_active_modes = 2500;
end
if ~exist('n_lambda_steps', 'var') || isempty(n_lambda_steps)
    n_lambda_steps = 12;
end

eta_lin = S.results.eta_results{akp_index, alpha_index, chi_index, spread_index};
x_vec = S.x_vec;
y_vec = S.y_vec;
akp_val = S.parameters.Akp_values(akp_index);
alpha_val = S.parameters.Alpha_values(alpha_index);
chi_deg = S.parameters.propagation_directions(chi_index);
spread_deg = S.parameters.spread_angles(spread_index);

validation_file = fullfile(repo_root, ...
    'directional_validation_results_mf12_linear_groups_linear_groups_kd50_mc1500_fixakp_20260409_190137.mat');
validation_case = load_validation_third_order_case( ...
    validation_file, alpha_val, chi_deg, spread_deg, 50, akp_val);

cfg_base = struct();
cfg_base.g = 9.81;
cfg_base.energy_fraction = 0.9999999999;
cfg_base.max_active_modes = max_active_modes;
cfg_base.n_lambda_steps = n_lambda_steps;
cfg_base.propagation_direction_deg = chi_deg;
cfg_base.preserve_mean = true;
cfg_base.verbose = true;

cfg_legacy = cfg_base;
cfg_legacy.lambda_flow_model = 'legacy_zeta_only';

cfg_pair = cfg_base;
cfg_pair.lambda_flow_model = 'canonical_pair';

t_legacy = tic;
sep_legacy = creamer_four_phase_separation(eta_lin, x_vec, y_vec, cfg_legacy);
runtime_legacy = toc(t_legacy);

t_pair = tic;
sep_pair = creamer_four_phase_separation(eta_lin, x_vec, y_vec, cfg_pair);
runtime_pair = toc(t_pair);

[X_target, Y_target] = meshgrid(x_vec, y_vec);
[X_val, Y_val] = meshgrid(validation_case.x(:).', validation_case.y(:));
eta33_val = interp2(X_val, Y_val, validation_case.eta33, X_target, Y_target, 'linear', 0);

[~, iy0] = min(abs(y_vec - S.parameters.y_focus));
eta_lin_analytic_x = hilbert(eta_lin.').';
row_envelope_peak = max(abs(eta_lin_analytic_x), [], 2);
target_off_centerline = 0.5 * max(row_envelope_peak);
candidate_mask = (y_vec(:) >= S.parameters.y_focus) & ((1:numel(y_vec)).' ~= iy0);
candidate_idx = find(candidate_mask);
[~, best_off_idx_local] = min(abs(row_envelope_peak(candidate_idx) - target_off_centerline));
iy_off = candidate_idx(best_off_idx_local);
y_off = y_vec(iy_off);

center_val = eta33_val(iy0, :);
center_legacy = sep_legacy.eta33(iy0, :);
center_pair = sep_pair.eta33(iy0, :);
off_val = eta33_val(iy_off, :);
off_legacy = sep_legacy.eta33(iy_off, :);
off_pair = sep_pair.eta33(iy_off, :);

dx = mean(diff(x_vec));
nx = numel(x_vec);
Lx = dx * nx;
kp = S.parameters.kp;
lambda_p = 2 * pi / kp;
x_over_lambda = x_vec / lambda_p;
x_focus_n = S.parameters.x_focus / lambda_p;
xlim_wave = [x_focus_n - 4.0, x_focus_n + 4.0];

mx = [0:(nx/2), (-nx/2 + 1):-1];
kx = (2*pi/Lx) * mx;
pos_mask = (mx >= 0);
kx_over_kp = kx(pos_mask) / kp;

spec_center_val = abs(fft(center_val) / nx);
spec_center_legacy = abs(fft(center_legacy) / nx);
spec_center_pair = abs(fft(center_pair) / nx);
spec_off_val = abs(fft(off_val) / nx);
spec_off_legacy = abs(fft(off_legacy) / nx);
spec_off_pair = abs(fft(off_pair) / nx);

out_dir = fullfile(repo_root, 'MATLAB', 'output');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

tag = sprintf('alpha%d_chi%d_spread%d_Akp%04d', ...
    round(alpha_val), round(chi_deg), round(spread_deg), round(10000 * akp_val));

fig = figure('Color', 'w', 'Position', [120 120 1500 980]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(x_over_lambda, center_val, 'k-.', 'LineWidth', 1.3);
hold on;
plot(x_over_lambda, center_legacy, '--', 'LineWidth', 1.2, 'Color', [0.85 0.33 0.10]);
plot(x_over_lambda, center_pair, '-', 'LineWidth', 1.2, 'Color', [0 0.45 0.74]);
grid on;
xlim(xlim_wave);
legend({'MF12 \eta_{33}', 'Legacy zeta-only', 'Canonical-pair'}, 'Location', 'best');
title(sprintf('Centerline at y = %.6g', y_vec(iy0)));
xlabel('x / \lambda_p');
ylabel('\eta_{33}');

nexttile;
plot(x_over_lambda, off_val, 'k-.', 'LineWidth', 1.3);
hold on;
plot(x_over_lambda, off_legacy, '--', 'LineWidth', 1.2, 'Color', [0.85 0.33 0.10]);
plot(x_over_lambda, off_pair, '-', 'LineWidth', 1.2, 'Color', [0 0.45 0.74]);
grid on;
xlim(xlim_wave);
legend({'MF12 \eta_{33}', 'Legacy zeta-only', 'Canonical-pair'}, 'Location', 'best');
title(sprintf('Off-center at y = %.6g', y_off));
xlabel('x / \lambda_p');
ylabel('\eta_{33}');

nexttile;
semilogy(kx_over_kp, spec_center_val(pos_mask) + eps, 'k-.', 'LineWidth', 1.2);
hold on;
semilogy(kx_over_kp, spec_center_legacy(pos_mask) + eps, '--', 'LineWidth', 1.2, 'Color', [0.85 0.33 0.10]);
semilogy(kx_over_kp, spec_center_pair(pos_mask) + eps, '-', 'LineWidth', 1.2, 'Color', [0 0.45 0.74]);
grid on;
xlim([0, 5]);
legend({'MF12 \eta_{33}', 'Legacy zeta-only', 'Canonical-pair'}, 'Location', 'best');
title('Centerline spectrum');
xlabel('k_x / k_p');
ylabel('Amplitude');

nexttile;
plot(x_over_lambda, center_pair - center_legacy, '-', 'LineWidth', 1.2, 'Color', [0.49 0.18 0.56]);
grid on;
xlim(xlim_wave);
title('Canonical-pair minus legacy, centerline');
xlabel('x / \lambda_p');
ylabel('\Delta \eta_{33}');

nexttile;
plot(x_over_lambda, off_pair - off_legacy, '-', 'LineWidth', 1.2, 'Color', [0.49 0.18 0.56]);
grid on;
xlim(xlim_wave);
title('Canonical-pair minus legacy, off-center');
xlabel('x / \lambda_p');
ylabel('\Delta \eta_{33}');

nexttile;
semilogy(kx_over_kp, spec_off_val(pos_mask) + eps, 'k-.', 'LineWidth', 1.2);
hold on;
semilogy(kx_over_kp, spec_off_legacy(pos_mask) + eps, '--', 'LineWidth', 1.2, 'Color', [0.85 0.33 0.10]);
semilogy(kx_over_kp, spec_off_pair(pos_mask) + eps, '-', 'LineWidth', 1.2, 'Color', [0 0.45 0.74]);
grid on;
xlim([0, 5]);
legend({'MF12 \eta_{33}', 'Legacy zeta-only', 'Canonical-pair'}, 'Location', 'best');
title('Off-center spectrum');
xlabel('k_x / k_p');
ylabel('Amplitude');

sgtitle(sprintf(['Flow-model comparison for \\eta_{33} | ', ...
    'A k_p = %.3f, \\chi = %g^\\circ, \\alpha = %g, spread = %g^\\circ, ', ...
    'modes = %d, N_\\lambda = %d'], ...
    akp_val, chi_deg, alpha_val, spread_deg, max_active_modes, n_lambda_steps));

out_png = fullfile(out_dir, ['directional_creamer_eta33_flow_compare_' tag '.png']);
exportgraphics(fig, out_png, 'Resolution', 160);

fprintf('\nSaved flow-model comparison figure to:\n  %s\n', out_png);
fprintf('Validation source                      = %s\n', validation_file);
fprintf('Configured max active modes            = %d\n', max_active_modes);
fprintf('Configured lambda steps                = %d\n', n_lambda_steps);
fprintf('Legacy runtime                         = %.3f s\n', runtime_legacy);
fprintf('Canonical-pair runtime                 = %.3f s\n', runtime_pair);
fprintf('Centerline max |MF12 eta33|            = %.6g\n', max(abs(center_val)));
fprintf('Centerline max |legacy eta33|          = %.6g\n', max(abs(center_legacy)));
fprintf('Centerline max |canonical eta33|       = %.6g\n', max(abs(center_pair)));
fprintf('Off-center max |MF12 eta33|            = %.6g\n', max(abs(off_val)));
fprintf('Off-center max |legacy eta33|          = %.6g\n', max(abs(off_legacy)));
fprintf('Off-center max |canonical eta33|       = %.6g\n', max(abs(off_pair)));
