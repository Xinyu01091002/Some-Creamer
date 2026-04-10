% Compare MF12 and Creamer second-order phi20/phi22 fields.

clearvars -except akp_index alpha_index chi_index spread_index;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
data_file = fullfile(repo_root, ...
    'directional_waves_multi_Akp0.020-0.180_Alpha1-8_Chi0-0_Spread5-30_adaptive100pc_20260409_154705.mat');

if ~isfile(data_file)
    error('run_directional_creamer_phi_case:MissingData', ...
        'Could not find data file: %s', data_file);
end

S = load(data_file, 'results', 'x_vec', 'y_vec', 'parameters');

if ~exist('akp_index', 'var') || isempty(akp_index)
    akp_index = 2;     % Akp = 0.12
end
if ~exist('alpha_index', 'var') || isempty(alpha_index)
    alpha_index = 2;   % alpha = 8
end
if ~exist('chi_index', 'var') || isempty(chi_index)
    chi_index = 1;     % chi = 0 deg
end
if ~exist('spread_index', 'var') || isempty(spread_index)
    spread_index = 1;  % spread = 5 deg
end

eta_lin = S.results.eta_results{akp_index, alpha_index, chi_index, spread_index};
x_vec = S.x_vec;
y_vec = S.y_vec;
akp_val = S.parameters.Akp_values(akp_index);
spread_deg = S.parameters.spread_angles(spread_index);
alpha_val = S.parameters.Alpha_values(alpha_index);
chi_deg = S.parameters.propagation_directions(chi_index);

cfg = struct();
cfg.g = 9.81;
cfg.energy_fraction = 0.99999999;
cfg.max_active_modes = 2000;
cfg.n_lambda_steps = 8;
cfg.propagation_direction_deg = chi_deg;
cfg.preserve_mean = true;
cfg.verbose = true;

validation_file = fullfile(repo_root, ...
    'directional_validation_results_mf12_linear_groups_linear_groups_kd50_mc1500_fixakp_20260409_190137.mat');
if ~isfile(validation_file)
    error('run_directional_creamer_phi_case:MissingValidationData', ...
        'Could not find validation data file: %s', validation_file);
end

[validation_case, validation_match] = match_validation_second_order_case( ...
    validation_file, S.parameters.Alpha_values, chi_deg, spread_deg, 50, ...
    eta_lin, x_vec, y_vec, akp_val);

[~, iy0] = min(abs(y_vec - S.parameters.y_focus));
eta_lin_analytic_x = hilbert(eta_lin.').';
row_envelope_peak = max(abs(eta_lin_analytic_x), [], 2);
target_off_centerline = 0.5 * max(row_envelope_peak);
candidate_mask = (y_vec(:) >= S.parameters.y_focus) & ((1:numel(y_vec)).' ~= iy0);
candidate_idx = find(candidate_mask);
[~, best_off_idx_local] = min(abs(row_envelope_peak(candidate_idx) - target_off_centerline));
iy_off = candidate_idx(best_off_idx_local);
y_off = y_vec(iy_off);

[X_target, Y_target] = meshgrid(x_vec, y_vec);
[X_val, Y_val] = meshgrid(validation_case.x(:).', validation_case.y(:));
phi20_val_interp = interp2(X_val, Y_val, validation_case.phi20, X_target, Y_target, 'linear', 0);
phi22_val_interp = interp2(X_val, Y_val, validation_case.phi22, X_target, Y_target, 'linear', 0);

t_creamer_sep = tic;
creamer_sep = creamer_four_phase_separation(eta_lin, x_vec, y_vec, cfg);
creamer_separation_runtime_s = toc(t_creamer_sep);
phi20_creamer = creamer_sep.phi20;
phi22_creamer = creamer_sep.phi22;

centerline_phi20_val = phi20_val_interp(iy0, :);
centerline_phi22_val = phi22_val_interp(iy0, :);
centerline_creamer_phi20 = phi20_creamer(iy0, :);
centerline_creamer_phi22 = phi22_creamer(iy0, :);

offline_phi20_val = phi20_val_interp(iy_off, :);
offline_phi22_val = phi22_val_interp(iy_off, :);
offline_creamer_phi20 = phi20_creamer(iy_off, :);
offline_creamer_phi22 = phi22_creamer(iy_off, :);

dx = mean(diff(x_vec));
nx = numel(x_vec);
Lx = dx * nx;
kp = S.parameters.kp;
lambda_p = 2 * pi / kp;
x_over_lambda = x_vec / lambda_p;
focus_half_window_lambda = 4.0;
x_focus_n = S.parameters.x_focus / lambda_p;
xlim_wave = [x_focus_n - focus_half_window_lambda, x_focus_n + focus_half_window_lambda];

if mod(nx, 2) ~= 0
    error('run_directional_creamer_phi_case:EvenGridRequired', ...
        'The phi spectrum helper assumes an even number of x points.');
end

mx = [0:(nx/2), (-nx/2 + 1):-1];
kx = (2 * pi / Lx) * mx;
pos_mask = (mx >= 0);
kx_over_kp = kx(pos_mask) / kp;

centerline_phi20_val_hat = fft(centerline_phi20_val) / nx;
centerline_phi22_val_hat = fft(centerline_phi22_val) / nx;
centerline_creamer_phi20_hat = fft(centerline_creamer_phi20) / nx;
centerline_creamer_phi22_hat = fft(centerline_creamer_phi22) / nx;
offline_phi20_val_hat = fft(offline_phi20_val) / nx;
offline_phi22_val_hat = fft(offline_phi22_val) / nx;
offline_creamer_phi20_hat = fft(offline_creamer_phi20) / nx;
offline_creamer_phi22_hat = fft(offline_creamer_phi22) / nx;

spec_phi20_val = abs(centerline_phi20_val_hat(pos_mask));
spec_phi22_val = abs(centerline_phi22_val_hat(pos_mask));
spec_creamer_phi20 = abs(centerline_creamer_phi20_hat(pos_mask));
spec_creamer_phi22 = abs(centerline_creamer_phi22_hat(pos_mask));
spec_offline_phi20_val = abs(offline_phi20_val_hat(pos_mask));
spec_offline_phi22_val = abs(offline_phi22_val_hat(pos_mask));
spec_offline_creamer_phi20 = abs(offline_creamer_phi20_hat(pos_mask));
spec_offline_creamer_phi22 = abs(offline_creamer_phi22_hat(pos_mask));

tag = sprintf('alpha%d_chi%d_spread%d_Akp%04d', ...
    round(alpha_val), round(chi_deg), round(spread_deg), round(10000 * akp_val));

out_dir = fullfile(repo_root, 'MATLAB', 'output');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

fig_phi = figure('Color', 'w', 'Position', [120 120 1450 980]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(x_over_lambda, centerline_phi22_val, '-.', 'LineWidth', 1.3);
hold on;
plot(x_over_lambda, centerline_creamer_phi22, '--', 'LineWidth', 1.3);
grid on;
xlim(xlim_wave);
legend({'MF12 \phi_{22}', 'Creamer \phi_{22}'}, 'Location', 'best');
title(sprintf('kd50 MF12 \\phi_{22} comparison at y = %.6g (matched \\alpha=%g)', ...
    y_vec(iy0), validation_case.alpha));
xlabel('x / \lambda_p');
ylabel('\phi_s');

nexttile;
plot(x_over_lambda, centerline_phi20_val, ':', 'LineWidth', 1.3);
hold on;
plot(x_over_lambda, centerline_creamer_phi20, '--', 'LineWidth', 1.3);
grid on;
xlim(xlim_wave);
legend({'MF12 \phi_{20}', 'Creamer \phi_{20}'}, 'Location', 'best');
title(sprintf('kd50 MF12 \\phi_{20} comparison on centerline (matched \\alpha=%g)', ...
    validation_case.alpha));
xlabel('x / \lambda_p');
ylabel('\phi_s');

nexttile;
semilogy(kx_over_kp, spec_phi20_val + eps, ':', 'LineWidth', 1.2);
hold on;
semilogy(kx_over_kp, spec_phi22_val + eps, '-.', 'LineWidth', 1.2);
semilogy(kx_over_kp, spec_creamer_phi20 + eps, ':', 'LineWidth', 1.2, 'Color', [0.2 0.6 0.2]);
semilogy(kx_over_kp, spec_creamer_phi22 + eps, '--', 'LineWidth', 1.2, 'Color', [0.6 0.2 0.6]);
grid on;
xlim([0, 5]);
legend({'MF12 \phi_{20}', 'MF12 \phi_{22}', 'Creamer \phi_{20}', 'Creamer \phi_{22}'}, ...
    'Location', 'best');
title('Centerline spatial spectrum');
xlabel('k_x / k_p');
ylabel('Amplitude');

nexttile;
plot(x_over_lambda, offline_phi22_val, '-.', 'LineWidth', 1.3);
hold on;
plot(x_over_lambda, offline_creamer_phi22, '--', 'LineWidth', 1.3);
grid on;
xlim(xlim_wave);
legend({'MF12 \phi_{22}', 'Creamer \phi_{22}'}, 'Location', 'best');
title(sprintf('kd50 MF12 \\phi_{22} off-center at y = %.6g', y_off));
xlabel('x / \lambda_p');
ylabel('\phi_s');

nexttile;
plot(x_over_lambda, offline_phi20_val, ':', 'LineWidth', 1.3);
hold on;
plot(x_over_lambda, offline_creamer_phi20, '--', 'LineWidth', 1.3);
grid on;
xlim(xlim_wave);
legend({'MF12 \phi_{20}', 'Creamer \phi_{20}'}, 'Location', 'best');
title(sprintf('kd50 MF12 \\phi_{20} off-center at y = %.6g', y_off));
xlabel('x / \lambda_p');
ylabel('\phi_s');

nexttile;
semilogy(kx_over_kp, spec_offline_phi20_val + eps, ':', 'LineWidth', 1.2);
hold on;
semilogy(kx_over_kp, spec_offline_phi22_val + eps, '-.', 'LineWidth', 1.2);
semilogy(kx_over_kp, spec_offline_creamer_phi20 + eps, ':', 'LineWidth', 1.2, 'Color', [0.2 0.6 0.2]);
semilogy(kx_over_kp, spec_offline_creamer_phi22 + eps, '--', 'LineWidth', 1.2, 'Color', [0.6 0.2 0.6]);
grid on;
xlim([0, 5]);
legend({'MF12 \phi_{20}', 'MF12 \phi_{22}', 'Creamer \phi_{20}', 'Creamer \phi_{22}'}, ...
    'Location', 'best');
title('Off-center spatial spectrum');
xlabel('k_x / k_p');
ylabel('Amplitude');

sgtitle(sprintf(['Centerline comparison: \\phi_{22}, \\phi_{20}, and spectrum | ', ...
    'A k_p = %.3f, \\chi = %g^\\circ, \\alpha = %g, spread = %g^\\circ'], ...
    akp_val, chi_deg, alpha_val, spread_deg));

out_phi_png = fullfile(out_dir, ['directional_creamer_phi_centerline_' tag '.png']);
exportgraphics(fig_phi, out_phi_png, 'Resolution', 160);

fprintf('\nSaved phi centerline figure to:\n  %s\n', out_phi_png);
fprintf('Validation source                      = %s\n', validation_file);
fprintf('Matched alpha                          = %g\n', validation_case.alpha);
fprintf('Match rel L2                           = %.6g\n', validation_match.selected_rel_l2);
fprintf('Match amp ratio                        = %.6g\n', validation_match.selected_amp_ratio);
fprintf('Match corr                             = %.12g\n', validation_match.selected_corr);
fprintf('Creamer four-phase runtime             = %.3f s\n', creamer_separation_runtime_s);
fprintf('Off-centerline y target (half env)     = %.6g\n', y_off);
fprintf('Off-centerline env / max env           = %.6g\n', row_envelope_peak(iy_off) / max(row_envelope_peak));
fprintf('Centerline max |MF12 phi20|            = %.6g\n', max(abs(centerline_phi20_val)));
fprintf('Centerline max |MF12 phi22|            = %.6g\n', max(abs(centerline_phi22_val)));
fprintf('Centerline max |Creamer phi20|         = %.6g\n', max(abs(centerline_creamer_phi20)));
fprintf('Centerline max |Creamer phi22|         = %.6g\n', max(abs(centerline_creamer_phi22)));
fprintf('Off-center max |MF12 phi20|            = %.6g\n', max(abs(offline_phi20_val)));
fprintf('Off-center max |MF12 phi22|            = %.6g\n', max(abs(offline_phi22_val)));
fprintf('Off-center max |Creamer phi20|         = %.6g\n', max(abs(offline_creamer_phi20)));
fprintf('Off-center max |Creamer phi22|         = %.6g\n', max(abs(offline_creamer_phi22)));
