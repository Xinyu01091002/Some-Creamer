% First-pass driver for the directional Creamer MATLAB prototype.
%
% Locked v1 case:
%   Akp   = 0.12
%   alpha = 8
%   chi   = 0 deg
%   spread= 30 deg
%
% Data source:
%   directional_waves_multi_Akp0.020-0.180_Alpha1-8_Chi0-0_Spread5-30_adaptive100pc_20260409_154705.mat

clearvars -except akp_index alpha_index chi_index spread_index;
clc;

this_dir = fileparts(mfilename('fullpath'));
matlab_dir = fileparts(this_dir);
repo_root = fileparts(matlab_dir);
addpath(this_dir);
addpath(fullfile(matlab_dir, 'core'));
addpath(fullfile(matlab_dir, 'validation'));
data_file = fullfile(repo_root, ...
    'directional_waves_multi_Akp0.020-0.180_Alpha1-8_Chi0-0_Spread5-30_adaptive100pc_20260409_154705.mat');

if ~isfile(data_file)
    error('run_directional_creamer_case:MissingData', ...
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
    spread_index = 6;  % spread = 30 deg
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
cfg.propagation_direction_deg = S.parameters.propagation_directions(chi_index);
cfg.preserve_mean = true;
cfg.verbose = true;

t_creamer = tic;
[eta_nl, diagnostics] = directional_creamer_transform(eta_lin, x_vec, y_vec, cfg);
creamer_runtime_s = toc(t_creamer);

validation_file = fullfile(repo_root, ...
    'directional_validation_results_mf12_linear_groups_linear_groups_kd50_mc1500_fixakp_20260409_190137.mat');
if ~isfile(validation_file)
    error('run_directional_creamer_case:MissingValidationData', ...
        'Could not find validation data file: %s', validation_file);
end

[validation_case, validation_match] = match_validation_second_order_case( ...
    validation_file, S.parameters.Alpha_values, chi_deg, spread_deg, 50, ...
    eta_lin, x_vec, y_vec, akp_val);

out_dir = fullfile(repo_root, 'MATLAB', 'output', 'directional');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

[~, ix0] = min(abs(x_vec - S.parameters.x_focus));
[~, iy0] = min(abs(y_vec - S.parameters.y_focus));

eta_lin_analytic_x = hilbert(eta_lin.').';
row_envelope_peak = max(abs(eta_lin_analytic_x), [], 2);
target_off_centerline = 0.5 * max(row_envelope_peak);
candidate_mask = (y_vec(:) >= S.parameters.y_focus) & ((1:numel(y_vec)).' ~= iy0);
candidate_idx = find(candidate_mask);
[~, best_off_idx_local] = min(abs(row_envelope_peak(candidate_idx) - target_off_centerline));
iy_off = candidate_idx(best_off_idx_local);
y_off = y_vec(iy_off);

centerline_lin = eta_lin(iy0, :);
centerline_nl = eta_nl(iy0, :);
centerline_delta = centerline_nl - centerline_lin;

[X_target, Y_target] = meshgrid(x_vec, y_vec);
[X_val, Y_val] = meshgrid(validation_case.x(:).', validation_case.y(:));
eta20_val_interp = interp2(X_val, Y_val, validation_case.eta20, X_target, Y_target, 'linear', 0);
eta22_val_interp = interp2(X_val, Y_val, validation_case.eta22, X_target, Y_target, 'linear', 0);
eta2_total_val_interp = interp2(X_val, Y_val, validation_case.eta2_total, X_target, Y_target, 'linear', 0);

t_creamer_sep = tic;
creamer_sep = creamer_four_phase_separation(eta_lin, x_vec, y_vec, cfg);
creamer_separation_runtime_s = toc(t_creamer_sep);
eta20_creamer = creamer_sep.eta20;
eta22_creamer = creamer_sep.eta22;

centerline_eta20_val = eta20_val_interp(iy0, :);
centerline_eta22_val = eta22_val_interp(iy0, :);
centerline_eta2_total_val = eta2_total_val_interp(iy0, :);
centerline_creamer_eta20 = eta20_creamer(iy0, :);
centerline_creamer_eta22 = eta22_creamer(iy0, :);

offline_lin = eta_lin(iy_off, :);
offline_nl = eta_nl(iy_off, :);
offline_delta = offline_nl - offline_lin;
offline_eta20_val = eta20_val_interp(iy_off, :);
offline_eta22_val = eta22_val_interp(iy_off, :);
offline_eta2_total_val = eta2_total_val_interp(iy_off, :);
offline_creamer_eta20 = eta20_creamer(iy_off, :);
offline_creamer_eta22 = eta22_creamer(iy_off, :);

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
    error('run_directional_creamer_case:EvenGridRequired', ...
        'The centerline spectrum helper assumes an even number of x points.');
end

mx = [0:(nx/2), (-nx/2 + 1):-1];
kx = (2 * pi / Lx) * mx;

centerline_lin_hat = fft(centerline_lin) / nx;
centerline_nl_hat = fft(centerline_nl) / nx;
centerline_delta_hat = fft(centerline_delta) / nx;
centerline_eta20_val_hat = fft(centerline_eta20_val) / nx;
centerline_eta22_val_hat = fft(centerline_eta22_val) / nx;
centerline_eta2_total_val_hat = fft(centerline_eta2_total_val) / nx;
centerline_creamer_eta20_hat = fft(centerline_creamer_eta20) / nx;
centerline_creamer_eta22_hat = fft(centerline_creamer_eta22) / nx;
offline_lin_hat = fft(offline_lin) / nx;
offline_nl_hat = fft(offline_nl) / nx;
offline_delta_hat = fft(offline_delta) / nx;
offline_eta20_val_hat = fft(offline_eta20_val) / nx;
offline_eta22_val_hat = fft(offline_eta22_val) / nx;
offline_eta2_total_val_hat = fft(offline_eta2_total_val) / nx;
offline_creamer_eta20_hat = fft(offline_creamer_eta20) / nx;
offline_creamer_eta22_hat = fft(offline_creamer_eta22) / nx;

pos_mask = (mx >= 0);
kx_pos = kx(pos_mask);
kx_over_kp = kx_pos / kp;
spec_lin = abs(centerline_lin_hat(pos_mask));
spec_nl = abs(centerline_nl_hat(pos_mask));
spec_delta = abs(centerline_delta_hat(pos_mask));
spec_eta20_val = abs(centerline_eta20_val_hat(pos_mask));
spec_eta22_val = abs(centerline_eta22_val_hat(pos_mask));
spec_eta2_total_val = abs(centerline_eta2_total_val_hat(pos_mask));
spec_creamer_eta20 = abs(centerline_creamer_eta20_hat(pos_mask));
spec_creamer_eta22 = abs(centerline_creamer_eta22_hat(pos_mask));
spec_offline_lin = abs(offline_lin_hat(pos_mask));
spec_offline_nl = abs(offline_nl_hat(pos_mask));
spec_offline_delta = abs(offline_delta_hat(pos_mask));
spec_offline_eta20_val = abs(offline_eta20_val_hat(pos_mask));
spec_offline_eta22_val = abs(offline_eta22_val_hat(pos_mask));
spec_offline_eta2_total_val = abs(offline_eta2_total_val_hat(pos_mask));
spec_offline_creamer_eta20 = abs(offline_creamer_eta20_hat(pos_mask));
spec_offline_creamer_eta22 = abs(offline_creamer_eta22_hat(pos_mask));

tag = sprintf('alpha%d_chi%d_spread%d_Akp%04d', ...
    round(alpha_val), round(chi_deg), round(spread_deg), round(10000 * akp_val));

out_mat = fullfile(out_dir, ['directional_creamer_' tag '.mat']);

fig = figure('Color', 'w', 'Position', [100 100 1500 900]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
imagesc(x_vec, y_vec, eta_lin);
axis image;
set(gca, 'YDir', 'normal');
title('Linear parent surface \eta_{lin}(x,y)');
xlabel('x');
ylabel('y');
colorbar;

nexttile;
imagesc(x_vec, y_vec, eta_nl);
axis image;
set(gca, 'YDir', 'normal');
title('Directional Creamer surface \eta_{nl}(x,y)');
xlabel('x');
ylabel('y');
colorbar;

nexttile;
plot(x_vec, eta_lin(iy0, :), 'LineWidth', 1.1);
hold on;
plot(x_vec, eta_nl(iy0, :), 'LineWidth', 1.1);
grid on;
legend({'linear', 'Creamer'}, 'Location', 'best');
title(sprintf('Centerline at y = %.6g', y_vec(iy0)));
xlabel('x');
ylabel('\eta');

nexttile;
plot(y_vec, eta_lin(:, ix0), 'LineWidth', 1.1);
hold on;
plot(y_vec, eta_nl(:, ix0), 'LineWidth', 1.1);
grid on;
legend({'linear', 'Creamer'}, 'Location', 'best');
title(sprintf('Crossline at x = %.6g', x_vec(ix0)));
xlabel('y');
ylabel('\eta');

sgtitle(sprintf('Directional Creamer prototype: alpha=%g, chi=%g, spread=%g, Akp=%.3f', ...
    alpha_val, chi_deg, spread_deg, akp_val));

out_png = fullfile(out_dir, ['directional_creamer_' tag '.png']);
exportgraphics(fig, out_png, 'Resolution', 160);

fig_centerline = figure('Color', 'w', 'Position', [120 120 1450 980]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(x_over_lambda, centerline_eta22_val, '-.', 'LineWidth', 1.3);
hold on;
plot(x_over_lambda, centerline_creamer_eta22, '--', 'LineWidth', 1.3);
grid on;
xlim(xlim_wave);
legend({'MF12 \eta_{22}', 'Creamer \eta_{22}'}, 'Location', 'best');
title(sprintf('kd50 MF12 \\eta_{22} comparison at y = %.6g (matched \\alpha=%g)', ...
    y_vec(iy0), validation_case.alpha));
xlabel('x / \lambda_p');
ylabel('\eta');

nexttile;
plot(x_over_lambda, centerline_eta20_val, ':', 'LineWidth', 1.3);
hold on;
plot(x_over_lambda, centerline_creamer_eta20, '--', 'LineWidth', 1.3);
grid on;
xlim(xlim_wave);
second_order_legend = {'MF12 \eta_{20}', 'Creamer \eta_{20}'};
legend(second_order_legend, 'Location', 'best');
title(sprintf('kd50 MF12 \\eta_{20} comparison on centerline (matched \\alpha=%g)', ...
    validation_case.alpha));
xlabel('x / \lambda_p');
ylabel('\eta');

nexttile;
semilogy(kx_over_kp, spec_delta + eps, '--', 'LineWidth', 1.1);
hold on;
semilogy(kx_over_kp, spec_eta20_val + eps, ':', 'LineWidth', 1.2);
semilogy(kx_over_kp, spec_eta22_val + eps, '-.', 'LineWidth', 1.2);
semilogy(kx_over_kp, spec_creamer_eta20 + eps, ':', 'LineWidth', 1.2, 'Color', [0.2 0.6 0.2]);
semilogy(kx_over_kp, spec_creamer_eta22 + eps, '-.', 'LineWidth', 1.2, 'Color', [0.6 0.2 0.6]);
grid on;
xlim([0, 5]);
legend({'Creamer total 2nd', 'MF12 \eta_{20}', 'MF12 \eta_{22}', ...
    'Creamer \eta_{20}', 'Creamer \eta_{22}'}, ...
    'Location', 'best');
title('Centerline spatial spectrum');
xlabel('k_x / k_p');
ylabel('Amplitude');

nexttile;
plot(x_over_lambda, offline_eta22_val, '-.', 'LineWidth', 1.3);
hold on;
plot(x_over_lambda, offline_creamer_eta22, '--', 'LineWidth', 1.3);
grid on;
xlim(xlim_wave);
legend({'MF12 \eta_{22}', 'Creamer \eta_{22}'}, 'Location', 'best');
title(sprintf('kd50 MF12 \\eta_{22} off-center at y = %.6g', y_off));
xlabel('x / \lambda_p');
ylabel('\eta');

nexttile;
plot(x_over_lambda, offline_eta20_val, ':', 'LineWidth', 1.3);
hold on;
plot(x_over_lambda, offline_creamer_eta20, '--', 'LineWidth', 1.3);
grid on;
xlim(xlim_wave);
legend({'MF12 \eta_{20}', 'Creamer \eta_{20}'}, 'Location', 'best');
title(sprintf('kd50 MF12 \\eta_{20} off-center at y = %.6g', y_off));
xlabel('x / \lambda_p');
ylabel('\eta');

nexttile;
semilogy(kx_over_kp, spec_offline_delta + eps, '--', 'LineWidth', 1.1);
hold on;
semilogy(kx_over_kp, spec_offline_eta20_val + eps, ':', 'LineWidth', 1.2);
semilogy(kx_over_kp, spec_offline_eta22_val + eps, '-.', 'LineWidth', 1.2);
semilogy(kx_over_kp, spec_offline_creamer_eta20 + eps, ':', 'LineWidth', 1.2, 'Color', [0.2 0.6 0.2]);
semilogy(kx_over_kp, spec_offline_creamer_eta22 + eps, '-.', 'LineWidth', 1.2, 'Color', [0.6 0.2 0.6]);
grid on;
xlim([0, 5]);
legend({'Creamer total 2nd', 'MF12 \eta_{20}', 'MF12 \eta_{22}', ...
    'Creamer \eta_{20}', 'Creamer \eta_{22}'}, ...
    'Location', 'best');
title('Off-center spatial spectrum');
xlabel('k_x / k_p');
ylabel('Amplitude');

sgtitle(sprintf(['Centerline comparison: \\eta_{22}, \\eta_{20}, and spectrum | ', ...
    'A k_p = %.3f, \\chi = %g^\\circ, \\alpha = %g, spread = %g^\\circ'], ...
    akp_val, S.parameters.propagation_directions(chi_index), ...
    S.parameters.Alpha_values(alpha_index), S.parameters.spread_angles(spread_index)));

out_centerline_png = fullfile(out_dir, 'directional_creamer_centerline_alpha8_chi0_spread5_Akp0020.png');
out_centerline_png = fullfile(out_dir, ['directional_creamer_centerline_' tag '.png']);
exportgraphics(fig_centerline, out_centerline_png, 'Resolution', 160);

[~, peak_idx_lin] = max(spec_lin(2:end));
[~, peak_idx_nl] = max(spec_nl(2:end));
peak_idx_lin = peak_idx_lin + 1;
peak_idx_nl = peak_idx_nl + 1;

fprintf('\nSkipped MATLAB .mat export for this run.\n');
fprintf('Figure saved to:\n  %s\n', out_png);
fprintf('Saved centerline figure to:\n  %s\n', out_centerline_png);
fprintf('Peak amplitude ratio max(|eta_nl|)/max(|eta_lin|) = %.6f\n', ...
    diagnostics.max_eta_nl / diagnostics.max_eta_lin);
fprintf('Creamer runtime (current prototype)    = %.3f s\n', creamer_runtime_s);
fprintf('Creamer four-phase separation runtime  = %.3f s\n', creamer_separation_runtime_s);
fprintf('Creamer active modes / cap             = %d / %d\n', diagnostics.active_mode_count, cfg.max_active_modes);
fprintf('Validation MF12 source                 = %s\n', validation_file);
fprintf('Validation kd target                   = %.6g\n', validation_case.kd);
fprintf('Validation matched alpha              = %.6g\n', validation_case.alpha);
fprintf('Validation match rel L2               = %.6g\n', validation_match.selected_rel_l2);
fprintf('Validation match amp ratio            = %.6g\n', validation_match.selected_amp_ratio);
fprintf('Validation match corr                 = %.12f\n', validation_match.selected_corr);
fprintf('Off-centerline y target (half env)    = %.6g\n', y_off);
fprintf('Off-centerline env / max env          = %.6g\n', row_envelope_peak(iy_off) / max(row_envelope_peak));
fprintf('Centerline peak wavenumber (linear)    = %.6f (kx/kp = %.6f)\n', ...
    kx_pos(peak_idx_lin), kx_over_kp(peak_idx_lin));
fprintf('Centerline peak wavenumber (Creamer)   = %.6f (kx/kp = %.6f)\n', ...
    kx_pos(peak_idx_nl), kx_over_kp(peak_idx_nl));
fprintf('Centerline max |delta eta|             = %.6g\n', max(abs(centerline_delta)));
fprintf('Centerline max |MF12 eta20|            = %.6g\n', max(abs(centerline_eta20_val)));
fprintf('Centerline max |MF12 eta22|            = %.6g\n', max(abs(centerline_eta22_val)));
fprintf('Centerline max |MF12 eta2_total|       = %.6g\n', max(abs(centerline_eta2_total_val)));
fprintf('Centerline max |Creamer eta20|         = %.6g\n', max(abs(centerline_creamer_eta20)));
fprintf('Centerline max |Creamer eta22|         = %.6g\n', max(abs(centerline_creamer_eta22)));
fprintf('Off-center max |delta eta|             = %.6g\n', max(abs(offline_delta)));
fprintf('Off-center max |MF12 eta20|            = %.6g\n', max(abs(offline_eta20_val)));
fprintf('Off-center max |MF12 eta22|            = %.6g\n', max(abs(offline_eta22_val)));
fprintf('Off-center max |MF12 eta2_total|       = %.6g\n', max(abs(offline_eta2_total_val)));
fprintf('Off-center max |Creamer eta20|         = %.6g\n', max(abs(offline_creamer_eta20)));
fprintf('Off-center max |Creamer eta22|         = %.6g\n', max(abs(offline_creamer_eta22)));
