clearvars -except akp_index alpha_index chi_index spread_index;
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

eta_lin = S.results.eta_results{akp_index, alpha_index, chi_index, spread_index};
x_vec = S.x_vec;
y_vec = S.y_vec;
akp_val = S.parameters.Akp_values(akp_index);
alpha_val = S.parameters.Alpha_values(alpha_index);
chi_deg = S.parameters.propagation_directions(chi_index);
spread_deg = S.parameters.spread_angles(spread_index);

cfg = struct();
cfg.g = 9.81;
cfg.energy_fraction = 0.99999999;
cfg.max_active_modes = 2000;
cfg.propagation_direction_deg = chi_deg;
cfg.preserve_mean = true;
cfg.verbose = false;

[eta_nl, diagnostics] = directional_creamer_transform(eta_lin, x_vec, y_vec, cfg);

[~, iy0] = min(abs(y_vec - S.parameters.y_focus));
centerline_lin = eta_lin(iy0, :);
centerline_nl = eta_nl(iy0, :);
centerline_delta = centerline_nl - centerline_lin;

nx = numel(x_vec);
dx = mean(diff(x_vec));
Lx = dx * nx;
mx = [0:(nx/2), (-nx/2 + 1):-1];
kx = (2*pi/Lx) * mx;
pos_mask = (mx >= 0);

kp = S.parameters.kp;
kx_over_kp = kx(pos_mask) / kp;
spec_lin = abs(fft(centerline_lin) / nx);
spec_nl = abs(fft(centerline_nl) / nx);
spec_delta = abs(fft(centerline_delta) / nx);

fig = figure('Color', 'w', 'Position', [140 140 1100 650]);
semilogy(kx_over_kp, spec_lin(pos_mask) + eps, '-', 'LineWidth', 1.2);
hold on;
semilogy(kx_over_kp, spec_nl(pos_mask) + eps, '-', 'LineWidth', 1.5);
semilogy(kx_over_kp, spec_delta(pos_mask) + eps, '--', 'LineWidth', 1.2);
grid on;
xlim([0, 5]);
xline(1, ':', '1k_p', 'LabelVerticalAlignment', 'middle');
xline(2, ':', '2k_p', 'LabelVerticalAlignment', 'middle');
xline(3, ':', '3k_p', 'LabelVerticalAlignment', 'middle');
legend({'linear', 'Creamer \eta_{nl}', 'Creamer-\eta_{lin}'}, 'Location', 'best');
xlabel('k_x / k_p');
ylabel('Amplitude');
title(sprintf('Creamer centerline spectrum | Akp=%.3f, alpha=%g, chi=%g^\\circ, spread=%g^\\circ', ...
    akp_val, alpha_val, chi_deg, spread_deg));

out_dir = fullfile(repo_root, 'MATLAB', 'output');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
tag = sprintf('alpha%d_chi%d_spread%d_Akp%04d', ...
    round(alpha_val), round(chi_deg), round(spread_deg), round(10000 * akp_val));
out_png = fullfile(out_dir, ['creamer_centerline_spectrum_' tag '.png']);
exportgraphics(fig, out_png, 'Resolution', 180);

fprintf('Saved Creamer centerline spectrum to:\n  %s\n', out_png);
fprintf('Creamer active modes = %d\n', diagnostics.active_mode_count);
