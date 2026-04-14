function out_png = plot_unidirectional_eta31_profile(result, out_dir)
%PLOT_UNIDIRECTIONAL_ETA31_PROFILE Save eta31/eta33 comparison plus kp-normalized spectra.

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

x_over_lambda = result.x_vec / result.lambda_p;
xlim_wave = result.xlim_wave;
k_over_kp = result.k_over_kp;
spec_mask = (k_over_kp >= 0) & (k_over_kp <= 6);
creamer_eta31_env = abs(hilbert(result.eta31));
vwa_eta31_env = abs(hilbert(result.vwa_eta31));

fig = figure('Color', 'w', 'Position', [120 120 1420 900]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(x_over_lambda, result.eta11_input, '-', 'LineWidth', 2.0, 'Color', [0.00 0.45 0.74]);
hold on;
plot(x_over_lambda, result.eta1_total, '--', 'LineWidth', 1.35, 'Color', [0.85 0.33 0.10]);
grid on;
xlim(xlim_wave);
xlabel('x / \lambda_p');
ylabel('\eta');
title('\eta_{11} input vs four-phase \eta_1 total');
legend({'\eta_{11} input', 'four-phase \eta_1 total'}, 'Location', 'best');

nexttile;
semilogy(k_over_kp(spec_mask), result.spec.eta11_input(spec_mask) + eps, '-', ...
    'LineWidth', 2.0, 'Color', [0.00 0.45 0.74]);
hold on;
semilogy(k_over_kp(spec_mask), result.spec.eta1_total(spec_mask) + eps, '--', ...
    'LineWidth', 1.35, 'Color', [0.85 0.33 0.10]);
for kv = 1:6
    xline(kv, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
end
grid on;
xlim([0, 6]);
ylim(local_semilogy_limits([
    result.spec.eta11_input(spec_mask), ...
    result.spec.eta1_total(spec_mask)]));
xlabel('k / k_p');
ylabel('Amplitude');
title('\eta_{11} input and \eta_1 total spectra');
legend({'\eta_{11} input', 'four-phase \eta_1 total'}, 'Location', 'best');

nexttile;
third_order_xlim = [-2, 2];
plot(x_over_lambda, result.eta31, '-', 'LineWidth', 1.45, 'Color', [0.49 0.18 0.56]);
hold on;
plot(x_over_lambda, result.vwa_eta31, '--', 'LineWidth', 1.8, 'Color', [0.47 0.67 0.19]);
plot(x_over_lambda, creamer_eta31_env, '-', 'LineWidth', 1.15, 'Color', [0.74 0.58 0.84]);
plot(x_over_lambda, vwa_eta31_env, '--', 'LineWidth', 1.25, 'Color', [0.73 0.83 0.53]);
plot(x_over_lambda, result.eta33, '-.', 'LineWidth', 1.45, 'Color', [0.85 0.33 0.10]);
plot(x_over_lambda, result.vwa_eta33, ':', 'LineWidth', 1.8, 'Color', [0.00 0.45 0.74]);
yline(0, ':', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
grid on;
xlim(third_order_xlim);
xlabel('x / \lambda_p');
ylabel('\eta');
title('Third-order comparison: \eta_{31} and \eta_{33}');
legend({'Creamer H3 \eta_{31}', 'VWA \eta_{31}', ...
    'Creamer H3 env(\eta_{31})', 'VWA env(\eta_{31})', ...
    'Creamer H3 \eta_{33}', 'VWA \eta_{33}'}, ...
    'Location', 'northeast');

nexttile;
semilogy(k_over_kp(spec_mask), result.spec.eta31(spec_mask) + eps, '-', ...
    'LineWidth', 1.45, 'Color', [0.49 0.18 0.56]);
hold on;
semilogy(k_over_kp(spec_mask), result.spec.vwa_eta31(spec_mask) + eps, '--', ...
    'LineWidth', 1.8, 'Color', [0.47 0.67 0.19]);
semilogy(k_over_kp(spec_mask), result.spec.eta33(spec_mask) + eps, '-.', ...
    'LineWidth', 1.45, 'Color', [0.85 0.33 0.10]);
semilogy(k_over_kp(spec_mask), result.spec.vwa_eta33(spec_mask) + eps, ':', ...
    'LineWidth', 1.8, 'Color', [0.00 0.45 0.74]);
for kv = 1:6
    xline(kv, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
end
grid on;
xlim([0, 6]);
ylim(local_semilogy_limits([
    result.spec.eta31(spec_mask), ...
    result.spec.eta33(spec_mask), ...
    result.spec.vwa_eta33(spec_mask), ...
    result.spec.vwa_eta31(spec_mask)]));
xlabel('k / k_p');
ylabel('Amplitude');
title('Third-order spectra: \eta_{31} and \eta_{33}');
legend({'Creamer H3 \eta_{31}', 'VWA \eta_{31}', 'Creamer H3 \eta_{33}', 'VWA \eta_{33}'}, ...
    'Location', 'northeast');

sgtitle(sprintf(['Unidirectional deep-water eta31 / eta33 | H3-only four-phase | ', ...
    '1989 kernel | A k_p = %.3f, \\alpha = %g, N_x = %d, k_p h(ref) = %.3g, ', ...
    'modes = %d, min = %d, N_\\lambda = %d'], ...
    result.Akp, result.Alpha, numel(result.x_vec), result.kp_h_ref, ...
    result.max_active_modes, result.min_active_modes, result.n_lambda_steps));

out_png = fullfile(out_dir, sprintf( ...
    'unidirectional_deepwater_eta31_eta33_alpha%d_Akp%04d_Nx%d_mc%d_min%d_nl%d_khref%03d.png', ...
    round(result.Alpha), round(10000 * result.Akp), numel(result.x_vec), ...
    round(result.max_active_modes), round(result.min_active_modes), ...
    round(result.n_lambda_steps), round(10 * result.kp_h_ref)));
exportgraphics(fig, out_png, 'Resolution', 170);
end

function lim = local_semilogy_limits(vals)
vals = vals(isfinite(vals) & vals > 0);
if isempty(vals)
    lim = [1e-12, 1];
    return;
end
peak_val = max(vals);
floor_val = max(min(vals), peak_val * 1e-6);
if floor_val >= peak_val
    floor_val = peak_val / 10;
end
lim = [floor_val, 1.2 * peak_val];
end
