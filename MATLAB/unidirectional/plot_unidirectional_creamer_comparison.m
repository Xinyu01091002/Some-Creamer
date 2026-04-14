function out_png = plot_unidirectional_creamer_comparison(result, out_dir)
%PLOT_UNIDIRECTIONAL_CREAMER_COMPARISON Save a 2x3 unidirectional comparison figure.

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

x_over_lambda = result.x_vec / result.lambda_p;
xlim_wave = result.xlim_wave;
k_over_kp = result.k_over_kp;
spec_mask = (k_over_kp >= 0) & (k_over_kp <= 5);

fig = figure('Color', 'w', 'Position', [120 120 1480 980]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(x_over_lambda, result.mf12.eta22, '-.', 'LineWidth', 1.35, 'Color', [0.00 0.45 0.74]);
hold on;
plot(x_over_lambda, result.creamer.eta22, '--', 'LineWidth', 1.35, 'Color', [0.85 0.33 0.10]);
grid on;
xlim(xlim_wave);
xlabel('x / \lambda_p');
ylabel('\eta');
title('\eta_{22}');
legend({'MF12 \eta_{22}', 'Creamer \eta_{22}'}, 'Location', 'best');

nexttile;
plot(x_over_lambda, result.mf12.eta20, ':', 'LineWidth', 1.55, 'Color', [0.00 0.45 0.74]);
hold on;
plot(x_over_lambda, result.creamer.eta20, '--', 'LineWidth', 1.35, 'Color', [0.85 0.33 0.10]);
grid on;
xlim(xlim_wave);
xlabel('x / \lambda_p');
ylabel('\eta');
title('\eta_{20}');
legend({'MF12 \eta_{20}', 'Creamer \eta_{20}'}, 'Location', 'best');

nexttile;
semilogy(k_over_kp(spec_mask), result.spec.mf12_eta20(spec_mask) + eps, ':', 'LineWidth', 1.35, 'Color', [0.00 0.45 0.74]);
hold on;
semilogy(k_over_kp(spec_mask), result.spec.creamer_eta20(spec_mask) + eps, '--', 'LineWidth', 1.2, 'Color', [0.85 0.33 0.10]);
semilogy(k_over_kp(spec_mask), result.spec.mf12_eta22(spec_mask) + eps, '-.', 'LineWidth', 1.35, 'Color', [0.00 0.60 0.50]);
semilogy(k_over_kp(spec_mask), result.spec.creamer_eta22(spec_mask) + eps, '--', 'LineWidth', 1.2, 'Color', [0.49 0.18 0.56]);
semilogy(k_over_kp(spec_mask), result.spec.mf12_eta33(spec_mask) + eps, '-', 'LineWidth', 1.35, 'Color', [0.47 0.67 0.19]);
semilogy(k_over_kp(spec_mask), result.spec.creamer_eta33(spec_mask) + eps, '--', 'LineWidth', 1.2, 'Color', [0.64 0.08 0.18]);
for kv = [1 2 3]
    xline(kv, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
end
grid on;
xlim([0, 5]);
xlabel('k / k_p');
ylabel('Amplitude');
title('Harmonic spectra');
legend({'MF12 \eta_{20}', 'Creamer \eta_{20}', 'MF12 \eta_{22}', 'Creamer \eta_{22}', 'MF12 \eta_{33}', 'Creamer \eta_{33}'}, ...
    'Location', 'best');

nexttile;
plot(x_over_lambda, result.mf12.eta33, '-.', 'LineWidth', 1.35, 'Color', [0.00 0.45 0.74]);
hold on;
plot(x_over_lambda, result.creamer.eta33, '--', 'LineWidth', 1.35, 'Color', [0.85 0.33 0.10]);
grid on;
xlim(xlim_wave);
xlabel('x / \lambda_p');
ylabel('\eta');
title('\eta_{33}');
legend({'MF12 \eta_{33}', 'Creamer \eta_{33}'}, 'Location', 'best');

nexttile;
plot(x_over_lambda, result.mf12.delta_total, '-.', 'LineWidth', 1.35, 'Color', [0.00 0.45 0.74]);
hold on;
plot(x_over_lambda, result.creamer.delta_total, '--', 'LineWidth', 1.35, 'Color', [0.85 0.33 0.10]);
grid on;
xlim(xlim_wave);
xlabel('x / \lambda_p');
ylabel('\Delta \eta');
title('\eta_{nl} - \eta_{lin}');
legend({'MF12 total correction', 'Creamer total correction'}, 'Location', 'best');

nexttile;
plot(x_over_lambda, result.creamer.eta11_lin, '-', 'LineWidth', 2.0, 'Color', [0.00 0.45 0.74]);
hold on;
plot(x_over_lambda, result.creamer.eta11, '--', 'LineWidth', 1.35, 'Color', [0.85 0.33 0.10]);
plot(x_over_lambda, result.creamer.eta11 - result.creamer.eta11_lin, ':', 'LineWidth', 1.2, 'Color', [0.50 0.50 0.50]);
yline(0, ':', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.6);
grid on;
xlim(xlim_wave);
xlabel('x / \lambda_p');
ylabel('\eta');
title('\eta_{11}: input vs output');
legend({'\eta_{lin} (input)', '\eta_{11} (Creamer out)', 'difference'}, 'Location', 'best');

sgtitle(sprintf(['Unidirectional Creamer vs MF12 | A k_p = %.3f, \\alpha = %g, N_x = %d, ', ...
    'k_p h = %.3g, ppw = %.2f, modes = %d, min = %d, N_\\lambda = %d, model = %s'], ...
    result.Akp, result.Alpha, numel(result.x_vec), result.kp_h, result.points_per_wavelength, ...
    result.max_active_modes, ...
    result.min_active_modes, ...
    result.n_lambda_steps, result.lambda_flow_model));

out_png = fullfile(out_dir, sprintf('unidirectional_creamer_alpha%d_Akp%04d_Nx%d_mc%d_nl%d_%s.png', ...
    round(result.Alpha), round(10000 * result.Akp), numel(result.x_vec), ...
    round(result.max_active_modes), round(result.n_lambda_steps), result.lambda_flow_model));
exportgraphics(fig, out_png, 'Resolution', 170);
end
