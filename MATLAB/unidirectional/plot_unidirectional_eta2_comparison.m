function out_png = plot_unidirectional_eta2_comparison(result, out_dir)
%PLOT_UNIDIRECTIONAL_ETA2_COMPARISON Save eta20/eta22 plus eta33 diagnostics.

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

x_over_lambda = result.x_vec / result.lambda_p;
k_over_kp = result.k_over_kp;
spec_mask = (k_over_kp >= 0) & (k_over_kp <= 6);

fig = figure('Color', 'w', 'Position', [120 120 1480 920]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(x_over_lambda, result.mf12.eta22, '-.', 'LineWidth', 1.35, 'Color', [0.00 0.45 0.74]);
hold on;
plot(x_over_lambda, result.creamer.eta22, '--', 'LineWidth', 1.35, 'Color', [0.85 0.33 0.10]);
grid on;
xlim(result.xlim_wave);
xlabel('x / \lambda_p');
ylabel('\eta');
title('\eta_{22}');
legend({'MF12 \eta_{22}', 'Creamer \eta_{22}'}, 'Location', 'best');

nexttile;
if isfield(result.mf12, 'eta20_filtered') && isfield(result.creamer, 'eta20_filtered')
    eta20_mf12_plot = result.mf12.eta20_filtered;
    eta20_creamer_plot = result.creamer.eta20_filtered;
    eta20_title = '\eta_{20}, 0 \leq |k|/k_p \leq 3';
else
    eta20_mf12_plot = result.mf12.eta20;
    eta20_creamer_plot = result.creamer.eta20;
    eta20_title = '\eta_{20}';
end
plot(x_over_lambda, eta20_mf12_plot, '-.', 'LineWidth', 1.35, 'Color', [0.00 0.45 0.74]);
hold on;
plot(x_over_lambda, eta20_creamer_plot, '--', 'LineWidth', 1.35, 'Color', [0.85 0.33 0.10]);
grid on;
xlim(result.xlim_wave);
xlabel('x / \lambda_p');
ylabel('\eta');
title(eta20_title);
legend({'MF12 \eta_{20}', 'Creamer \eta_{20}'}, 'Location', 'best');

nexttile;
if isfield(result.spec, 'mf12_eta20_filtered') && isfield(result.spec, 'creamer_eta20_filtered')
    eta20_mf12_spec = result.spec.mf12_eta20_filtered;
    eta20_creamer_spec = result.spec.creamer_eta20_filtered;
else
    eta20_mf12_spec = result.spec.mf12_eta20;
    eta20_creamer_spec = result.spec.creamer_eta20;
end
semilogy(k_over_kp(spec_mask), eta20_mf12_spec(spec_mask) + eps, '-.', 'LineWidth', 1.25, 'Color', [0.00 0.45 0.74]);
hold on;
semilogy(k_over_kp(spec_mask), eta20_creamer_spec(spec_mask) + eps, '--', 'LineWidth', 1.25, 'Color', [0.85 0.33 0.10]);
semilogy(k_over_kp(spec_mask), result.spec.mf12_eta22(spec_mask) + eps, '-.', 'LineWidth', 1.25, 'Color', [0.00 0.60 0.50]);
semilogy(k_over_kp(spec_mask), result.spec.creamer_eta22(spec_mask) + eps, '--', 'LineWidth', 1.25, 'Color', [0.49 0.18 0.56]);
for kv = [0.5 1 2 3 4 5 6]
    xline(kv, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
end
grid on;
xlim([0, 6]);
xlabel('k / k_p');
ylabel('Amplitude');
title('Second-order spectra, \eta_{20} filtered to 0-3k_p');
legend({'MF12 \eta_{20}, 0-3k_p', 'Creamer \eta_{20}, 0-3k_p', 'MF12 \eta_{22}', 'Creamer \eta_{22}'}, ...
    'Location', 'best');

nexttile;
plot(x_over_lambda, result.mf12.eta33, '-.', 'LineWidth', 1.35, 'Color', [0.00 0.45 0.74]);
hold on;
plot(x_over_lambda, result.creamer.eta33, '--', 'LineWidth', 1.35, 'Color', [0.85 0.33 0.10]);
grid on;
xlim(result.xlim_wave);
xlabel('x / \lambda_p');
ylabel('\eta');
title('\eta_{33}');
legend({'MF12 \eta_{33}', 'Creamer \eta_{33}'}, 'Location', 'best');

nexttile;
plot(x_over_lambda, result.creamer.eta33 - result.mf12.eta33, '--', 'LineWidth', 1.35, 'Color', [0.64 0.08 0.18]);
hold on;
yline(0, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8);
grid on;
xlim(result.xlim_wave);
xlabel('x / \lambda_p');
ylabel('\Delta \eta_{33}');
title('\eta_{33} residual');
legend({'Creamer - MF12'}, 'Location', 'best');

nexttile;
if isfield(result.spec, 'mf12_eta33') && isfield(result.spec, 'creamer_eta33')
    semilogy(k_over_kp(spec_mask), result.spec.mf12_eta33(spec_mask) + eps, '-.', 'LineWidth', 1.35, 'Color', [0.00 0.45 0.74]);
    hold on;
    semilogy(k_over_kp(spec_mask), result.spec.creamer_eta33(spec_mask) + eps, '--', 'LineWidth', 1.35, 'Color', [0.85 0.33 0.10]);
    for kv = [1 2 3 4 5 6]
        xline(kv, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
    end
    grid on;
    xlim([0, 6]);
    xlabel('k / k_p');
    ylabel('Amplitude');
    title('\eta_{33} spectrum');
    legend({'MF12 \eta_{33}', 'Creamer \eta_{33}'}, 'Location', 'best');
else
    axis off;
    text(0.1, 0.5, '\eta_{33} spectrum not available');
end

sgtitle(sprintf(['Finite-depth unidirectional Creamer vs MF12 | A k_p = %.3f, ', ...
    '\\alpha = %g, k_p h = %.3g, h = %.4g m, N_x = %d, ppw = %.2f, ', ...
    'modes = %d, min = %d, N_\\lambda = %d, model = %s, backend = %s'], ...
    result.Akp, result.Alpha, result.kp_h, result.h, numel(result.x_vec), ...
    result.points_per_wavelength, result.max_active_modes, result.min_active_modes, ...
    result.n_lambda_steps, result.lambda_flow_model, result.creamer_backend));

out_png = fullfile(out_dir, sprintf('unidirectional_finite_depth_eta2_eta33_alpha%d_Akp%04d_kph%03d_Nx%d_mc%d_nl%d_%s_%s.png', ...
    round(result.Alpha), round(10000 * result.Akp), round(100 * result.kp_h), ...
    numel(result.x_vec), round(result.max_active_modes), round(result.n_lambda_steps), ...
    result.lambda_flow_model, result.creamer_backend));
exportgraphics(fig, out_png, 'Resolution', 170);
end
