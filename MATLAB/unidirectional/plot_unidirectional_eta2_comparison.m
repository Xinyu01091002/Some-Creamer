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
triad_mode_note = sprintf('MF12 active modes = %d | Creamer H3/H3+H4 parent modes = %d (same set)', ...
    local_get_field_or_default(result, {'mf12_active_mode_count', 'mf12ActiveModeCount'}, NaN), ...
    result.max_triad_active_modes);

nexttile;
plot(x_over_lambda, result.mf12.eta22, '-.', 'LineWidth', 1.35, 'Color', [0.00 0.45 0.74]);
hold on;
plot(x_over_lambda, result.creamer.eta22_h3, '--', 'LineWidth', 1.35, 'Color', [0.85 0.33 0.10]);
plot(x_over_lambda, result.creamer.eta22_h3h4, '-', 'LineWidth', 1.20, 'Color', [0.47 0.67 0.19]);
grid on;
xlim(result.xlim_wave);
xlabel('x / \lambda_p');
ylabel('\eta');
title('\eta_{22}');
legend({'MF12 \eta_{22}', 'Creamer(H3) \eta_{22}', 'Creamer(H3+H4) \eta_{22}'}, 'Location', 'best');

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
ylim(local_semilogy_limits([
    eta20_mf12_spec(spec_mask), ...
    eta20_creamer_spec(spec_mask), ...
    result.spec.mf12_eta22(spec_mask), ...
    result.spec.creamer_eta22(spec_mask)]));
xlabel('k / k_p');
ylabel('Amplitude');
title('Second-order spectra, \eta_{20} filtered to 0-3k_p');
legend({'MF12 \eta_{20}, 0-3k_p', 'Creamer \eta_{20}, 0-3k_p', 'MF12 \eta_{22}', 'Creamer \eta_{22}'}, ...
    'Location', 'best');

nexttile;
plot(x_over_lambda, result.mf12.eta33_filtered, '-.', 'LineWidth', 1.35, 'Color', [0.00 0.45 0.74]);
hold on;
plot(x_over_lambda, result.creamer.eta33_filtered, '--', 'LineWidth', 1.35, 'Color', [0.85 0.33 0.10]);
if isfield(result.creamer, 'eta33_h3h4')
    plot(x_over_lambda, result.creamer.eta33_h3h4_filtered, '-', 'LineWidth', 1.25, 'Color', [0.47 0.67 0.19]);
end
grid on;
xlim(result.xlim_wave);
xlabel('x / \lambda_p');
ylabel('\eta');
title('\eta_{33}, 2 \leq |k|/k_p \leq 4');
if isfield(result.creamer, 'eta33_h3h4')
    legend({'MF12 \eta_{33}', 'Creamer(H3) \eta_{33}', 'Creamer(H3+H4) \eta_{33}'}, 'Location', 'best');
else
    legend({'MF12 \eta_{33}', 'Creamer \eta_{33}'}, 'Location', 'best');
end

nexttile;
plot(x_over_lambda, result.creamer.eta33_filtered - result.mf12.eta33_filtered, '--', 'LineWidth', 1.35, 'Color', [0.64 0.08 0.18]);
hold on;
if isfield(result.creamer, 'eta33_h3h4')
    plot(x_over_lambda, result.creamer.eta33_h3h4_filtered - result.mf12.eta33_filtered, '-', 'LineWidth', 1.25, 'Color', [0.47 0.67 0.19]);
end
yline(0, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8);
grid on;
xlim(result.xlim_wave);
xlabel('x / \lambda_p');
ylabel('\Delta \eta_{33}');
title('\eta_{33} residual, 2 \leq |k|/k_p \leq 4');
if isfield(result.creamer, 'eta33_h3h4')
    legend({'Creamer(H3) - MF12', 'Creamer(H3+H4) - MF12'}, 'Location', 'best');
else
    legend({'Creamer - MF12'}, 'Location', 'best');
end

nexttile;
if isfield(result.spec, 'mf12_eta33') && isfield(result.spec, 'creamer_eta33')
    eta33_spec_mask = spec_mask;
    semilogy(k_over_kp(eta33_spec_mask), result.spec.mf12_eta33_filtered(eta33_spec_mask) + eps, '-.', 'LineWidth', 1.35, 'Color', [0.00 0.45 0.74]);
    hold on;
    semilogy(k_over_kp(eta33_spec_mask), result.spec.creamer_eta33_filtered(eta33_spec_mask) + eps, '--', 'LineWidth', 1.35, 'Color', [0.85 0.33 0.10]);
    if isfield(result.spec, 'creamer_eta33_h3h4')
        semilogy(k_over_kp(eta33_spec_mask), result.spec.creamer_eta33_h3h4_filtered(eta33_spec_mask) + eps, '-', 'LineWidth', 1.25, 'Color', [0.47 0.67 0.19]);
    end
    for kv = [1 2 3 4 5 6]
        xline(kv, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
    end
    grid on;
    xlim([0, 6]);
    vals_for_ylim = [result.spec.mf12_eta33_filtered(eta33_spec_mask), result.spec.creamer_eta33_filtered(eta33_spec_mask)];
    if isfield(result.spec, 'creamer_eta33_h3h4')
        vals_for_ylim = [vals_for_ylim, result.spec.creamer_eta33_h3h4_filtered(eta33_spec_mask)];
    end
    ylim(local_semilogy_limits(vals_for_ylim));
    xlabel('k / k_p');
    ylabel('Amplitude');
    title('\eta_{33} spectrum, filtered to 2 \leq |k|/k_p \leq 4');
    if isfield(result.spec, 'creamer_eta33_h3h4')
        legend({'MF12 \eta_{33}', 'Creamer(H3) \eta_{33}', 'Creamer(H3+H4) \eta_{33}'}, 'Location', 'best');
    else
        legend({'MF12 \eta_{33}', 'Creamer \eta_{33}'}, 'Location', 'best');
    end
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
annotation(fig, 'textbox', [0.54 0.905 0.42 0.032], ...
    'String', triad_mode_note, ...
    'FitBoxToText', 'off', 'EdgeColor', 'none', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'middle', 'FontSize', 9, 'Color', [0.2 0.2 0.2]);

out_png = fullfile(out_dir, sprintf(['unidirectional_finite_depth_eta2_eta33_alpha%d_Akp%04d_kph%03d_Nx%d_', ...
    'sep%d_triad%d_nl%d_%s_%s.png'], ...
    round(result.Alpha), round(10000 * result.Akp), round(100 * result.kp_h), ...
    numel(result.x_vec), round(result.max_active_modes), round(result.max_triad_active_modes), ...
    round(result.n_lambda_steps), result.lambda_flow_model, result.creamer_backend));
exportgraphics(fig, out_png, 'Resolution', 170);

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

function val = local_get_field_or_default(s, names, fallback)
val = fallback;
for i = 1:numel(names)
    if isfield(s, names{i}) && ~isempty(s.(names{i}))
        val = s.(names{i});
        return;
    end
end
if isfield(s, 'mf12') && isfield(s.mf12, 'active_mode_count')
    val = s.mf12.active_mode_count;
end
end
end
