% Bichromatic finite-depth eta33 comparison:
%   Linear wave with exactly two cosine modes at wavenumbers k1=mk1*dk and k2=mk2*dk.
%   Compare MF12, Creamer(H3), Creamer(H3+H4) at the four sum-frequency output modes
%     3k1        = 3*mk1 * dk
%     2k1+k2     = (2*mk1+mk2) * dk
%     k1+2k2     = (mk1+2*mk2) * dk
%     3k2        = 3*mk2 * dk
%   Validates the symbolic result from Mathematica/finite_depth/20b+20c:
%     H3+H4 = MF12 (ratio = 1.000) at all sum-freq triads for kh = 0.5, 1, 2.

clearvars -except Akp1 Akp2 mk1 mk2 k1h kp n_periods Nx triad_backend;
clc;

this_dir = fileparts(mfilename('fullpath'));
matlab_dir = fileparts(this_dir);
repo_root = fileparts(matlab_dir);
addpath(this_dir);
addpath(fullfile(matlab_dir, 'core'));
addpath(fullfile(matlab_dir, 'validation'));

% --- Parameters -----------------------------------------------------------
if ~exist('Akp1', 'var') || isempty(Akp1),   Akp1 = 0.05;          end
if ~exist('Akp2', 'var') || isempty(Akp2),   Akp2 = 0.03;          end
if ~exist('mk1',  'var') || isempty(mk1),    mk1  = 1;              end
if ~exist('mk2',  'var') || isempty(mk2),    mk2  = 2;              end
if ~exist('k1h',  'var') || isempty(k1h),    k1h  = 1.0;            end
if ~exist('kp',   'var') || isempty(kp),     kp   = 0.0279;         end
if ~exist('n_periods', 'var') || isempty(n_periods), n_periods = 64; end
if ~exist('Nx',   'var') || isempty(Nx),     Nx   = 4096;           end
if ~exist('triad_backend', 'var') || isempty(triad_backend)
    triad_backend = 'matlab';
end

g = 9.81;

% Grid: dk = kp, so k1 = mk1*kp, k2 = mk2*kp
% Lx = n_periods * 2*pi/kp  (n_periods wavelengths at the base wavenumber kp)
dk    = kp;                        % base wavenumber = grid spacing
k1    = mk1 * dk;
k2    = mk2 * dk;
h     = k1h / k1;                  % depth: k1*h = k1h
k2h   = k2 * h;
Lx    = n_periods * 2*pi / kp;
x_vec = (0:Nx-1) * (Lx / Nx);
dx    = Lx / Nx;

% Physical amplitudes
A1 = Akp1 / k1;
A2 = Akp2 / k2;

% Bichromatic linear surface
eta_lin = A1 * cos(k1 * x_vec) + A2 * cos(k2 * x_vec);

fprintf('=== Bichromatic eta33 test ===\n');
fprintf('  mk1=%d, mk2=%d, k1h=%.4g, k2h=%.4g\n', mk1, mk2, k1h, k2h);
fprintf('  A1=%.6g (Akp1=%.4g), A2=%.6g (Akp2=%.4g)\n', A1, Akp1, A2, Akp2);
fprintf('  Nx=%d, n_periods=%d, ppw_k1=%.2f\n', Nx, n_periods, 2*pi/k1/dx);

% --- MF12 reference -------------------------------------------------------
mf12_cfg = struct('g', g, 'h', h, 't', 0, ...
    'energy_fraction', 1, 'max_active_modes', 2, 'max_order', 3);
t0 = tic;
mf12 = mf12_from_linear_focus_1d(eta_lin, x_vec, mf12_cfg);
fprintf('MF12 runtime: %.3f s  (active modes: %d)\n', toc(t0), mf12.active_mode_count);

% --- Creamer H3/H3+H4 triad -----------------------------------------------
% NOTE: 'targeted' evaluation model has a known bug (overcounts H4 correction).
%       Use 'dictionary' for correctness.  The targeted mode is faster for large
%       mode sets but is only calibrated for single-frequency self-interactions.
triad_cfg = struct('g', g, 'max_triad_active_modes', 2, 'evaluation_model', 'dictionary');
if strcmpi(triad_backend, 'cpp')
    triad_cfg.cpp_exe = fullfile(repo_root, 'cpp', 'creamer_flow', 'build', ...
        'creamer_eta33_triad_1d.exe');
    triad_cfg.cpp_job_dir = fullfile(tempdir, 'creamer_bichromatic_triad');
    t0 = tic;
    creamer = creamer_eta33_h3h4_triad_1d_cpp(eta_lin, x_vec, h, triad_cfg);
else
    t0 = tic;
    creamer = creamer_eta33_h3h4_triad_1d(eta_lin, x_vec, h, triad_cfg);
end
fprintf('Creamer triad runtime: %.3f s  (backend: %s)\n', toc(t0), triad_backend);

% --- Spectral extraction --------------------------------------------------
mx      = [0:(Nx/2), (-Nx/2 + 1):-1];
k_grid  = (2*pi/Lx) * mx;
k_over_dk = k_grid / dk;

% Expected output mode indices (integer multiples of dk)
output_modes = [3*mk1, 2*mk1+mk2, mk1+2*mk2, 3*mk2];
output_labels = {sprintf('3k1  (m=%d)',   3*mk1), ...
                 sprintf('2k1+k2 (m=%d)', 2*mk1+mk2), ...
                 sprintf('k1+2k2 (m=%d)', mk1+2*mk2), ...
                 sprintf('3k2  (m=%d)',   3*mk2)};

% FFT spectra
spec_mf12  = abs(fft(mf12.eta33)      / Nx);
spec_h3    = abs(creamer.eta33_h3_hat);
spec_h3h4  = abs(creamer.eta33_h3h4_hat);

fprintf('\n--- Spectral amplitudes at sum-freq output bins ---\n');
fprintf('%-18s  %12s  %12s  %12s  %8s  %8s\n', ...
    'Mode', 'MF12', 'H3-only', 'H3+H4', 'H3/MF12', 'H3+H4/MF12');
fprintf('%s\n', repmat('-', 1, 76));

amps = struct();
for j = 1:numel(output_modes)
    m_target = output_modes(j);
    [~, idx] = min(abs(k_over_dk - m_target));
    if abs(k_over_dk(idx) - m_target) > 0.5
        fprintf('  WARNING: exact bin for m=%d not found in spectrum\n', m_target);
    end
    a_mf12  = spec_mf12(idx);
    a_h3    = spec_h3(idx);
    a_h3h4  = spec_h3h4(idx);
    amps(j).mf12   = a_mf12;
    amps(j).h3     = a_h3;
    amps(j).h3h4   = a_h3h4;
    amps(j).label  = output_labels{j};
    ratio_h3    = local_safe_ratio(a_h3,   a_mf12);
    ratio_h3h4  = local_safe_ratio(a_h3h4, a_mf12);
    fprintf('%-18s  %12.6g  %12.6g  %12.6g  %8.5f  %8.5f\n', ...
        output_labels{j}, a_mf12, a_h3, a_h3h4, ratio_h3, ratio_h3h4);
end
fprintf('\n');
fprintf('Symbolic prediction (Mathematica 20b+20c): H3+H4/MF12 = 1.000 for all triads.\n');

% --- Plot -----------------------------------------------------------------
out_dir = fullfile(repo_root, 'MATLAB', 'output', 'unidirectional');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

pos_mask    = mx >= 0;
k_plot      = k_grid(pos_mask) / k1;    % k/k1
spec_mask   = (k_plot >= 0) & (k_plot <= 7);

fig = figure('Color', 'w', 'Position', [120 120 1400 900]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

x_over_lambda1 = x_vec * k1 / (2*pi);
xlim_wave = [0, 4];

% Panel 1: eta_lin
nexttile;
plot(x_over_lambda1, eta_lin, '-', 'LineWidth', 1.2, 'Color', [0.4 0.4 0.4]);
grid on; xlim(xlim_wave);
xlabel('x / \lambda_1'); ylabel('\eta_{lin}');
title('Linear surface');

% Panel 2: eta33 time domain
nexttile;
plot(x_over_lambda1, mf12.eta33,         '-',  'LineWidth', 2.2, 'Color', [0.00 0.45 0.74]);
hold on;
plot(x_over_lambda1, creamer.eta33_h3,   '--', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
plot(x_over_lambda1, creamer.eta33_h3h4, ':',  'LineWidth', 2.0, 'Color', [0.47 0.67 0.19]);
grid on; xlim(xlim_wave);
xlabel('x / \lambda_1'); ylabel('\eta_{33}');
title('\eta_{33} profiles');
legend({'MF12', 'Creamer(H3)', 'Creamer(H3+H4)'}, 'Location', 'best');

% Panel 3: residuals
nexttile;
plot(x_over_lambda1, creamer.eta33_h3   - mf12.eta33, '--', 'LineWidth', 1.4, 'Color', [0.85 0.33 0.10]);
hold on;
plot(x_over_lambda1, creamer.eta33_h3h4 - mf12.eta33, '-',  'LineWidth', 1.4, 'Color', [0.47 0.67 0.19]);
yline(0, ':', 'Color', [0.4 0.4 0.4]);
grid on; xlim(xlim_wave);
xlabel('x / \lambda_1'); ylabel('\Delta\eta_{33}');
title('Residual vs MF12');
legend({'Creamer(H3)-MF12', 'Creamer(H3+H4)-MF12'}, 'Location', 'best');

% Panel 4: spectrum (eta33)
nexttile;
mf12_pos  = spec_mf12(pos_mask);
h3_pos    = spec_h3(pos_mask);
h3h4_pos  = spec_h3h4(pos_mask);
semilogy(k_plot(spec_mask), mf12_pos(spec_mask)  + eps, '-',  'LineWidth', 2.2, 'Color', [0.00 0.45 0.74]);
hold on;
semilogy(k_plot(spec_mask), h3_pos(spec_mask)   + eps, '--', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
semilogy(k_plot(spec_mask), h3h4_pos(spec_mask) + eps, ':',  'LineWidth', 2.0, 'Color', [0.47 0.67 0.19]);
for mv = output_modes
    xline(mv/mk1, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
end
grid on; xlim([0, 7]);
ylim(local_semilogy_limits([mf12_pos(spec_mask), h3_pos(spec_mask), h3h4_pos(spec_mask)]));
xlabel('k / k_1'); ylabel('Amplitude');
title('\eta_{33} spectrum');
legend({'MF12', 'Creamer(H3)', 'Creamer(H3+H4)'}, 'Location', 'best');

% Panel 5: bar chart at output modes
nexttile;
n_out  = numel(output_modes);
bar_data = zeros(n_out, 3);
for j = 1:n_out
    bar_data(j, :) = [amps(j).mf12, amps(j).h3, amps(j).h3h4];
end
b = bar(bar_data);
b(1).FaceColor = [0.00 0.45 0.74];
b(2).FaceColor = [0.85 0.33 0.10];
b(3).FaceColor = [0.47 0.67 0.19];
xticks(1:n_out);
xticklabels(output_labels);
xtickangle(20);
grid on;
ylabel('Spectral amplitude');
title('Output bin amplitudes');
legend({'MF12', 'H3', 'H3+H4'}, 'Location', 'best');

% Panel 6: H3+H4 / MF12 ratio per bin
nexttile;
ratio_h3h4_vals = zeros(1, n_out);
ratio_h3_vals   = zeros(1, n_out);
for j = 1:n_out
    ratio_h3h4_vals(j) = local_safe_ratio(amps(j).h3h4, amps(j).mf12);
    ratio_h3_vals(j)   = local_safe_ratio(amps(j).h3,   amps(j).mf12);
end
br = bar([ratio_h3_vals; ratio_h3h4_vals].');
br(1).FaceColor = [0.85 0.33 0.10];
br(2).FaceColor = [0.47 0.67 0.19];
yline(1, '-k', 'LineWidth', 1.5);
xticks(1:n_out);
xticklabels(output_labels);
xtickangle(20);
grid on; ylim([0, 2]);
ylabel('ratio / MF12');
title('H3-only / H3+H4 vs MF12  (target = 1)');
legend({'H3-only', 'H3+H4'}, 'Location', 'best');

sgtitle(sprintf(['Bichromatic eta33 | k_1h=%.3g, k_2h=%.3g, mk1=%d, mk2=%d, ', ...
    'Ak_1=%.4g, Ak_2=%.4g, Nx=%d, backend=%s'], ...
    k1h, k2h, mk1, mk2, Akp1, Akp2, Nx, triad_backend));

out_png = fullfile(out_dir, sprintf( ...
    'unidirectional_bichromatic_eta33_mk%d_mk%d_k1h%03d_Akp1_%04d_Akp2_%04d_Nx%d_%s.png', ...
    mk1, mk2, round(100*k1h), round(10000*Akp1), round(10000*Akp2), Nx, triad_backend));
exportgraphics(fig, out_png, 'Resolution', 170);
fprintf('Saved figure to:\n  %s\n', out_png);

% -------------------------------------------------------------------------
function r = local_safe_ratio(a, b)
if abs(b) > 1e-15
    r = a / b;
else
    r = NaN;
end
end

function lim = local_semilogy_limits(vals)
vals = vals(isfinite(vals) & vals > 0);
if isempty(vals)
    lim = [1e-12, 1];
    return;
end
peak_val  = max(vals);
floor_val = max(min(vals), peak_val * 1e-4);
if floor_val >= peak_val, floor_val = peak_val / 10; end
lim = [floor_val, 1.2 * peak_val];
end
