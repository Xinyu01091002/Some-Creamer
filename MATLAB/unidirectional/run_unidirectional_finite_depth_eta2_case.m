% Compare finite-depth unidirectional Creamer eta20/eta22 against live MF12.

clearvars -except Akp Alpha Nx x_vec h kph min_active_modes max_active_modes mf12_max_active_modes max_triad_active_modes triad_backend n_lambda_steps lambda_flow_model n_picard_iters energy_fraction t phishift creamer_backend;
clc;

this_dir = fileparts(mfilename('fullpath'));
matlab_dir = fileparts(this_dir);
repo_root = fileparts(matlab_dir);
addpath(this_dir);
addpath(fullfile(matlab_dir, 'core'));
addpath(fullfile(matlab_dir, 'validation'));

if ~exist('Akp', 'var') || isempty(Akp)
    Akp = 0.12;
end
if ~exist('Alpha', 'var') || isempty(Alpha)
    Alpha = 8;
end
if ~exist('Nx', 'var') || isempty(Nx)
    Nx = 4096;
end
if ~exist('max_active_modes', 'var') || isempty(max_active_modes)
    max_active_modes = 2500;
end
if ~exist('max_triad_active_modes', 'var') || isempty(max_triad_active_modes)
    max_triad_active_modes = min(60, max_active_modes);
end
if ~exist('triad_backend', 'var') || isempty(triad_backend)
    triad_backend = 'cpp';
end
if ~exist('mf12_max_active_modes', 'var') || isempty(mf12_max_active_modes)
    mf12_max_active_modes = max_active_modes;
end
if ~exist('min_active_modes', 'var') || isempty(min_active_modes)
    min_active_modes = 0;
end
if ~exist('n_lambda_steps', 'var') || isempty(n_lambda_steps)
    n_lambda_steps = 12;
end
if ~exist('lambda_flow_model', 'var') || isempty(lambda_flow_model)
    lambda_flow_model = 'canonical_pair';
end
if ~exist('creamer_backend', 'var') || isempty(creamer_backend)
    creamer_backend = 'matlab';
end
if ~exist('n_picard_iters', 'var') || isempty(n_picard_iters)
    n_picard_iters = 4;
end
if ~exist('energy_fraction', 'var') || isempty(energy_fraction)
    energy_fraction = 0.9999999999;
end
if ~exist('t', 'var') || isempty(t)
    t = 0;
end
if ~exist('phishift', 'var') || isempty(phishift)
    phishift = 0;
end

kp = 0.0279;
lambda_p = 2 * pi / kp;
if ~exist('h', 'var') || isempty(h)
    h = 2 / kp;
end
g = 9.81;

if ~exist('x_vec', 'var') || isempty(x_vec)
    x_vec = linspace(-64 * lambda_p, 64 * lambda_p, Nx);
else
    x_vec = x_vec(:).';
    Nx = numel(x_vec);
end

eta_lin = real(Linear_focus_envelope(Akp, Alpha, x_vec, t, phishift));
eta_lin = reshape(eta_lin, 1, []);

cfg = struct();
cfg.g = g;
cfg.depth_h = h;
cfg.energy_fraction = energy_fraction;
cfg.min_active_modes = min_active_modes;
cfg.max_active_modes = max_active_modes;
cfg.n_lambda_steps = n_lambda_steps;
cfg.n_picard_iters = n_picard_iters;
cfg.lambda_flow_model = lambda_flow_model;
cfg.creamer_backend = creamer_backend;
cfg.propagation_direction_deg = 0;
cfg.preserve_mean = true;
cfg.verbose = true;
cfg.cpp_exe = fullfile(repo_root, 'cpp', 'creamer_flow', 'build', 'creamer_flow_plan.exe');
cfg.cpp_job_dir = fullfile(tempdir, 'creamer_cpp_unidirectional_finite_depth_eta2');

mf12_cfg = struct();
mf12_cfg.g = g;
mf12_cfg.h = h;
mf12_cfg.t = t;
mf12_cfg.energy_fraction = energy_fraction;
mf12_cfg.max_active_modes = mf12_max_active_modes;
mf12_cfg.max_order = 3;

t_mf12 = tic;
mf12 = mf12_from_linear_focus_1d(eta_lin, x_vec, mf12_cfg);
mf12_runtime_s = toc(t_mf12);

t_creamer = tic;
creamer_sep = creamer_four_phase_separation_1d(eta_lin, x_vec, cfg);
creamer_runtime_s = toc(t_creamer);

t_triad = tic;
triad_cfg = struct('g', g, 'max_triad_active_modes', max_triad_active_modes);
if strcmpi(triad_backend, 'cpp')
    triad_cfg.cpp_exe = fullfile(repo_root, 'cpp', 'creamer_flow', 'build', 'creamer_eta33_triad_1d.exe');
    triad_cfg.cpp_job_dir = fullfile(tempdir, 'creamer_cpp_unidirectional_eta33_triad');
    creamer_triad_eta33 = creamer_eta33_h3h4_triad_1d_cpp(eta_lin, x_vec, h, triad_cfg);
else
    creamer_triad_eta33 = creamer_eta33_h3h4_triad_1d(eta_lin, x_vec, h, triad_cfg);
end
triad_runtime_s = toc(t_triad);

creamer_self_eta33 = creamer_eta33_h3h4_single_frequency_1d(eta_lin, x_vec, h);

dx = mean(diff(x_vec));
Lx = dx * Nx;
mx = [0:(Nx/2), (-Nx/2 + 1):-1];
k = (2 * pi / Lx) * mx;
k_over_kp = k / kp;
pos_mask = (mx >= 0);
points_per_wavelength = lambda_p / dx;
domain_wavelengths = Lx / lambda_p;

spec = struct();
spec.mf12_eta20 = abs(fft(mf12.eta20) / Nx);
spec.creamer_eta20 = abs(fft(creamer_sep.eta20) / Nx);
spec.mf12_eta22 = abs(fft(mf12.eta22) / Nx);
spec.creamer_eta22 = abs(fft(creamer_sep.eta22) / Nx);
spec.mf12_eta33 = abs(fft(mf12.eta33) / Nx);
spec.creamer_eta33 = abs(fft(creamer_triad_eta33.eta33_h3) / Nx);
spec.creamer_eta33_h3h4 = abs(fft(creamer_triad_eta33.eta33_h3h4) / Nx);
spec.creamer_eta33_h4_delta = abs(fft(creamer_triad_eta33.eta33_h4_delta) / Nx);
spec.creamer_eta33_flow = abs(fft(creamer_sep.eta33) / Nx);
spec.creamer_eta33_h3h4_self = abs(fft(creamer_self_eta33.eta33_h3h4) / Nx);

eta20_filter = abs(k_over_kp) <= 3;
mf12_eta20_filtered = local_filter_by_mask(mf12.eta20, eta20_filter);
creamer_eta20_filtered = local_filter_by_mask(creamer_sep.eta20, eta20_filter);
spec.mf12_eta20_filtered = abs(fft(mf12_eta20_filtered) / Nx);
spec.creamer_eta20_filtered = abs(fft(creamer_eta20_filtered) / Nx);

eta33_filter = (abs(k_over_kp) >= 2) & (abs(k_over_kp) <= 4);
mf12_eta33_filtered = local_filter_by_mask(mf12.eta33, eta33_filter);
creamer_eta33_h3_filtered = local_filter_by_mask(creamer_triad_eta33.eta33_h3, eta33_filter);
creamer_eta33_h3h4_filtered = local_filter_by_mask(creamer_triad_eta33.eta33_h3h4, eta33_filter);
spec.mf12_eta33_filtered = abs(fft(mf12_eta33_filtered) / Nx);
spec.creamer_eta33_filtered = abs(fft(creamer_eta33_h3_filtered) / Nx);
spec.creamer_eta33_h3h4_filtered = abs(fft(creamer_eta33_h3h4_filtered) / Nx);

result = struct();
result.Akp = Akp;
result.Alpha = Alpha;
result.kp = kp;
result.lambda_p = lambda_p;
result.h = h;
result.kp_h = kp * h;
result.g = g;
result.x_vec = x_vec;
result.xlim_wave = [-2, 2];
result.k_over_kp = k_over_kp(pos_mask);
result.dx = dx;
result.points_per_wavelength = points_per_wavelength;
result.domain_wavelengths = domain_wavelengths;
result.max_active_modes = max_active_modes;
result.mf12_max_active_modes = mf12_max_active_modes;
result.max_triad_active_modes = max_triad_active_modes;
result.min_active_modes = min_active_modes;
result.n_lambda_steps = n_lambda_steps;
result.lambda_flow_model = lambda_flow_model;
result.creamer_backend = creamer_backend;
result.triad_backend = triad_backend;
result.depth_model = creamer_sep.diagnostics{1}.depth_model;
result.mf12 = struct( ...
    'active_mode_count', mf12.active_mode_count, ...
    'eta20', mf12.eta20, ...
    'eta20_filtered', mf12_eta20_filtered, ...
    'eta22', mf12.eta22, ...
    'eta33', mf12.eta33, ...
    'eta33_filtered', mf12_eta33_filtered);
result.creamer = struct( ...
    'eta20', creamer_sep.eta20, ...
    'eta20_filtered', creamer_eta20_filtered, ...
    'eta22', creamer_sep.eta22, ...
    'eta22_h3', creamer_sep.eta22, ...
    'eta22_h3h4', creamer_sep.eta22, ...
    'eta33', creamer_triad_eta33.eta33_h3, ...
    'eta33_h3_triad', creamer_triad_eta33.eta33_h3, ...
    'eta33_h3h4', creamer_triad_eta33.eta33_h3h4, ...
    'eta33_filtered', creamer_eta33_h3_filtered, ...
    'eta33_h3h4_filtered', creamer_eta33_h3h4_filtered, ...
    'eta33_h3h4_triad', creamer_triad_eta33.eta33_h3h4, ...
    'eta33_h4_delta_triad', creamer_triad_eta33.eta33_h4_delta, ...
    'eta33_flow_separated', creamer_sep.eta33, ...
    'eta33_h3_self', creamer_self_eta33.eta33_h3, ...
    'eta33_h3h4_self', creamer_self_eta33.eta33_h3h4);
result.spec = struct( ...
    'mf12_eta20', spec.mf12_eta20(pos_mask), ...
    'creamer_eta20', spec.creamer_eta20(pos_mask), ...
    'mf12_eta20_filtered', spec.mf12_eta20_filtered(pos_mask), ...
    'creamer_eta20_filtered', spec.creamer_eta20_filtered(pos_mask), ...
    'mf12_eta22', spec.mf12_eta22(pos_mask), ...
    'creamer_eta22', spec.creamer_eta22(pos_mask), ...
    'mf12_eta33', spec.mf12_eta33(pos_mask), ...
    'mf12_eta33_filtered', spec.mf12_eta33_filtered(pos_mask), ...
    'creamer_eta33', spec.creamer_eta33(pos_mask), ...
    'creamer_eta33_filtered', spec.creamer_eta33_filtered(pos_mask), ...
    'creamer_eta33_h3h4', spec.creamer_eta33_h3h4(pos_mask), ...
    'creamer_eta33_h3h4_filtered', spec.creamer_eta33_h3h4_filtered(pos_mask), ...
    'creamer_eta33_h4_delta', spec.creamer_eta33_h4_delta(pos_mask), ...
    'creamer_eta33_flow', spec.creamer_eta33_flow(pos_mask), ...
    'creamer_eta33_h3h4_self', spec.creamer_eta33_h3h4_self(pos_mask));

[~, idx_3kp] = min(abs(k_over_kp(pos_mask) - 3));
mf12_eta33_3kp = result.spec.mf12_eta33(idx_3kp);
creamer_h3_eta33_3kp = result.spec.creamer_eta33(idx_3kp);
creamer_h3h4_eta33_3kp = result.spec.creamer_eta33_h3h4(idx_3kp);

out_dir = fullfile(repo_root, 'MATLAB', 'output', 'unidirectional');
out_png = plot_unidirectional_eta2_comparison(result, out_dir);

fprintf('\nSaved finite-depth second-order comparison figure to:\n  %s\n', out_png);
fprintf('MF12 runtime                            = %.3f s\n', mf12_runtime_s);
fprintf('Creamer four-phase runtime              = %.3f s\n', creamer_runtime_s);
fprintf('Creamer H3+H4 triad runtime             = %.3f s\n', triad_runtime_s);
fprintf('Depth h                                = %.12g m\n', h);
fprintf('k_p h                                  = %.12g\n', result.kp_h);
fprintf('Creamer depth model                    = %s\n', result.depth_model);
fprintf('Creamer backend                        = %s\n', creamer_backend);
fprintf('Creamer triad backend                  = %s\n', triad_backend);
fprintf('Grid spacing dx                         = %.6g m\n', dx);
fprintf('Points per wavelength                   = %.6g\n', points_per_wavelength);
fprintf('Domain length / lambda_p                = %.6g\n', domain_wavelengths);
fprintf('MF12 active modes                       = %d\n', mf12.active_mode_count);
fprintf('MF12 max active modes request           = %d\n', mf12_max_active_modes);
fprintf('Creamer active modes                    = %d\n', creamer_sep.diagnostics{1}.active_mode_count);
if isfield(creamer_triad_eta33.diagnostics, 'positive_modes')
    triad_active_count = numel(creamer_triad_eta33.diagnostics.positive_modes);
elseif isfield(creamer_triad_eta33.diagnostics, 'n_positive')
    triad_active_count = creamer_triad_eta33.diagnostics.n_positive;
else
    triad_active_count = NaN;
end
fprintf('Creamer triad active parent modes       = %g\n', triad_active_count);
fprintf('Max |MF12 eta20|                        = %.6g\n', max(abs(mf12.eta20)));
fprintf('Max |Creamer eta20|                     = %.6g\n', max(abs(creamer_sep.eta20)));
fprintf('Max |MF12 eta20|, 0-3kp filtered        = %.6g\n', max(abs(mf12_eta20_filtered)));
fprintf('Max |Creamer eta20|, 0-3kp filtered     = %.6g\n', max(abs(creamer_eta20_filtered)));
fprintf('Max |MF12 eta22|                        = %.6g\n', max(abs(mf12.eta22)));
fprintf('Max |Creamer eta22|                     = %.6g\n', max(abs(creamer_sep.eta22)));
fprintf('Max |MF12 eta33|                        = %.6g\n', max(abs(mf12.eta33)));
fprintf('Max |Creamer(H3) triad eta33|           = %.6g\n', max(abs(creamer_triad_eta33.eta33_h3)));
fprintf('Max |Creamer(H3+H4) triad eta33|        = %.6g\n', max(abs(creamer_triad_eta33.eta33_h3h4)));
fprintf('Max |Creamer flow-separated eta33|      = %.6g\n', max(abs(creamer_sep.eta33)));
fprintf('Peak-bin eta33 amplitude near 3kp:\n');
fprintf('  MF12 eta33                            = %.12g\n', mf12_eta33_3kp);
fprintf('  Creamer(H3) triad eta33               = %.12g\n', creamer_h3_eta33_3kp);
fprintf('  Creamer(H3+H4) triad eta33            = %.12g\n', creamer_h3h4_eta33_3kp);
fprintf('  Creamer flow-separated eta33 diag     = %.12g\n', result.spec.creamer_eta33_flow(idx_3kp));
fprintf('  Creamer(H3+H4) self-only diag         = %.12g\n', result.spec.creamer_eta33_h3h4_self(idx_3kp));
if isfield(creamer_triad_eta33.diagnostics, 'pruning_stats') && ~isempty(creamer_triad_eta33.diagnostics.pruning_stats)
    pstats = creamer_triad_eta33.diagnostics.pruning_stats;
    fprintf('Triad pruning diagnostics:\n');
    fprintf('  Terms seen                            = %g\n', pstats.terms_seen);
    fprintf('  Rejected by target gate               = %g\n', pstats.terms_rejected_target_gate);
    fprintf('  Derivative slots seen                 = %g\n', pstats.derivative_slots_seen);
    fprintf('  Rejected by non-linear support gate   = %g\n', pstats.derivative_slots_rejected_nonlinear_support);
    fprintf('  Normal forms seen                     = %g\n', pstats.normal_forms_seen);
    fprintf('  Zero-coefficient normal forms         = %g\n', pstats.normal_forms_zero_coeff);
    fprintf('  Resonant normal forms                 = %g\n', pstats.normal_forms_resonant);
    fprintf('  Used normal forms                     = %g\n', pstats.normal_forms_used);
    fprintf('  Targeted H4 combos tried              = %g\n', pstats.targeted_h4_combo_total);
    fprintf('  Targeted H4 consistency pass          = %g\n', pstats.targeted_h4_consistency_pass);
    fprintf('  Targeted H4 nonzero pass              = %g\n', pstats.targeted_h4_nonzero_pass);
    fprintf('  Targeted H4 unique rows               = %g\n', pstats.targeted_h4_unique_rows);
    fprintf('  Targeted H4 support-pruned combos     = %g\n', pstats.targeted_h4_support_pruned);
    fprintf('  Targeted bracket A-terms seen         = %g\n', pstats.targeted_bracket_terms_a_seen);
    fprintf('  Targeted bracket keys requested       = %g\n', pstats.targeted_bracket_needed_keys);
    fprintf('  Targeted bracket keys matched         = %g\n', pstats.targeted_bracket_matched_keys);
    fprintf('  Targeted bracket products formed      = %g\n', pstats.targeted_bracket_products_formed);
    fprintf('  Targeted bracket support-pruned       = %g\n', pstats.targeted_bracket_support_pruned);
end

function eta_filtered = local_filter_by_mask(eta, mask)
eta_hat = fft(eta);
eta_hat(~mask) = 0;
eta_filtered = real(ifft(eta_hat));
end
