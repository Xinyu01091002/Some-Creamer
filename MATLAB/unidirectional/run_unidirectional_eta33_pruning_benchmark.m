% Focused 1D finite-depth eta33 pruning benchmark against MF12.

clearvars -except Akp Alpha Nx h kph max_triad_active_modes evaluation_modes t phishift energy_fraction;
clc;

this_dir = fileparts(mfilename('fullpath'));
matlab_dir = fileparts(this_dir);
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
    Nx = 768;
end
if ~exist('max_triad_active_modes', 'var') || isempty(max_triad_active_modes)
    max_triad_active_modes = 6;
end
if ~exist('evaluation_modes', 'var') || isempty(evaluation_modes)
    evaluation_modes = {'targeted', 'semi_pruned', 'less_pruned'};
end
if ~exist('t', 'var') || isempty(t)
    t = 0;
end
if ~exist('phishift', 'var') || isempty(phishift)
    phishift = 0;
end
if ~exist('energy_fraction', 'var') || isempty(energy_fraction)
    energy_fraction = 0.9999999999;
end

kp = 0.0279;
lambda_p = 2 * pi / kp;
if ~exist('h', 'var') || isempty(h)
    if ~exist('kph', 'var') || isempty(kph)
        kph = 0.5;
    end
    h = kph / kp;
else
    kph = kp * h;
end

x_vec = linspace(-64 * lambda_p, 64 * lambda_p, Nx);
eta_lin = real(Linear_focus_envelope(Akp, Alpha, x_vec, t, phishift));
eta_lin = reshape(eta_lin, 1, []);

mf12_cfg = struct();
mf12_cfg.g = 9.81;
mf12_cfg.h = h;
mf12_cfg.t = t;
mf12_cfg.energy_fraction = energy_fraction;
mf12_cfg.max_active_modes = 2500;
mf12_cfg.max_order = 3;

t_mf12 = tic;
mf12 = mf12_from_linear_focus_1d(eta_lin, x_vec, mf12_cfg);
mf12_runtime_s = toc(t_mf12);

dx = mean(diff(x_vec));
Lx = dx * Nx;
mx = [0:(Nx/2), (-Nx/2 + 1):-1];
k = (2*pi/Lx) * mx;
k_over_kp = k / kp;
eta33_filter = (abs(k_over_kp) >= 2) & (abs(k_over_kp) <= 4);
mf12_eta33_filtered = local_filter_by_mask(mf12.eta33, eta33_filter);
mf12_eta33_spec = abs(fft(mf12.eta33) / Nx);
[~, idx_3kp] = min(abs(k_over_kp - 3));
mf12_3kp = mf12_eta33_spec(idx_3kp);
mf12_peak = max(abs(mf12_eta33_filtered));

rows = repmat(struct( ...
    'mode', "", ...
    'support_kind', "", ...
    'runtime_s', NaN, ...
    'filtered_peak', NaN, ...
    'peak_rel_to_mf12', NaN, ...
    'amp_3kp', NaN, ...
    'amp_3kp_rel_to_mf12', NaN, ...
    'candidate_terms', NaN, ...
    'used_terms', NaN, ...
    'h4_unique_rows', NaN, ...
    'h4_support_pruned', NaN, ...
    'bracket_products', NaN, ...
    'bracket_support_pruned', NaN), 1, numel(evaluation_modes));

for i = 1:numel(evaluation_modes)
    triad_cfg = struct( ...
        'g', 9.81, ...
        'max_triad_active_modes', max_triad_active_modes, ...
        'evaluation_model', evaluation_modes{i});
    t_triad = tic;
    triad = creamer_eta33_h3h4_triad_1d(eta_lin, x_vec, h, triad_cfg);
    runtime_s = toc(t_triad);

    triad_eta33_filtered = local_filter_by_mask(triad.eta33_h3h4, eta33_filter);
    triad_eta33_spec = abs(fft(triad.eta33_h3h4) / Nx);
    peak_val = max(abs(triad_eta33_filtered));
    amp_3kp = triad_eta33_spec(idx_3kp);
    pstats = triad.diagnostics.pruning_stats;

    rows(i).mode = string(evaluation_modes{i});
    rows(i).support_kind = string(triad.diagnostics.support_modes_kind);
    rows(i).runtime_s = runtime_s;
    rows(i).filtered_peak = peak_val;
    rows(i).peak_rel_to_mf12 = peak_val / max(mf12_peak, eps);
    rows(i).amp_3kp = amp_3kp;
    rows(i).amp_3kp_rel_to_mf12 = amp_3kp / max(mf12_3kp, eps);
    rows(i).candidate_terms = triad.diagnostics.candidate_terms;
    rows(i).used_terms = triad.diagnostics.used_terms;
    rows(i).h4_unique_rows = pstats.targeted_h4_unique_rows;
    rows(i).h4_support_pruned = pstats.targeted_h4_support_pruned;
    rows(i).bracket_products = pstats.targeted_bracket_products_formed;
    rows(i).bracket_support_pruned = pstats.targeted_bracket_support_pruned;
end

summary_table = struct2table(rows);

fprintf('Focused eta33 pruning benchmark\n');
fprintf('  Alpha = %.12g, Akp = %.12g, kph = %.12g, Nx = %d, max_triad_active_modes = %d\n', ...
    Alpha, Akp, kph, Nx, max_triad_active_modes);
fprintf('MF12 reference\n');
fprintf('  runtime_s      = %.3f\n', mf12_runtime_s);
fprintf('  filtered_peak  = %.12g\n', mf12_peak);
fprintf('  amp_3kp        = %.12g\n', mf12_3kp);
disp(summary_table);

benchmark_result = struct();
benchmark_result.Alpha = Alpha;
benchmark_result.Akp = Akp;
benchmark_result.h = h;
benchmark_result.kph = kph;
benchmark_result.Nx = Nx;
benchmark_result.max_triad_active_modes = max_triad_active_modes;
benchmark_result.k_over_kp = k_over_kp;
benchmark_result.x_vec = x_vec;
benchmark_result.eta_lin = eta_lin;
benchmark_result.mf12_runtime_s = mf12_runtime_s;
benchmark_result.mf12_filtered_peak = mf12_peak;
benchmark_result.mf12_3kp = mf12_3kp;
benchmark_result.summary_table = summary_table;

function eta_filtered = local_filter_by_mask(eta, mask)
eta_hat = fft(eta);
eta_hat(~mask) = 0;
eta_filtered = real(ifft(eta_hat));
end
