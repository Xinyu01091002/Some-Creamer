% Compare unidirectional four-phase-separated Creamer harmonics against live MF12.

clearvars -except Akp Alpha Nx x_vec h min_active_modes max_active_modes n_lambda_steps lambda_flow_model n_picard_iters energy_fraction t phishift;
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
if ~exist('min_active_modes', 'var') || isempty(min_active_modes)
    min_active_modes = 0;
end
if ~exist('n_lambda_steps', 'var') || isempty(n_lambda_steps)
    n_lambda_steps = 12;
end
if ~exist('lambda_flow_model', 'var') || isempty(lambda_flow_model)
    lambda_flow_model = 'canonical_pair';
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
    h = 1500;
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
cfg.propagation_direction_deg = 0;
cfg.preserve_mean = true;
cfg.verbose = true;

mf12_cfg = struct();
mf12_cfg.g = g;
mf12_cfg.h = h;
mf12_cfg.t = t;
mf12_cfg.energy_fraction = energy_fraction;
mf12_cfg.max_active_modes = max_active_modes;

t_mf12 = tic;
mf12 = mf12_from_linear_focus_1d(eta_lin, x_vec, mf12_cfg);
mf12_runtime_s = toc(t_mf12);

t_creamer = tic;
creamer_sep = creamer_four_phase_separation_1d(eta_lin, x_vec, cfg);
creamer_runtime_s = toc(t_creamer);

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
spec.creamer_eta33 = abs(fft(creamer_sep.eta33) / Nx);
spec.delta_mf12 = abs(fft(mf12.eta_total_order3 - eta_lin) / Nx);
spec.delta_creamer = abs(fft(creamer_sep.eta_nl_phase(1, :) - eta_lin) / Nx);

result = struct();
result.Akp = Akp;
result.Alpha = Alpha;
result.kp = kp;
result.lambda_p = lambda_p;
result.h = h;
result.kp_h = kp * h;
result.g = g;
result.x_vec = x_vec;
result.xlim_wave = [-4, 4];
result.k_over_kp = k_over_kp(pos_mask);
result.dx = dx;
result.points_per_wavelength = points_per_wavelength;
result.domain_wavelengths = domain_wavelengths;
result.max_active_modes = max_active_modes;
result.min_active_modes = min_active_modes;
result.n_lambda_steps = n_lambda_steps;
result.lambda_flow_model = lambda_flow_model;
result.mf12 = struct( ...
    'eta20', mf12.eta20, ...
    'eta22', mf12.eta22, ...
    'eta33', mf12.eta33, ...
    'delta_total', mf12.eta_total_order3 - eta_lin);
result.creamer = struct( ...
    'eta11', creamer_sep.eta1_total, ...
    'eta11_lin', eta_lin, ...
    'eta20', creamer_sep.eta20, ...
    'eta22', creamer_sep.eta22, ...
    'eta33', creamer_sep.eta33, ...
    'delta_total', creamer_sep.eta_nl_phase(1, :) - eta_lin);
result.spec = struct( ...
    'mf12_eta20', spec.mf12_eta20(pos_mask), ...
    'creamer_eta20', spec.creamer_eta20(pos_mask), ...
    'mf12_eta22', spec.mf12_eta22(pos_mask), ...
    'creamer_eta22', spec.creamer_eta22(pos_mask), ...
    'mf12_eta33', spec.mf12_eta33(pos_mask), ...
    'creamer_eta33', spec.creamer_eta33(pos_mask), ...
    'delta_mf12', spec.delta_mf12(pos_mask), ...
    'delta_creamer', spec.delta_creamer(pos_mask));

out_dir = fullfile(repo_root, 'MATLAB', 'output', 'unidirectional');
out_png = plot_unidirectional_creamer_comparison(result, out_dir);

fprintf('\nSaved unidirectional comparison figure to:\n  %s\n', out_png);
fprintf('MF12 runtime                            = %.3f s\n', mf12_runtime_s);
fprintf('Creamer four-phase runtime              = %.3f s\n', creamer_runtime_s);
fprintf('Grid spacing dx                         = %.6g m\n', dx);
fprintf('Points per wavelength                   = %.6g\n', points_per_wavelength);
fprintf('Domain length / lambda_p                = %.6g\n', domain_wavelengths);
fprintf('Creamer min active modes request        = %d\n', min_active_modes);
fprintf('MF12 active modes                       = %d\n', mf12.active_mode_count);
fprintf('Creamer active modes                    = %d\n', creamer_sep.diagnostics{1}.active_mode_count);
fprintf('Center max |MF12 eta20|                 = %.6g\n', max(abs(mf12.eta20)));
fprintf('Center max |Creamer eta20|              = %.6g\n', max(abs(creamer_sep.eta20)));
fprintf('Center max |MF12 eta22|                 = %.6g\n', max(abs(mf12.eta22)));
fprintf('Center max |Creamer eta22|              = %.6g\n', max(abs(creamer_sep.eta22)));
fprintf('Center max |MF12 eta33|                 = %.6g\n', max(abs(mf12.eta33)));
fprintf('Center max |Creamer eta33|              = %.6g\n', max(abs(creamer_sep.eta33)));
