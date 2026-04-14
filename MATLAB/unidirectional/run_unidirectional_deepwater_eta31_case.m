% Spatial-only runner for the unidirectional deep-water eta31 workflow.
%
% This runner keeps a finite reference depth k_p h(ref) = 10 for reporting,
% but explicitly locks the reconstruction kernel to the 1989 deep-water H3
% model through cfg.kernel_model = 'deep_water_1989'.

clearvars -except Akp Alpha Nx x_vec kp_h_ref min_active_modes max_active_modes n_lambda_steps n_picard_iters energy_fraction t phishift;
clc;

this_dir = fileparts(mfilename('fullpath'));
matlab_dir = fileparts(this_dir);
repo_root = fileparts(matlab_dir);
addpath(this_dir);
addpath(fullfile(matlab_dir, 'core'));

if ~exist('Akp', 'var') || isempty(Akp)
    Akp = 0.12;
end
if ~exist('Alpha', 'var') || isempty(Alpha)
    Alpha = 8;
end
if ~exist('Nx', 'var') || isempty(Nx)
    Nx = 4096;
end
if ~exist('kp_h_ref', 'var') || isempty(kp_h_ref)
    kp_h_ref = 10;
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
h_ref = kp_h_ref / kp;
g = 9.81;

if ~exist('x_vec', 'var') || isempty(x_vec)
    x_vec = linspace(-64 * lambda_p, 64 * lambda_p, Nx);
else
    x_vec = x_vec(:).';
    Nx = numel(x_vec);
end

eta11_input = real(Linear_focus_envelope(Akp, Alpha, x_vec, t, phishift));
eta11_input = reshape(eta11_input, 1, []);

cfg = struct();
cfg.g = g;
cfg.depth_h = h_ref;
cfg.kernel_model = 'deep_water_1989';
cfg.energy_fraction = energy_fraction;
cfg.min_active_modes = min_active_modes;
cfg.max_active_modes = max_active_modes;
cfg.n_lambda_steps = n_lambda_steps;
cfg.n_picard_iters = n_picard_iters;
cfg.lambda_flow_model = 'canonical_pair';
cfg.propagation_direction_deg = 0;
cfg.preserve_mean = true;
cfg.verbose = true;

t_creamer = tic;
creamer_sep = creamer_four_phase_separation_1d(eta11_input, x_vec, cfg);
creamer_runtime_s = toc(t_creamer);
t_vwa = tic;
vwa_eta33 = vwa_broadband_eta33_1d(eta11_input, x_vec, h_ref, g);
vwa_runtime_s = toc(t_vwa);

closure_err = max(abs(creamer_sep.eta1_total - (eta11_input + creamer_sep.eta31)));
dx = mean(diff(x_vec));
Lx = dx * Nx;
mx = [0:(Nx/2), (-Nx/2 + 1):-1];
k = (2 * pi / Lx) * mx;
pos_mask = (mx >= 0);

spec = struct();
spec.eta11_input = abs(fft(eta11_input) / Nx);
spec.eta1_total = abs(fft(creamer_sep.eta1_total) / Nx);
spec.eta31 = abs(fft(creamer_sep.eta31) / Nx);
spec.eta33 = abs(fft(creamer_sep.eta33) / Nx);
spec.vwa_eta33 = abs(vwa_eta33.eta33_hat);
eta11_phase = angle(hilbert(eta11_input));
vwa_eta31 = -real((abs(hilbert(creamer_sep.eta33)) / 3) .* exp(1i * eta11_phase));
spec.vwa_eta31 = abs(fft(vwa_eta31) / Nx);

result = struct();
result.Akp = Akp;
result.Alpha = Alpha;
result.kp = kp;
result.lambda_p = lambda_p;
result.k_over_kp = k(pos_mask) / kp;
result.kp_h_ref = kp_h_ref;
result.h_ref = h_ref;
result.g = g;
result.x_vec = x_vec;
result.xlim_wave = [-4, 4];
result.min_active_modes = min_active_modes;
result.max_active_modes = max_active_modes;
result.n_lambda_steps = n_lambda_steps;
result.kernel_model = 'deep_water_1989';
result.lambda_flow_model = 'canonical_pair';
result.eta11_input = eta11_input;
result.eta1_total = creamer_sep.eta1_total;
result.eta31 = creamer_sep.eta31;
result.eta33 = creamer_sep.eta33;
result.vwa_eta33 = vwa_eta33.eta33;
result.vwa_eta31 = vwa_eta31;
result.spec = struct( ...
    'eta11_input', spec.eta11_input(pos_mask), ...
    'eta1_total', spec.eta1_total(pos_mask), ...
    'eta31', spec.eta31(pos_mask), ...
    'eta33', spec.eta33(pos_mask), ...
    'vwa_eta33', spec.vwa_eta33(pos_mask), ...
    'vwa_eta31', spec.vwa_eta31(pos_mask));
result.creamer_runtime_s = creamer_runtime_s;
result.vwa_runtime_s = vwa_runtime_s;
result.eta31_closure_err = closure_err;
result.phase_diagnostics = creamer_sep.diagnostics;
result.diagnostics = creamer_sep.diagnostics{1};

out_dir = fullfile(repo_root, 'MATLAB', 'output', 'unidirectional', 'eta31');
out_png = plot_unidirectional_eta31_profile(result, out_dir);

fprintf('\nSaved deep-water eta31 figure to:\n  %s\n', out_png);
fprintf('Creamer four-phase runtime              = %.3f s\n', creamer_runtime_s);
fprintf('Broadband VWA eta33 runtime             = %.3f s\n', vwa_runtime_s);
fprintf('Kernel model                            = %s\n', result.diagnostics.kernel_model);
fprintf('Depth model                             = %s\n', result.diagnostics.depth_model);
fprintf('Reference k_p h                         = %.6g\n', result.kp_h_ref);
fprintf('Reference h                             = %.6g m\n', result.h_ref);
fprintf('Active modes                            = %d\n', result.diagnostics.active_mode_count);
fprintf('Max |eta11 input|                       = %.6g\n', max(abs(result.eta11_input)));
fprintf('Max |eta1 total|                        = %.6g\n', max(abs(result.eta1_total)));
fprintf('Max |eta31|                             = %.6g\n', max(abs(result.eta31)));
fprintf('Max |eta33|                             = %.6g\n', max(abs(result.eta33)));
fprintf('Max |Broadband VWA eta33|               = %.6g\n', max(abs(result.vwa_eta33)));
fprintf('Max |Phase-modulated VWA eta31|         = %.6g\n', max(abs(result.vwa_eta31)));
fprintf('Closure max |eta1-(eta11+eta31)|        = %.6g\n', result.eta31_closure_err);
