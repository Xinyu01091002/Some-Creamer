function out = creamer_four_phase_separation_1d(eta_lin, x_vec, cfg)
%CREAMER_FOUR_PHASE_SEPARATION_1D Separate 1D harmonic sectors using four
%global phase shifts: 0, pi/2, pi, 3pi/2.
%
% In addition to eta20/eta22/eta33, this helper returns:
%   eta1_total = odd-fundamental content recovered from four-phase data
%   eta31      = eta1_total - eta_lin
% The eta31 naming is kept for continuity with the project notes, but the
% quantity should be read as an odd-fundamental correction that may include
% higher odd orders such as eta51.

phase_list = (0:3) * (pi / 2);
n_phase = numel(phase_list);
nx = numel(eta_lin);

eta_phase_shifted = zeros(n_phase, nx);
eta_nl_phase = zeros(n_phase, nx);
phi_nl_phase = zeros(n_phase, nx);
eta2_phase = zeros(n_phase, nx);
diagnostics = cell(1, n_phase);

for idx = 1:n_phase
    eta_phase_shifted(idx, :) = apply_global_phase_shift_1d(eta_lin, x_vec, phase_list(idx));
    if isfield(cfg, 'creamer_backend') && strcmpi(cfg.creamer_backend, 'cpp')
        phase_cfg = cfg;
        if isfield(cfg, 'cpp_job_dir')
            phase_cfg.cpp_job_dir = fullfile(cfg.cpp_job_dir, sprintf('phase_%d', idx));
        end
        [eta_nl_i, diagnostics{idx}, phi_nl_i] = directional_creamer_transform_cpp( ...
            reshape(eta_phase_shifted(idx, :), 1, []), x_vec, 0, phase_cfg);
    else
        [eta_nl_i, diagnostics{idx}, phi_nl_i] = directional_creamer_transform( ...
            reshape(eta_phase_shifted(idx, :), 1, []), x_vec, 0, cfg);
    end
    eta_nl_phase(idx, :) = eta_nl_i;
    phi_nl_phase(idx, :) = phi_nl_i;
    eta2_phase(idx, :) = eta_nl_phase(idx, :) - eta_phase_shifted(idx, :);
end

coef = [
    0.25  0    -0.25  0     0    -0.25  0     0.25
    0.25 -0.25  0.25 -0.25  0     0     0     0
    0.25  0    -0.25  0     0     0.25  0    -0.25
    0.25  0.25  0.25  0.25  0     0     0     0
];

phase_hilbert = zeros(n_phase, nx);
phi_phase_hilbert = zeros(n_phase, nx);
for idx = 1:n_phase
    phase_hilbert(idx, :) = -imag(hilbert(eta_nl_phase(idx, :)));
    phi_phase_hilbert(idx, :) = -imag(hilbert(phi_nl_phase(idx, :)));
end

all_eta = [eta_nl_phase; phase_hilbert];
all_phi = [phi_nl_phase; phi_phase_hilbert];

eta1_total = zeros(1, nx);
eta22 = zeros(1, nx);
eta33 = zeros(1, nx);
eta20 = zeros(1, nx);
phi1_total = zeros(1, nx);
phi22 = zeros(1, nx);
phi33 = zeros(1, nx);
phi20 = zeros(1, nx);
for idx = 1:8
    eta1_total = eta1_total + coef(1, idx) * all_eta(idx, :);
    eta22 = eta22 + coef(2, idx) * all_eta(idx, :);
    eta33 = eta33 + coef(3, idx) * all_eta(idx, :);
    eta20 = eta20 + coef(4, idx) * all_eta(idx, :);
    phi1_total = phi1_total + coef(1, idx) * all_phi(idx, :);
    phi22 = phi22 + coef(2, idx) * all_phi(idx, :);
    phi33 = phi33 + coef(3, idx) * all_phi(idx, :);
    phi20 = phi20 + coef(4, idx) * all_phi(idx, :);
end

eta31 = eta1_total - eta_lin;

out = struct();
out.coef = coef;
out.phase_list = phase_list;
out.eta_phase_shifted = eta_phase_shifted;
out.eta_nl_phase = eta_nl_phase;
out.phi_nl_phase = phi_nl_phase;
out.eta2_phase = eta2_phase;
out.phase_hilbert = phase_hilbert;
out.phi_phase_hilbert = phi_phase_hilbert;
out.eta1_total = eta1_total;
out.eta31 = eta31;
out.eta20 = eta20;
out.eta22 = eta22;
out.eta33 = eta33;
out.phi1_total = phi1_total;
out.phi20 = phi20;
out.phi22 = phi22;
out.phi33 = phi33;
out.diagnostics = diagnostics;
end
