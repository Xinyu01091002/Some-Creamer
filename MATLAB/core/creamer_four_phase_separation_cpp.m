function out = creamer_four_phase_separation_cpp(eta_lin, x_vec, y_vec, cfg)
%CREAMER_FOUR_PHASE_SEPARATION_CPP Four-phase separation using C++ backend.

phase_list = (0:3) * (pi / 2);
n_phase = numel(phase_list);
[ny, nx] = size(eta_lin);

eta_phase_shifted = zeros(ny, nx, n_phase);
eta_nl_phase = zeros(ny, nx, n_phase);
phi_nl_phase = zeros(ny, nx, n_phase);
diagnostics = cell(1, n_phase);

base_job_dir = cfg.cpp_job_dir;
for idx = 1:n_phase
    phase_cfg = cfg;
    phase_cfg.cpp_job_dir = fullfile(base_job_dir, sprintf('phase_%d', idx));
    eta_phase_shifted(:, :, idx) = apply_global_phase_shift_2d(eta_lin, x_vec, y_vec, phase_list(idx));
    [eta_nl_phase(:, :, idx), diagnostics{idx}, phi_nl_phase(:, :, idx)] = directional_creamer_transform_cpp( ...
        eta_phase_shifted(:, :, idx), x_vec, y_vec, phase_cfg);
end

coef = [
    0.25  0    -0.25  0     0    -0.25  0     0.25
    0.25 -0.25  0.25 -0.25  0     0     0     0
    0.25  0    -0.25  0     0     0.25  0    -0.25
    0.25  0.25  0.25  0.25  0     0     0     0
];

phase_hilbert = zeros(ny, nx, n_phase);
phi_phase_hilbert = zeros(ny, nx, n_phase);
for idx = 1:n_phase
    phase_hilbert(:, :, idx) = -imag(hilbert(eta_nl_phase(:, :, idx).').');
    phi_phase_hilbert(:, :, idx) = -imag(hilbert(phi_nl_phase(:, :, idx).').');
end

all_eta = cat(3, eta_nl_phase, phase_hilbert);
all_phi = cat(3, phi_nl_phase, phi_phase_hilbert);

eta1_total = zeros(ny, nx);
eta22 = zeros(ny, nx);
eta33 = zeros(ny, nx);
eta20 = zeros(ny, nx);
phi1_total = zeros(ny, nx);
phi22 = zeros(ny, nx);
phi33 = zeros(ny, nx);
phi20 = zeros(ny, nx);
for idx = 1:8
    eta1_total = eta1_total + coef(1, idx) * all_eta(:, :, idx);
    eta22 = eta22 + coef(2, idx) * all_eta(:, :, idx);
    eta33 = eta33 + coef(3, idx) * all_eta(:, :, idx);
    eta20 = eta20 + coef(4, idx) * all_eta(:, :, idx);
    phi1_total = phi1_total + coef(1, idx) * all_phi(:, :, idx);
    phi22 = phi22 + coef(2, idx) * all_phi(:, :, idx);
    phi33 = phi33 + coef(3, idx) * all_phi(:, :, idx);
    phi20 = phi20 + coef(4, idx) * all_phi(:, :, idx);
end

out = struct();
out.coef = coef;
out.phase_list = phase_list;
out.eta_phase_shifted = eta_phase_shifted;
out.eta_nl_phase = eta_nl_phase;
out.phi_nl_phase = phi_nl_phase;
out.phase_hilbert = phase_hilbert;
out.phi_phase_hilbert = phi_phase_hilbert;
out.eta1_total = eta1_total;
out.eta20 = eta20;
out.eta22 = eta22;
out.eta33 = eta33;
out.phi1_total = phi1_total;
out.phi20 = phi20;
out.phi22 = phi22;
out.phi33 = phi33;
out.diagnostics = diagnostics;
end
