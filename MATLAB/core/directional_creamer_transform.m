function [eta_nl, diagnostics, phi_nl] = directional_creamer_transform(eta_lin, x_vec, y_vec, cfg)
% DIRECTIONAL_CREAMER_TRANSFORM
% Directional Creamer prototype for deep-water and optional finite-depth waves.
%
% This implementation uses the general 2D kernels B and D from Creamer et al.
% (1989), equations (3.4), (3.5), for deep water. By default the kernel
% choice follows the historical depth_h convention:
%   - depth_h = inf   -> deep-water 1989 kernels
%   - finite depth_h  -> finite-depth Wright & Creamer (1994) kernels
% cfg.kernel_model can override that inference explicitly and force either
% 'deep_water_1989' or 'finite_depth_1994'. The inverse map is implied by
% solving the Lie-transform equations from lambda = 1 to lambda = 0,
% truncated at the first non-trivial order:
%
%   zeta ~= zeta_bar - [ 3 D phi_bar phi_bar + B zeta_bar zeta_bar ].
%
% The goal of this routine is to produce a usable nonlinear surface
% eta_nl(x,y) from one linear directional parent surface eta_lin(x,y).
%
% By default (n_lambda_steps = 0) the routine returns the original
% first-nontrivial-order truncation. When n_lambda_steps > 0, it instead
% integrates a lambda-flow with RK4 from lambda=1 to lambda=0.
%
% Two lambda-flow closures are supported:
%   - 'canonical_pair'  : evolve (Z,Phi) together using the paper-level
%                         coupled equations (3.12a,b)
%   - 'backward_picard_316'
%                       : solve the backward integral equations (3.16a,b)
%                         by Picard iteration on a lambda grid
%   - 'legacy_zeta_only': evolve Z only and rebuild Phi from Z via the
%                         linear deep-water relation at each stage
%
% The canonical-pair flow is the default because it is structurally closer
% to the Lie-transform formulation in the paper.
%
% Important modelling choice in v1:
% The input only specifies eta_bar(x,y), not the corresponding linear
% surface potential phi_bar(x,y). We therefore reconstruct phi_bar in
% Fourier space using a positive-frequency convention aligned with
% cfg.propagation_direction_deg. In finite depth this uses
% theta(k)=|k| tanh(|k| h) rather than the deep-water theta(k)=|k|.
% This is the main additional closure assumption in the current prototype.

if nargin < 4 || isempty(cfg)
    cfg = struct();
end

cfg = local_apply_defaults(cfg);
[kernel_model, kernel_depth_h] = local_resolve_kernel_model(cfg);

if ~isvector(x_vec) || ~isvector(y_vec)
    error('directional_creamer_transform:VectorRequired', ...
        'x_vec and y_vec must both be vectors.');
end

[ny, nx] = size(eta_lin);

if numel(x_vec) ~= nx
    error('directional_creamer_transform:SizeMismatch', ...
        'Length of x_vec (%d) must match size(eta_lin,2) (%d).', numel(x_vec), nx);
end

if numel(y_vec) ~= ny
    error('directional_creamer_transform:SizeMismatch', ...
        'Length of y_vec (%d) must match size(eta_lin,1) (%d).', numel(y_vec), ny);
end

x_vec = x_vec(:).';
y_vec = y_vec(:).';

dx = mean(diff(x_vec));
if ny > 1
    dy = mean(diff(y_vec));
else
    dy = 1;
end
Lx = dx * nx;
Ly = dy * ny;
n_total = nx * ny;

[kx, mx] = local_fft_wavenumbers(nx, dx);
[ky, my] = local_fft_wavenumbers(ny, dy);
[KX, KY] = meshgrid(kx, ky);
Kmag = hypot(KX, KY);
Theta = local_theta_symbol(Kmag, kernel_model, kernel_depth_h);

zeta_hat = fft2(eta_lin) / n_total;
phi_bar_hat = local_linear_phi_from_eta(zeta_hat, Theta, KX, KY, cfg.g, cfg.propagation_direction_deg);

if ~cfg.preserve_mean
    zeta_hat(1, 1) = 0;
    phi_bar_hat(1, 1) = 0;
end

active_mask = local_build_active_mask(zeta_hat, cfg.energy_fraction, cfg.min_active_modes, cfg.max_active_modes, cfg.preserve_mean);
n_active = nnz(active_mask);
active_plan = local_build_active_plan(active_mask, Kmag, Theta, mx, ny, nx, cfg.g, cfg.plan_chunk_size, kernel_model, kernel_depth_h);

if cfg.n_lambda_steps <= 0
    delta_zeta_hat = local_inverse_zeta_correction(zeta_hat, phi_bar_hat, active_plan, ny, nx);
    eta_nl_hat = zeta_hat + delta_zeta_hat;
    phi_nl_hat = phi_bar_hat + local_inverse_phi_correction(zeta_hat, phi_bar_hat, active_plan, ny, nx);
else
    zeta_state = zeta_hat;
    phi_state = phi_bar_hat;
    if strcmpi(cfg.lambda_flow_model, 'backward_picard_316')
        [zeta_state, phi_state] = local_backward_picard_316( ...
            zeta_hat, phi_bar_hat, active_mask, active_plan, Theta, kx, mx, ny, nx, cfg);
    elseif strcmpi(cfg.lambda_stepper, 'adaptive_rk4')
        [zeta_state, phi_state] = local_adaptive_lambda_flow( ...
            zeta_state, phi_state, active_plan, Theta, KX, KY, kx, mx, ny, nx, cfg);
    else
        h_lambda = -1 / cfg.n_lambda_steps;
        for step = 1:cfg.n_lambda_steps
            [zeta_state, phi_state] = local_take_fixed_lambda_step( ...
                zeta_state, phi_state, h_lambda, active_plan, Theta, KX, KY, kx, mx, ny, nx, cfg);

            if ~cfg.preserve_mean
                zeta_state(1,1) = 0;
                phi_state(1,1) = 0;
            end
        end
    end
    eta_nl_hat = zeta_state;
    phi_nl_hat = phi_state;
end

if cfg.preserve_mean
    eta_nl_hat(1, 1) = eta_nl_hat(1, 1);
else
    eta_nl_hat(1, 1) = 0;
end

eta_nl = real(ifft2(eta_nl_hat * n_total));
phi_nl = real(ifft2(phi_nl_hat * n_total));

diagnostics = struct();
diagnostics.dx = dx;
diagnostics.dy = dy;
diagnostics.Lx = Lx;
diagnostics.Ly = Ly;
diagnostics.nx = nx;
diagnostics.ny = ny;
diagnostics.active_mode_count = n_active;
diagnostics.energy_fraction = cfg.energy_fraction;
diagnostics.n_lambda_steps = cfg.n_lambda_steps;
diagnostics.lambda_flow_model = cfg.lambda_flow_model;
diagnostics.lambda_stepper = cfg.lambda_stepper;
diagnostics.max_eta_lin = max(abs(eta_lin(:)));
diagnostics.max_eta_nl = max(abs(eta_nl(:)));
diagnostics.max_imag_eta_nl = max(abs(imag(ifft2(eta_nl_hat * n_total))), [], 'all');
diagnostics.max_imag_phi_nl = max(abs(imag(ifft2(phi_nl_hat * n_total))), [], 'all');
diagnostics.max_abs_phi_hat = max(abs(phi_nl_hat(:)));
diagnostics.input_mean = mean(eta_lin(:));
diagnostics.output_mean = mean(eta_nl(:));
diagnostics.propagation_direction_deg = cfg.propagation_direction_deg;
diagnostics.depth_h = cfg.depth_h;
diagnostics.kernel_depth_h = kernel_depth_h;
diagnostics.kernel_model = kernel_model;
diagnostics.depth_model = active_plan.depth_model;

if cfg.verbose
    fprintf('Directional Creamer prototype summary:\n');
    fprintf('  Grid                : ny=%d, nx=%d\n', ny, nx);
    fprintf('  Spacing             : dx=%.6g, dy=%.6g\n', dx, dy);
    fprintf('  Active modes        : %d\n', n_active);
    fprintf('  Lambda RK steps     : %d\n', cfg.n_lambda_steps);
    fprintf('  Lambda flow model   : %s\n', cfg.lambda_flow_model);
    fprintf('  Lambda stepper      : %s\n', cfg.lambda_stepper);
    fprintf('  Kernel model        : %s\n', diagnostics.kernel_model);
    fprintf('  Depth model         : %s\n', diagnostics.depth_model);
    if isfinite(cfg.depth_h)
        fprintf('  Depth input h       : %.6g\n', cfg.depth_h);
    end
    fprintf('  Energy kept         : %.2f%%\n', 100 * cfg.energy_fraction);
    fprintf('  Max |eta_lin|       : %.6g\n', diagnostics.max_eta_lin);
    fprintf('  Max |eta_nl|        : %.6g\n', diagnostics.max_eta_nl);
    fprintf('  Mean eta_lin        : %.6g\n', diagnostics.input_mean);
    fprintf('  Mean eta_nl         : %.6g\n', diagnostics.output_mean);
    fprintf('  Max imag residual   : %.6g\n', diagnostics.max_imag_eta_nl);
    fprintf('  Max imag phi resid  : %.6g\n', diagnostics.max_imag_phi_nl);
end

end

function cfg = local_apply_defaults(cfg)
defaults = struct( ...
    'g', 9.81, ...
    'energy_fraction', 0.99, ...
    'min_active_modes', 0, ...
    'max_active_modes', inf, ...
    'n_lambda_steps', 0, ...
    'n_picard_iters', 4, ...
    'lambda_flow_model', 'canonical_pair', ...
    'lambda_stepper', 'fixed_rk4', ...
    'lambda_rtol', 1e-6, ...
    'lambda_atol', 1e-9, ...
    'lambda_initial_step', [], ...
    'lambda_min_step', 1e-4, ...
    'lambda_max_step', 0.25, ...
    'depth_h', inf, ...
    'kernel_model', '', ...
    'plan_chunk_size', 128, ...
    'propagation_direction_deg', 0, ...
    'preserve_mean', true, ...
    'verbose', true);

names = fieldnames(defaults);
for n = 1:numel(names)
    name = names{n};
    if ~isfield(cfg, name) || isempty(cfg.(name))
        cfg.(name) = defaults.(name);
    end
end
end

function [kernel_model, kernel_depth_h] = local_resolve_kernel_model(cfg)
if ~isfield(cfg, 'kernel_model') || isempty(cfg.kernel_model)
    if isfinite(cfg.depth_h)
        kernel_model = 'finite_depth_1994';
    else
        kernel_model = 'deep_water_1989';
    end
else
    kernel_model = char(cfg.kernel_model);
end

if strcmpi(kernel_model, 'deep_water_1989')
    kernel_model = 'deep_water_1989';
    kernel_depth_h = inf;
elseif strcmpi(kernel_model, 'finite_depth_1994')
    if ~isfinite(cfg.depth_h)
        error('directional_creamer_transform:FiniteDepthRequired', ...
            'cfg.depth_h must be finite when cfg.kernel_model is finite_depth_1994.');
    end
    kernel_model = 'finite_depth_1994';
    kernel_depth_h = cfg.depth_h;
else
    error('directional_creamer_transform:UnknownKernelModel', ...
        'Unknown cfg.kernel_model: %s', kernel_model);
end
end

function theta = local_theta_symbol(Kmag, kernel_model, depth_h)
if strcmpi(kernel_model, 'finite_depth_1994')
    theta = Kmag .* tanh(Kmag * depth_h);
else
    theta = Kmag;
end
theta(Kmag == 0) = 0;
end

function active_mask = local_build_active_mask(zeta_hat, energy_fraction, min_active_modes, max_active_modes, preserve_mean)
energy_density = abs(zeta_hat).^2;
flat_energy = energy_density(:);
[sorted_energy, order] = sort(flat_energy, 'descend');
cumulative_fraction = cumsum(sorted_energy) / sum(sorted_energy);

keep_count = find(cumulative_fraction >= energy_fraction, 1, 'first');
if isempty(keep_count)
    keep_count = numel(order);
end
if isfinite(min_active_modes)
    keep_count = max(keep_count, min(min_active_modes, numel(order)));
end
if isfinite(max_active_modes)
    keep_count = min(keep_count, max_active_modes);
end

active_mask = false(size(zeta_hat));
active_mask(order(1:keep_count)) = true;
if preserve_mean && abs(zeta_hat(1,1)) > 0
    active_mask(1,1) = true;
else
    active_mask(1,1) = false;
end
end

function plan = local_build_active_plan(active_mask, Kmag, Theta, mx, ny, nx, g, chunk_size, kernel_model, depth_h)
[row_idx, col_idx] = find(active_mask);
active_idx = sub2ind([ny, nx], row_idx, col_idx);
n_active = numel(active_idx);

my = local_fft_mode_numbers(ny);
finite_depth = strcmpi(kernel_model, 'finite_depth_1994');
plan = struct();
plan.active_idx = active_idx;
plan.mx = mx(col_idx).';
plan.my = my(row_idx).';
plan.kmag = Kmag(active_idx);
plan.theta = Theta(active_idx);
plan.n_active = n_active;
plan.chunk_size = max(1, round(chunk_size));
plan.depth_h = depth_h;
plan.kernel_model = kernel_model;
if finite_depth
    plan.depth_model = 'finite_depth_1994';
else
    plan.depth_model = 'deep_water_1989';
end
plan.dest_idx = zeros(n_active, n_active, 'uint32');
plan.Bz = zeros(n_active, n_active);
plan.Bphi = zeros(n_active, n_active);
plan.D = zeros(n_active, n_active);

for a = 1:n_active
    mx_a = plan.mx(a);
    my_a = plan.my(a);
    kmag_a = plan.kmag(a);
    theta_a = plan.theta(a);
    ksq_a = kmag_a^2;

    for b = 1:n_active
        kmag_b = plan.kmag(b);
        theta_b = plan.theta(b);
        ksq_b = kmag_b^2;
        col_k = local_mode_to_index(mx_a + plan.mx(b), nx);
        row_k = local_mode_to_index(my_a + plan.my(b), ny);
        kmag_k = Kmag(row_k, col_k);
        theta_k = Theta(row_k, col_k);
        ksq_k = kmag_k^2;

        plan.dest_idx(a, b) = uint32(sub2ind([ny, nx], row_k, col_k));
        if finite_depth
            plan.D(a, b) = local_kernel_D_finite(theta_k, theta_a, theta_b, ksq_k, ksq_a, ksq_b, g);
            plan.Bz(a, b) = local_kernel_B_finite(theta_a, theta_b, theta_k, ksq_a, ksq_b, ksq_k);
            plan.Bphi(a, b) = local_kernel_B_finite(theta_k, theta_a, theta_b, ksq_k, ksq_a, ksq_b);
        else
            plan.D(a, b) = local_kernel_D_deep(kmag_k, kmag_a, kmag_b, g);
            plan.Bz(a, b) = local_kernel_B_deep(kmag_a, kmag_b, kmag_k);
            plan.Bphi(a, b) = local_kernel_B_deep(kmag_k, kmag_a, kmag_b);
        end
    end
end
end

function delta_zeta_hat = local_inverse_zeta_correction(zeta_hat, phi_hat, plan, ny, nx)
n_active = plan.n_active;
zeta_active = zeta_hat(plan.active_idx);
phi_active = phi_hat(plan.active_idx);
delta_zeta_hat = zeros(ny, nx);
for a = 1:n_active
    zeta_a = zeta_active(a);
    phi_a = phi_active(a);

    for b = 1:n_active
        dest = double(plan.dest_idx(a, b));
        delta_zeta_hat(dest) = delta_zeta_hat(dest) ...
            - (3 * plan.D(a, b) * phi_a * phi_active(b) + plan.Bz(a, b) * zeta_a * zeta_active(b));
    end
end
end

function delta_phi_hat = local_inverse_phi_correction(zeta_hat, phi_hat, plan, ny, nx)
n_active = plan.n_active;
zeta_active = zeta_hat(plan.active_idx);
phi_active = phi_hat(plan.active_idx);
delta_phi_hat = zeros(ny, nx);
for a = 1:n_active
    zeta_a = zeta_active(a);

    for b = 1:n_active
        dest = double(plan.dest_idx(a, b));
        delta_phi_hat(dest) = delta_phi_hat(dest) + 2 * plan.Bphi(a, b) * zeta_a * phi_active(b);
    end
end
end

function [zeta_next, phi_next] = local_take_fixed_lambda_step(zeta_state, phi_state, h_lambda, plan, Theta, KX, KY, kx, mx, ny, nx, cfg)
if strcmpi(cfg.lambda_flow_model, 'canonical_pair')
    [k1z, k1p] = local_coupled_lambda_rhs(zeta_state, phi_state, plan, ny, nx);
    [k2z, k2p] = local_coupled_lambda_rhs( ...
        zeta_state + 0.5 * h_lambda * k1z, ...
        phi_state + 0.5 * h_lambda * k1p, ...
        plan, ny, nx);
    [k3z, k3p] = local_coupled_lambda_rhs( ...
        zeta_state + 0.5 * h_lambda * k2z, ...
        phi_state + 0.5 * h_lambda * k2p, ...
        plan, ny, nx);
    [k4z, k4p] = local_coupled_lambda_rhs( ...
        zeta_state + h_lambda * k3z, ...
        phi_state + h_lambda * k3p, ...
        plan, ny, nx);

    zeta_next = zeta_state + (h_lambda / 6) * (k1z + 2*k2z + 2*k3z + k4z);
    phi_next = phi_state + (h_lambda / 6) * (k1p + 2*k2p + 2*k3p + k4p);
else
    k1z = local_legacy_lambda_rhs(zeta_state, plan, Theta, KX, KY, cfg);
    k2z = local_legacy_lambda_rhs(zeta_state + 0.5 * h_lambda * k1z, plan, Theta, KX, KY, cfg);
    k3z = local_legacy_lambda_rhs(zeta_state + 0.5 * h_lambda * k2z, plan, Theta, KX, KY, cfg);
    k4z = local_legacy_lambda_rhs(zeta_state + h_lambda * k3z, plan, Theta, KX, KY, cfg);
    zeta_next = zeta_state + (h_lambda / 6) * (k1z + 2*k2z + 2*k3z + k4z);
    phi_next = local_linear_phi_from_eta(zeta_next, Theta, KX, KY, cfg.g, cfg.propagation_direction_deg);
end
end

function [zeta_final, phi_final] = local_adaptive_lambda_flow(zeta_state, phi_state, plan, Theta, KX, KY, kx, mx, ny, nx, cfg)
lambda_curr = 1;
if isempty(cfg.lambda_initial_step)
    h = min(cfg.lambda_max_step, max(cfg.lambda_min_step, 1 / cfg.n_lambda_steps));
else
    h = min(cfg.lambda_max_step, max(cfg.lambda_min_step, abs(cfg.lambda_initial_step)));
end

while lambda_curr > 0
    h = min(h, lambda_curr);
    h_step = -h;

    [zeta_full, phi_full] = local_take_fixed_lambda_step( ...
        zeta_state, phi_state, h_step, plan, Theta, KX, KY, kx, mx, ny, nx, cfg);

    [zeta_half, phi_half] = local_take_fixed_lambda_step( ...
        zeta_state, phi_state, 0.5 * h_step, plan, Theta, KX, KY, kx, mx, ny, nx, cfg);
    [zeta_half2, phi_half2] = local_take_fixed_lambda_step( ...
        zeta_half, phi_half, 0.5 * h_step, plan, Theta, KX, KY, kx, mx, ny, nx, cfg);

    err = local_relative_pair_error(zeta_full, phi_full, zeta_half2, phi_half2, cfg.lambda_atol, cfg.lambda_rtol);
    if err <= 1
        zeta_state = zeta_half2;
        phi_state = phi_half2;
        lambda_curr = lambda_curr - h;
        if ~cfg.preserve_mean
            zeta_state(1,1) = 0;
            phi_state(1,1) = 0;
        end
        if err > 0
            h = min(cfg.lambda_max_step, max(cfg.lambda_min_step, 0.9 * h * err^(-0.2)));
        else
            h = min(cfg.lambda_max_step, 2.0 * h);
        end
    else
        h = max(cfg.lambda_min_step, 0.9 * h * err^(-0.2));
        if h <= cfg.lambda_min_step + eps
            zeta_state = zeta_half2;
            phi_state = phi_half2;
            lambda_curr = lambda_curr - h;
        end
    end
end

zeta_final = zeta_state;
phi_final = phi_state;
end

function err = local_relative_pair_error(zeta_a, phi_a, zeta_b, phi_b, atol, rtol)
scale_z = atol + rtol * max(abs(zeta_a), abs(zeta_b));
scale_p = atol + rtol * max(abs(phi_a), abs(phi_b));
err_z = max(abs(zeta_a - zeta_b) ./ max(scale_z, eps), [], 'all');
err_p = max(abs(phi_a - phi_b) ./ max(scale_p, eps), [], 'all');
err = max(err_z, err_p);
end

function delta_zeta_hat = local_legacy_lambda_rhs(zeta_hat, plan, Theta, KX, KY, cfg)
phi_hat = local_linear_phi_from_eta(zeta_hat, Theta, KX, KY, cfg.g, cfg.propagation_direction_deg);
ny = size(zeta_hat, 1);
nx = size(zeta_hat, 2);
n_active = plan.n_active;
zeta_active = zeta_hat(plan.active_idx);
phi_active = phi_hat(plan.active_idx);
delta_zeta_hat = zeros(ny, nx);
for a = 1:n_active
    zeta_a = zeta_active(a);
    phi_a = phi_active(a);

    for b = 1:n_active
        dest = double(plan.dest_idx(a, b));
        delta_zeta_hat(dest) = delta_zeta_hat(dest) ...
            - (3 * plan.D(a, b) * phi_a * phi_active(b) + plan.Bz(a, b) * zeta_a * zeta_active(b));
    end
end
end

function [dzeta_hat, dphi_hat] = local_coupled_lambda_rhs(zeta_hat, phi_hat, plan, ny, nx)
n_active = plan.n_active;
zeta_active = zeta_hat(plan.active_idx);
phi_active = phi_hat(plan.active_idx);

dzeta_hat = zeros(ny, nx);
dphi_hat = zeros(ny, nx);

for a = 1:n_active
    zeta_a = zeta_active(a);
    phi_a = phi_active(a);

    for b = 1:n_active
        dest = double(plan.dest_idx(a, b));
        dzeta_hat(dest) = dzeta_hat(dest) ...
            + 3 * plan.D(a, b) * phi_a * phi_active(b) ...
            + plan.Bz(a, b) * zeta_a * zeta_active(b);
        dphi_hat(dest) = dphi_hat(dest) - 2 * plan.Bphi(a, b) * zeta_a * phi_active(b);
    end
end
end

function [zeta_final, phi_final] = local_backward_picard_316(zeta_bar_hat, phi_bar_hat, active_mask, plan, Theta, kx, mx, ny, nx, cfg)
active_idx = plan.active_idx;
n_active = numel(active_idx);
n_nodes = cfg.n_lambda_steps + 1;
h = 1 / cfg.n_lambda_steps;

zeta_bar_active = zeta_bar_hat(active_idx);
phi_bar_active = phi_bar_hat(active_idx);

zeta_path = repmat(zeta_bar_active, 1, n_nodes);
phi_path = repmat(phi_bar_active, 1, n_nodes);

for iter = 1:cfg.n_picard_iters
    rhs_z_active = zeros(n_active, n_nodes - 1);
    rhs_p_active = zeros(n_active, n_nodes - 1);

    for node = 1:(n_nodes - 1)
        zeta_node = zeros(ny, nx);
        phi_node = zeros(ny, nx);
        zeta_node(active_idx) = zeta_path(:, node);
        phi_node(active_idx) = phi_path(:, node);

        [rhs_z_full, rhs_p_full] = local_coupled_lambda_rhs( ...
            zeta_node, phi_node, plan, ny, nx);

        rhs_z_active(:, node) = rhs_z_full(active_idx);
        rhs_p_active(:, node) = rhs_p_full(active_idx);
    end

    new_zeta_path = zeros(n_active, n_nodes);
    new_phi_path = zeros(n_active, n_nodes);
    new_zeta_path(:, 1) = zeta_bar_active;
    new_phi_path(:, 1) = phi_bar_active;

    cum_z = zeros(n_active, 1);
    cum_p = zeros(n_active, 1);
    for node = 2:n_nodes
        cum_z = cum_z + h * rhs_z_active(:, node - 1);
        cum_p = cum_p + h * rhs_p_active(:, node - 1);
        new_zeta_path(:, node) = zeta_bar_active - cum_z;
        new_phi_path(:, node) = phi_bar_active - cum_p;
    end

    zeta_path = new_zeta_path;
    phi_path = new_phi_path;
end

zeta_final = zeta_bar_hat;
phi_final = phi_bar_hat;
for node = 1:(n_nodes - 1)
    zeta_node = zeros(ny, nx);
    phi_node = zeros(ny, nx);
    zeta_node(active_idx) = zeta_path(:, node);
    phi_node(active_idx) = phi_path(:, node);

    [rhs_z_full, rhs_p_full] = local_coupled_lambda_rhs( ...
        zeta_node, phi_node, plan, ny, nx);

    zeta_final = zeta_final - h * rhs_z_full;
    phi_final = phi_final - h * rhs_p_full;
end
end

function [kvec, mvec] = local_fft_wavenumbers(n, d)
if n == 1
    mvec = 0;
    kvec = 0;
    return;
end
L = n * d;
if mod(n, 2) ~= 0
    error('directional_creamer_transform:EvenGridRequired', ...
        'Current helper assumes an even grid size. Got n=%d.', n);
end
mvec = local_fft_mode_numbers(n);
kvec = (2 * pi / L) * mvec;
end

function mvec = local_fft_mode_numbers(n)
if n == 1
    mvec = 0;
    return;
end
if mod(n, 2) ~= 0
    error('directional_creamer_transform:EvenGridRequired', ...
        'Current helper assumes an even grid size. Got n=%d.', n);
end
mvec = [0:(n/2), (-n/2 + 1):-1];
end

function idx = local_mode_to_index(m, n)
m_wrapped = mod(m + n/2, n) - n/2;
idx = m_wrapped + 1;
if idx < 1
    idx = idx + n;
end
end

function phi_hat = local_linear_phi_from_eta(zeta_hat, Theta, KX, KY, g, propagation_direction_deg)
phi_hat = zeros(size(zeta_hat));

dir_vec = [cosd(propagation_direction_deg), sind(propagation_direction_deg)];
selector = KX * dir_vec(1) + KY * dir_vec(2);
eps_dir = 1e-12;

positive_mask = selector > eps_dir;
tie_mask = abs(selector) <= eps_dir;
positive_mask = positive_mask | (tie_mask & (KY > eps_dir));
positive_mask = positive_mask | (tie_mask & abs(KY) <= eps_dir & KX > eps_dir);

nonzero_mask = Theta > 0;
positive_mask = positive_mask & nonzero_mask;

phase_factor = sqrt(g ./ Theta(positive_mask));
phi_hat(positive_mask) = -1i * phase_factor .* zeta_hat(positive_mask);

neg_rows = size(zeta_hat, 1);
neg_cols = size(zeta_hat, 2);
[rows, cols] = find(positive_mask);
for n = 1:numel(rows)
    row = rows(n);
    col = cols(n);
    row_neg = mod(neg_rows - row + 1, neg_rows) + 1;
    col_neg = mod(neg_cols - col + 1, neg_cols) + 1;
    phi_hat(row_neg, col_neg) = conj(phi_hat(row, col));
end

phi_hat(~nonzero_mask) = 0;
end

function d_val = local_kernel_D_deep(theta1, theta2, theta3, g)
den_base = theta1 * (theta2 + theta3 - theta1) ...
         + theta2 * (theta1 + theta3 - theta2) ...
         + theta3 * (theta1 + theta2 - theta3);

if abs(den_base) < 1e-14 || theta1 == 0 || theta2 == 0 || theta3 == 0
    d_val = 0;
    return;
end

num = (theta1 + theta2 + theta3) ...
    * (theta1 + theta2 - theta3) ...
    * (theta1 - theta2 - theta3) ...
    * (theta1 + theta3 - theta2);

d_val = num / (12 * g * den_base);
end

function b_val = local_kernel_B_deep(theta1, theta2, theta3)
den_base = theta1 * (theta2 + theta3 - theta1) ...
         + theta2 * (theta1 + theta3 - theta2) ...
         + theta3 * (theta1 + theta2 - theta3);

if abs(den_base) < 1e-14
    b_val = 0;
    return;
end

num = (theta1 - theta2)^2 * (theta1 + theta2) ...
    - theta3^2 * (2 * theta3 - theta1 - theta2);

b_val = num / (2 * den_base);
end

function d_val = local_kernel_D_finite(theta1, theta2, theta3, k1sq, k2sq, k3sq, g)
den_base = theta1 * (theta2 + theta3 - theta1) ...
         + theta2 * (theta1 + theta3 - theta2) ...
         + theta3 * (theta1 + theta2 - theta3);

if abs(den_base) < 1e-14 || theta1 == 0 || theta2 == 0 || theta3 == 0
    d_val = 0;
    return;
end

num = k1sq * (theta1^2 - (theta2 - theta3)^2) ...
    + k2sq * (theta2^2 - (theta1 - theta3)^2) ...
    + k3sq * (theta3^2 - (theta1 - theta2)^2) ...
    - 2 * theta1 * theta2 * theta3 * (theta1 + theta2 + theta3);

d_val = num / (12 * g * den_base);
end

function b_val = local_kernel_B_finite(theta1, theta2, theta3, k1sq, k2sq, k3sq)
den_base = theta1 * (theta2 + theta3 - theta1) ...
         + theta2 * (theta1 + theta3 - theta2) ...
         + theta3 * (theta1 + theta2 - theta3);

if abs(den_base) < 1e-14
    b_val = 0;
    return;
end

num = theta3 * (theta3 * (theta1 + theta2) - theta1^2 - theta2^2) ...
    + theta3 * (k1sq + k2sq - 2 * k3sq) ...
    + (theta1 - theta2) * (k1sq - k2sq);

b_val = num / (2 * den_base);
end
