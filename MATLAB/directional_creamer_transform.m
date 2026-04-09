function [eta_nl, diagnostics, phi_nl] = directional_creamer_transform(eta_lin, x_vec, y_vec, cfg)
% DIRECTIONAL_CREAMER_TRANSFORM
% First-pass directional Creamer prototype for deep-water waves.
%
% This implementation uses the general 2D deep-water kernels B and D from
% Creamer et al. (1989), equations (3.4), (3.5), and the inverse map
% implied by solving the Lie-transform equations from lambda = 1 to
% lambda = 0, truncated at the first non-trivial order:
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
% Fourier space using a positive-frequency deep-water convention aligned
% with cfg.propagation_direction_deg. This is the main additional closure
% assumption in the current prototype.

if nargin < 4 || isempty(cfg)
    cfg = struct();
end

cfg = local_apply_defaults(cfg);

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
dy = mean(diff(y_vec));
Lx = dx * nx;
Ly = dy * ny;
n_total = nx * ny;

[kx, mx] = local_fft_wavenumbers(nx, dx);
[ky, my] = local_fft_wavenumbers(ny, dy);
[KX, KY] = meshgrid(kx, ky);
Kmag = hypot(KX, KY);

zeta_hat = fft2(eta_lin) / n_total;
phi_bar_hat = local_linear_phi_from_eta(zeta_hat, Kmag, KX, KY, cfg.g, cfg.propagation_direction_deg);

if ~cfg.preserve_mean
    zeta_hat(1, 1) = 0;
    phi_bar_hat(1, 1) = 0;
end

active_mask = local_build_active_mask(zeta_hat, cfg.energy_fraction, cfg.max_active_modes, cfg.preserve_mean);
n_active = nnz(active_mask);

if cfg.n_lambda_steps <= 0
    delta_zeta_hat = local_inverse_zeta_correction(zeta_hat, phi_bar_hat, active_mask, Kmag, kx, mx, ny, nx, cfg);
    eta_nl_hat = zeta_hat + delta_zeta_hat;
    phi_nl_hat = phi_bar_hat + local_inverse_phi_correction(zeta_hat, phi_bar_hat, active_mask, Kmag, kx, mx, ny, nx);
else
    zeta_state = zeta_hat;
    phi_state = phi_bar_hat;
    h_lambda = -1 / cfg.n_lambda_steps;
    for step = 1:cfg.n_lambda_steps
        if strcmpi(cfg.lambda_flow_model, 'canonical_pair')
            [k1z, k1p] = local_coupled_lambda_rhs(zeta_state, phi_state, active_mask, Kmag, kx, mx, ny, nx, cfg);
            [k2z, k2p] = local_coupled_lambda_rhs( ...
                zeta_state + 0.5 * h_lambda * k1z, ...
                phi_state + 0.5 * h_lambda * k1p, ...
                active_mask, Kmag, kx, mx, ny, nx, cfg);
            [k3z, k3p] = local_coupled_lambda_rhs( ...
                zeta_state + 0.5 * h_lambda * k2z, ...
                phi_state + 0.5 * h_lambda * k2p, ...
                active_mask, Kmag, kx, mx, ny, nx, cfg);
            [k4z, k4p] = local_coupled_lambda_rhs( ...
                zeta_state + h_lambda * k3z, ...
                phi_state + h_lambda * k3p, ...
                active_mask, Kmag, kx, mx, ny, nx, cfg);

            zeta_state = zeta_state + (h_lambda / 6) * (k1z + 2*k2z + 2*k3z + k4z);
            phi_state = phi_state + (h_lambda / 6) * (k1p + 2*k2p + 2*k3p + k4p);
        elseif strcmpi(cfg.lambda_flow_model, 'backward_picard_316')
            [zeta_state, phi_state] = local_backward_picard_316( ...
                zeta_hat, phi_bar_hat, active_mask, Kmag, kx, mx, ny, nx, cfg);
            break;
        else
            k1z = local_legacy_lambda_rhs(zeta_state, active_mask, Kmag, KX, KY, kx, mx, ny, nx, cfg);
            k2z = local_legacy_lambda_rhs(zeta_state + 0.5 * h_lambda * k1z, active_mask, Kmag, KX, KY, kx, mx, ny, nx, cfg);
            k3z = local_legacy_lambda_rhs(zeta_state + 0.5 * h_lambda * k2z, active_mask, Kmag, KX, KY, kx, mx, ny, nx, cfg);
            k4z = local_legacy_lambda_rhs(zeta_state + h_lambda * k3z, active_mask, Kmag, KX, KY, kx, mx, ny, nx, cfg);
            zeta_state = zeta_state + (h_lambda / 6) * (k1z + 2*k2z + 2*k3z + k4z);
            phi_state = local_linear_phi_from_eta(zeta_state, Kmag, KX, KY, cfg.g, cfg.propagation_direction_deg);
        end

        if ~cfg.preserve_mean
            zeta_state(1,1) = 0;
            phi_state(1,1) = 0;
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
diagnostics.max_eta_lin = max(abs(eta_lin(:)));
diagnostics.max_eta_nl = max(abs(eta_nl(:)));
diagnostics.max_imag_eta_nl = max(abs(imag(ifft2(eta_nl_hat * n_total))), [], 'all');
diagnostics.max_imag_phi_nl = max(abs(imag(ifft2(phi_nl_hat * n_total))), [], 'all');
diagnostics.max_abs_phi_hat = max(abs(phi_nl_hat(:)));
diagnostics.input_mean = mean(eta_lin(:));
diagnostics.output_mean = mean(eta_nl(:));
diagnostics.propagation_direction_deg = cfg.propagation_direction_deg;

if cfg.verbose
    fprintf('Directional Creamer prototype summary:\n');
    fprintf('  Grid                : ny=%d, nx=%d\n', ny, nx);
    fprintf('  Spacing             : dx=%.6g, dy=%.6g\n', dx, dy);
    fprintf('  Active modes        : %d\n', n_active);
    fprintf('  Lambda RK steps     : %d\n', cfg.n_lambda_steps);
    fprintf('  Lambda flow model   : %s\n', cfg.lambda_flow_model);
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
    'max_active_modes', inf, ...
    'n_lambda_steps', 0, ...
    'n_picard_iters', 4, ...
    'lambda_flow_model', 'canonical_pair', ...
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

function active_mask = local_build_active_mask(zeta_hat, energy_fraction, max_active_modes, preserve_mean)
energy_density = abs(zeta_hat).^2;
flat_energy = energy_density(:);
[sorted_energy, order] = sort(flat_energy, 'descend');
cumulative_fraction = cumsum(sorted_energy) / sum(sorted_energy);

keep_count = find(cumulative_fraction >= energy_fraction, 1, 'first');
if isempty(keep_count)
    keep_count = numel(order);
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

function delta_zeta_hat = local_inverse_zeta_correction(zeta_hat, phi_hat, active_mask, Kmag, kx, mx, ny, nx, cfg)
[row_idx, col_idx] = find(active_mask);
active_idx = sub2ind([ny, nx], row_idx, col_idx);
n_active = numel(active_idx);

mode_struct = struct();
mode_struct.mx = mx(col_idx).';
my = [0:(ny/2), (-ny/2 + 1):-1];
mode_struct.my = my(row_idx).';
mode_struct.kmag = Kmag(active_idx);
mode_struct.zeta = zeta_hat(active_idx);
mode_struct.phi = phi_hat(active_idx);

delta_zeta_hat = zeros(ny, nx);
for a = 1:n_active
    mx_a = mode_struct.mx(a);
    my_a = mode_struct.my(a);
    kmag_a = mode_struct.kmag(a);
    zeta_a = mode_struct.zeta(a);
    phi_a = mode_struct.phi(a);

    for b = 1:n_active
        mx_b = mode_struct.mx(b);
        my_b = mode_struct.my(b);
        kmag_b = mode_struct.kmag(b);
        zeta_b = mode_struct.zeta(b);
        phi_b = mode_struct.phi(b);

        mx_k = mx_a + mx_b;
        my_k = my_a + my_b;

        col_k = local_mode_to_index(mx_k, nx);
        row_k = local_mode_to_index(my_k, ny);

        kmag_k = Kmag(row_k, col_k);

        d_term = 0;
        if kmag_a > 0 && kmag_b > 0 && abs(phi_a) > 0 && abs(phi_b) > 0
            d_val = local_kernel_D(kmag_k, kmag_a, kmag_b, cfg.g);
            d_term = 3 * d_val * phi_a * phi_b;
        end

        b_val = local_kernel_B(kmag_a, kmag_b, kmag_k);
        b_term = b_val * zeta_a * zeta_b;

        delta_zeta_hat(row_k, col_k) = delta_zeta_hat(row_k, col_k) - (d_term + b_term);
    end
end
end

function delta_phi_hat = local_inverse_phi_correction(zeta_hat, phi_hat, active_mask, Kmag, kx, mx, ny, nx)
[row_idx, col_idx] = find(active_mask);
active_idx = sub2ind([ny, nx], row_idx, col_idx);
n_active = numel(active_idx);

mode_struct = struct();
mode_struct.mx = mx(col_idx).';
my = [0:(ny/2), (-ny/2 + 1):-1];
mode_struct.my = my(row_idx).';
mode_struct.kmag = Kmag(active_idx);
mode_struct.zeta = zeta_hat(active_idx);
mode_struct.phi = phi_hat(active_idx);

delta_phi_hat = zeros(ny, nx);
for a = 1:n_active
    mx_a = mode_struct.mx(a);
    my_a = mode_struct.my(a);
    zeta_a = mode_struct.zeta(a);

    for b = 1:n_active
        mx_b = mode_struct.mx(b);
        my_b = mode_struct.my(b);
        phi_b = mode_struct.phi(b);

        mx_k = mx_a + mx_b;
        my_k = my_a + my_b;

        col_k = local_mode_to_index(mx_k, nx);
        row_k = local_mode_to_index(my_k, ny);
        kmag_k = Kmag(row_k, col_k);
        kmag_a = mode_struct.kmag(a);
        kmag_b = mode_struct.kmag(b);

        b_val = local_kernel_B(kmag_k, kmag_a, kmag_b);
        delta_phi_hat(row_k, col_k) = delta_phi_hat(row_k, col_k) + 2 * b_val * zeta_a * phi_b;
    end
end
end

function delta_zeta_hat = local_legacy_lambda_rhs(zeta_hat, active_mask, Kmag, KX, KY, kx, mx, ny, nx, cfg)
phi_hat = local_linear_phi_from_eta(zeta_hat, Kmag, KX, KY, cfg.g, cfg.propagation_direction_deg);

[row_idx, col_idx] = find(active_mask);
active_idx = sub2ind([ny, nx], row_idx, col_idx);
n_active = numel(active_idx);

mode_struct = struct();
mode_struct.mx = mx(col_idx).';
my = [0:(ny/2), (-ny/2 + 1):-1];
mode_struct.my = my(row_idx).';
mode_struct.kmag = Kmag(active_idx);
mode_struct.zeta = zeta_hat(active_idx);
mode_struct.phi = phi_hat(active_idx);

delta_zeta_hat = zeros(ny, nx);
for a = 1:n_active
    mx_a = mode_struct.mx(a);
    my_a = mode_struct.my(a);
    kmag_a = mode_struct.kmag(a);
    zeta_a = mode_struct.zeta(a);
    phi_a = mode_struct.phi(a);

    for b = 1:n_active
        mx_b = mode_struct.mx(b);
        my_b = mode_struct.my(b);
        kmag_b = mode_struct.kmag(b);
        zeta_b = mode_struct.zeta(b);
        phi_b = mode_struct.phi(b);

        mx_k = mx_a + mx_b;
        my_k = my_a + my_b;

        col_k = local_mode_to_index(mx_k, nx);
        row_k = local_mode_to_index(my_k, ny);

        kmag_k = Kmag(row_k, col_k);

        d_term = 0;
        if kmag_a > 0 && kmag_b > 0 && abs(phi_a) > 0 && abs(phi_b) > 0
            d_val = local_kernel_D(kmag_k, kmag_a, kmag_b, cfg.g);
            d_term = 3 * d_val * phi_a * phi_b;
        end

        b_val = local_kernel_B(kmag_a, kmag_b, kmag_k);
        b_term = b_val * zeta_a * zeta_b;

        delta_zeta_hat(row_k, col_k) = delta_zeta_hat(row_k, col_k) - (d_term + b_term);
    end
end
end

function [dzeta_hat, dphi_hat] = local_coupled_lambda_rhs(zeta_hat, phi_hat, active_mask, Kmag, kx, mx, ny, nx, cfg)
[row_idx, col_idx] = find(active_mask);
active_idx = sub2ind([ny, nx], row_idx, col_idx);
n_active = numel(active_idx);

mode_struct = struct();
mode_struct.mx = mx(col_idx).';
my = [0:(ny/2), (-ny/2 + 1):-1];
mode_struct.my = my(row_idx).';
mode_struct.kmag = Kmag(active_idx);
mode_struct.zeta = zeta_hat(active_idx);
mode_struct.phi = phi_hat(active_idx);

dzeta_hat = zeros(ny, nx);
dphi_hat = zeros(ny, nx);

for a = 1:n_active
    mx_a = mode_struct.mx(a);
    my_a = mode_struct.my(a);
    kmag_a = mode_struct.kmag(a);
    zeta_a = mode_struct.zeta(a);
    phi_a = mode_struct.phi(a);

    for b = 1:n_active
        mx_b = mode_struct.mx(b);
        my_b = mode_struct.my(b);
        kmag_b = mode_struct.kmag(b);
        zeta_b = mode_struct.zeta(b);
        phi_b = mode_struct.phi(b);

        mx_k = mx_a + mx_b;
        my_k = my_a + my_b;

        col_k = local_mode_to_index(mx_k, nx);
        row_k = local_mode_to_index(my_k, ny);
        kmag_k = Kmag(row_k, col_k);

        d_term = 0;
        if kmag_a > 0 && kmag_b > 0 && abs(phi_a) > 0 && abs(phi_b) > 0
            d_val = local_kernel_D(kmag_k, kmag_a, kmag_b, cfg.g);
            d_term = 3 * d_val * phi_a * phi_b;
        end

        b_val_z = local_kernel_B(kmag_a, kmag_b, kmag_k);
        dzeta_hat(row_k, col_k) = dzeta_hat(row_k, col_k) + d_term + b_val_z * zeta_a * zeta_b;

        b_val_phi = local_kernel_B(kmag_k, kmag_a, kmag_b);
        dphi_hat(row_k, col_k) = dphi_hat(row_k, col_k) - 2 * b_val_phi * zeta_a * phi_b;
    end
end
end

function [zeta_final, phi_final] = local_backward_picard_316(zeta_bar_hat, phi_bar_hat, active_mask, Kmag, kx, mx, ny, nx, cfg)
active_idx = find(active_mask(:));
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
            zeta_node, phi_node, active_mask, Kmag, kx, mx, ny, nx, cfg);

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
        zeta_node, phi_node, active_mask, Kmag, kx, mx, ny, nx, cfg);

    zeta_final = zeta_final - h * rhs_z_full;
    phi_final = phi_final - h * rhs_p_full;
end
end

function [kvec, mvec] = local_fft_wavenumbers(n, d)
L = n * d;
if mod(n, 2) ~= 0
    error('directional_creamer_transform:EvenGridRequired', ...
        'Current helper assumes an even grid size. Got n=%d.', n);
end
mvec = [0:(n/2), (-n/2 + 1):-1];
kvec = (2 * pi / L) * mvec;
end

function idx = local_mode_to_index(m, n)
m_wrapped = mod(m + n/2, n) - n/2;
idx = m_wrapped + 1;
if idx < 1
    idx = idx + n;
end
end

function phi_hat = local_linear_phi_from_eta(zeta_hat, Kmag, KX, KY, g, propagation_direction_deg)
phi_hat = zeros(size(zeta_hat));

dir_vec = [cosd(propagation_direction_deg), sind(propagation_direction_deg)];
selector = KX * dir_vec(1) + KY * dir_vec(2);
eps_dir = 1e-12;

positive_mask = selector > eps_dir;
tie_mask = abs(selector) <= eps_dir;
positive_mask = positive_mask | (tie_mask & (KY > eps_dir));
positive_mask = positive_mask | (tie_mask & abs(KY) <= eps_dir & KX > eps_dir);

nonzero_mask = Kmag > 0;
positive_mask = positive_mask & nonzero_mask;

phase_factor = sqrt(g ./ Kmag(positive_mask));
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

function d_val = local_kernel_D(theta1, theta2, theta3, g)
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

function b_val = local_kernel_B(theta1, theta2, theta3)
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
