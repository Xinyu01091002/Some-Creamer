function mf12 = mf12_from_linear_focus_1d(eta_lin, x_vec, cfg)
%MF12_FROM_LINEAR_FOCUS_1D Build MF12 eta20/eta22/eta33 from a 1D linear surface.

if nargin < 3 || isempty(cfg)
    cfg = struct();
end
cfg = local_apply_defaults(cfg);

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
mf12_matlab_dir = fullfile(repo_root, 'external', ...
    'spectral-domain-implementation-of-wave-interaction-theory', 'matlab');
setup_file = fullfile(mf12_matlab_dir, 'setup_paths.m');
if ~isfile(setup_file)
    error('mf12_from_linear_focus_1d:MissingMF12', ...
        'Could not find cloned MF12 MATLAB repo at %s', mf12_matlab_dir);
end
addpath(fullfile(mf12_matlab_dir, 'irregularWavesMF12', 'Source'));

eta_lin = reshape(real(eta_lin), 1, []);
x_vec = x_vec(:).';
nx = numel(eta_lin);

if numel(x_vec) ~= nx
    error('mf12_from_linear_focus_1d:SizeMismatch', ...
        'Length of x_vec (%d) must match numel(eta_lin) (%d).', numel(x_vec), nx);
end
if mod(nx, 2) ~= 0
    error('mf12_from_linear_focus_1d:EvenGridRequired', ...
        'Current helper assumes an even grid size. Got nx=%d.', nx);
end

dx = mean(diff(x_vec));
Lx = dx * nx;

mx = [0:(nx/2), (-nx/2 + 1):-1];
kx_fft = (2 * pi / Lx) * mx;
eta_hat = fft(eta_lin) / nx;
energy_density = abs(eta_hat).^2;

positive_mask = (mx > 0);
positive_idx = find(positive_mask);
positive_energy = energy_density(positive_idx);

[sorted_energy, order] = sort(positive_energy, 'descend');
energy_fraction = cumsum(sorted_energy) / sum(sorted_energy);
keep_count = find(energy_fraction >= cfg.energy_fraction, 1, 'first');
if isempty(keep_count)
    keep_count = numel(order);
end
if isfinite(cfg.max_active_modes)
    keep_count = min(keep_count, cfg.max_active_modes);
end

selected_idx = positive_idx(order(1:keep_count));
selected_coeffs = 2 * eta_hat(selected_idx);
a = real(selected_coeffs(:)).';
b = imag(selected_coeffs(:)).';
kx = kx_fft(selected_idx);
ky = zeros(size(kx));

c1 = mf12_spectral_coefficients(1, cfg.g, cfg.h, a, b, kx, ky, cfg.Ux, cfg.Uy);
c2 = mf12_spectral_coefficients(2, cfg.g, cfg.h, a, b, kx, ky, cfg.Ux, cfg.Uy);
c3 = mf12_spectral_coefficients(3, cfg.g, cfg.h, a, b, kx, ky, cfg.Ux, cfg.Uy);

eta1 = local_eta_component_1d('first', c1, Lx, nx, cfg.t);
eta20 = local_eta_component_1d('second_sub', c2, Lx, nx, cfg.t);
eta22 = local_eta_component_1d('second_super', c2, Lx, nx, cfg.t);
eta33 = local_eta_component_1d('third_super', c3, Lx, nx, cfg.t);
eta_total_order3 = eta1 + eta20 + eta22 + eta33;

mf12 = struct();
mf12.a = a;
mf12.b = b;
mf12.kx = kx;
mf12.active_mode_count = numel(a);
mf12.energy_fraction = cfg.energy_fraction;
mf12.h = cfg.h;
mf12.g = cfg.g;
mf12.Ux = cfg.Ux;
mf12.Uy = cfg.Uy;
mf12.t = cfg.t;
mf12.eta1 = eta1;
mf12.eta20 = eta20;
mf12.eta22 = eta22;
mf12.eta33 = eta33;
mf12.eta_total_order3 = eta_total_order3;
end

function cfg = local_apply_defaults(cfg)
defaults = struct( ...
    'g', 9.81, ...
    'h', 35, ...
    'Ux', 0, ...
    'Uy', 0, ...
    't', 0, ...
    'energy_fraction', 0.9999999999, ...
    'max_active_modes', inf);

names = fieldnames(defaults);
for n = 1:numel(names)
    name = names{n};
    if ~isfield(cfg, name) || isempty(cfg.(name))
        cfg.(name) = defaults.(name);
    end
end
end

function eta = local_eta_component_1d(component, coeffs, Lx, Nx, t)
dkx = 2*pi/Lx;
spec_eta = complex(zeros(1, Nx));

switch component
    case 'first'
        Z = (coeffs.a(:) + 1i*coeffs.b(:)) .* exp(-1i * coeffs.omega(:) * t);
        spec_eta = deposit_1d(spec_eta, coeffs.kx(:), Z, dkx, Nx);

    case 'second_super'
        Z2 = (coeffs.A_2(:) + 1i*coeffs.B_2(:)) .* exp(-1i * (2 * coeffs.omega(:)) * t);
        spec_eta = deposit_1d(spec_eta, 2 * coeffs.kx(:), Z2 .* coeffs.G_2(:), dkx, Nx);
        if isfield(coeffs, 'kx_npm')
            plus_mask = false(numel(coeffs.kx_npm), 1);
            plus_mask(1:2:end) = true;
            spec_eta = deposit_1d(spec_eta, coeffs.kx_npm(plus_mask), ...
                ((coeffs.A_npm(plus_mask) + 1i*coeffs.B_npm(plus_mask)) .* exp(-1i * coeffs.omega_npm(plus_mask) * t)) .* coeffs.G_npm(plus_mask), ...
                dkx, Nx);
        end

    case 'second_sub'
        if isfield(coeffs, 'kx_npm')
            minus_mask = false(numel(coeffs.kx_npm), 1);
            minus_mask(2:2:end) = true;
            kx_vals = coeffs.kx_npm(minus_mask);
            Z = ((coeffs.A_npm(minus_mask) + 1i*coeffs.B_npm(minus_mask)) .* exp(-1i * coeffs.omega_npm(minus_mask) * t)) .* coeffs.G_npm(minus_mask);
            spec_eta = deposit_1d(spec_eta, kx_vals, Z, dkx, Nx);
        end

    case 'third_super'
        if isfield(coeffs, 'G_3')
            Z3 = (coeffs.A_3(:) + 1i*coeffs.B_3(:)) .* exp(-1i * (3 * coeffs.omega(:)) * t);
            spec_eta = deposit_1d(spec_eta, 3 * coeffs.kx(:), Z3 .* coeffs.G_3(:), dkx, Nx);
        end
        if isfield(coeffs, 'G_np2m')
            Z = (coeffs.A_np2m(:) + 1i*coeffs.B_np2m(:)) .* exp(-1i * coeffs.omega_np2m(:) * t);
            spec_eta = deposit_1d(spec_eta, coeffs.kx_np2m(:), Z .* coeffs.G_np2m(:), dkx, Nx);
        end
        if isfield(coeffs, 'G_2npm')
            Z = (coeffs.A_2npm(:) + 1i*coeffs.B_2npm(:)) .* exp(-1i * coeffs.omega_2npm(:) * t);
            spec_eta = deposit_1d(spec_eta, coeffs.kx_2npm(:), Z .* coeffs.G_2npm(:), dkx, Nx);
        end
        if isfield(coeffs, 'G_npmpp')
            Z = 2 * (coeffs.A_npmpp(:) + 1i*coeffs.B_npmpp(:)) .* exp(-1i * coeffs.omega_npmpp(:) * t);
            spec_eta = deposit_1d(spec_eta, coeffs.kx_npmpp(:), Z .* coeffs.G_npmpp(:), dkx, Nx);
        end

    otherwise
        error('mf12_from_linear_focus_1d:UnknownComponent', ...
            'Unknown MF12 component: %s', component);
end

eta = real(ifft(spec_eta)) * Nx;
end

function spec = deposit_1d(spec, kx_in, values, dkx, Nx)
ux = (kx_in(:) / dkx);
vals = values(:);

valid = isfinite(ux) & isfinite(vals) & (abs(kx_in(:)) > 0);
ux = ux(valid);
vals = vals(valid);
if isempty(vals)
    return;
end

ix0 = floor(ux);
fx = ux - ix0;
tol = 1e-12;
fx(abs(fx) < tol) = 0;
fx(abs(fx-1) < tol) = 1;

ix1 = ix0 + 1;
idx_x0 = mod(ix0, Nx) + 1;
idx_x1 = mod(ix1, Nx) + 1;

w0 = 1 - fx;
w1 = fx;

spec = accum_add_1d(spec, idx_x0, vals .* w0);
spec = accum_add_1d(spec, idx_x1, vals .* w1);
end

function spec = accum_add_1d(spec, ix, vals)
spec = spec + reshape(accumarray(ix, vals, [numel(spec), 1], @sum, 0), 1, []);
end
