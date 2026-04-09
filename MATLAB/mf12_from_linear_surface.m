function mf12 = mf12_from_linear_surface(eta_lin, x_vec, y_vec, cfg)
%MF12_FROM_LINEAR_SURFACE Build MF12 second-order fields from a linear surface.
%
% This helper bridges the current Creamer workflow to the cloned MF12 MATLAB
% implementation. It:
%   1. extracts a positive-half-plane modal set (a,b,kx,ky) from eta_lin
%   2. calls mf12_spectral_coefficients / mf12_spectral_surface
%   3. reconstructs first-order, second-superharmonic, second-subharmonic,
%      and total second-order surface fields on the same sample grid

if nargin < 4 || isempty(cfg)
    cfg = struct();
end
cfg = local_apply_defaults(cfg);

repo_root = fileparts(fileparts(mfilename('fullpath')));
mf12_matlab_dir = fullfile(repo_root, 'external', ...
    'spectral-domain-implementation-of-wave-interaction-theory', 'matlab');

setup_file = fullfile(mf12_matlab_dir, 'setup_paths.m');
if ~isfile(setup_file)
    error('mf12_from_linear_surface:MissingMF12', ...
        'Could not find cloned MF12 MATLAB repo at %s', mf12_matlab_dir);
end

run(setup_file);

[ny, nx] = size(eta_lin);
x_vec = x_vec(:).';
y_vec = y_vec(:).';
dx = mean(diff(x_vec));
dy = mean(diff(y_vec));
Lx = dx * nx;
Ly = dy * ny;
n_total = nx * ny;

[kx_fft, mx] = local_fft_wavenumbers(nx, dx);
[ky_fft, my] = local_fft_wavenumbers(ny, dy);
[KX, KY] = meshgrid(kx_fft, ky_fft);

eta_hat = fft2(eta_lin) / n_total;
energy_density = abs(eta_hat).^2;

positive_mask = (KX > 0) | (KX == 0 & KY > 0);
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

selected_linear_idx = positive_idx(order(1:keep_count));

% eta_hat stores the complex coefficient of the full Fourier expansion.
% The MF12 a/b convention describes the real surface using only one member
% of each conjugate pair:
%   eta = a cos(theta) + b sin(theta)
% Therefore, when we keep only the positive half plane, we must fold the
% conjugate partner into the retained mode and multiply the complex
% amplitude by 2.
selected_coeffs = 2 * eta_hat(selected_linear_idx);
selected_kx = KX(selected_linear_idx);
selected_ky = KY(selected_linear_idx);

a = real(selected_coeffs(:)).';
b = imag(selected_coeffs(:)).';
kx = selected_kx(:).';
ky = selected_ky(:).';

c1 = mf12_spectral_coefficients(1, cfg.g, cfg.h, a, b, kx, ky, cfg.Ux, cfg.Uy);
c2 = mf12_spectral_coefficients(2, cfg.g, cfg.h, a, b, kx, ky, cfg.Ux, cfg.Uy);

eta1 = mf12_eta_component_spectral('first', c1, Lx, Ly, nx, ny, cfg.t);
eta20 = mf12_eta_component_spectral('second_sub', c2, Lx, Ly, nx, ny, cfg.t);
eta22 = mf12_eta_component_spectral('second_super', c2, Lx, Ly, nx, ny, cfg.t);
eta2_total = eta20 + eta22;
eta_total_order2 = eta1 + eta2_total;

mf12 = struct();
mf12.a = a;
mf12.b = b;
mf12.kx = kx;
mf12.ky = ky;
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
mf12.eta2_total = eta2_total;
mf12.eta_total_order2 = eta_total_order2;
end

function cfg = local_apply_defaults(cfg)
defaults = struct( ...
    'g', 9.81, ...
    'h', [], ...
    'Ux', 0, ...
    'Uy', 0, ...
    't', 0, ...
    'energy_fraction', 0.99999, ...
    'max_active_modes', inf);

names = fieldnames(defaults);
for n = 1:numel(names)
    name = names{n};
    if ~isfield(cfg, name) || isempty(cfg.(name))
        cfg.(name) = defaults.(name);
    end
end

if isempty(cfg.h)
    error('mf12_from_linear_surface:MissingDepth', ...
        'cfg.h must be supplied. For the current deep-water comparison use h = 5/kp.');
end
end

function [kvec, mvec] = local_fft_wavenumbers(n, d)
L = n * d;
mvec = [0:(n/2), (-n/2 + 1):-1];
kvec = (2 * pi / L) * mvec;
end

function eta = mf12_eta_component_spectral(component, coeffs, Lx, Ly, Nx, Ny, t)
dkx = 2*pi/Lx;
dky = 2*pi/Ly;
spec_eta = complex(zeros(Ny, Nx));
N = coeffs.N;

switch component
    case 'first'
        Z = (coeffs.a(:) + 1i*coeffs.b(:)) .* exp(-1i * coeffs.omega(:) * t);
        spec_eta = deposit(spec_eta, coeffs.kx(:), coeffs.ky(:), Z, dkx, dky, Nx, Ny);

    case 'second_super'
        Z2 = (coeffs.A_2(:) + 1i*coeffs.B_2(:)) .* exp(-1i * (2*coeffs.omega(:)) * t);
        spec_eta = deposit(spec_eta, 2*coeffs.kx(:), 2*coeffs.ky(:), Z2 .* coeffs.G_2(:), dkx, dky, Nx, Ny);
        cnm = 0;
        for n = 1:N
            for m = n+1:N
                for pm = [1 -1]
                    cnm = cnm + 1;
                    if pm ~= 1
                        continue;
                    end
                    kx_eff = coeffs.kx(n) + pm*coeffs.kx(m);
                    ky_eff = coeffs.ky(n) + pm*coeffs.ky(m);
                    om_eff = coeffs.omega(n) + pm*coeffs.omega(m);
                    Z = (coeffs.A_npm(cnm) + 1i*coeffs.B_npm(cnm)) * exp(-1i * om_eff * t);
                    spec_eta = deposit(spec_eta, kx_eff, ky_eff, Z * coeffs.G_npm(cnm), dkx, dky, Nx, Ny);
                end
            end
        end

    case 'second_sub'
        cnm = 0;
        for n = 1:N
            for m = n+1:N
                for pm = [1 -1]
                    cnm = cnm + 1;
                    if pm ~= -1
                        continue;
                    end
                    kx_eff = coeffs.kx(n) + pm*coeffs.kx(m);
                    ky_eff = coeffs.ky(n) + pm*coeffs.ky(m);
                    om_eff = coeffs.omega(n) + pm*coeffs.omega(m);
                    Z = (coeffs.A_npm(cnm) + 1i*coeffs.B_npm(cnm)) * exp(-1i * om_eff * t);
                    spec_eta = deposit(spec_eta, kx_eff, ky_eff, Z * coeffs.G_npm(cnm), dkx, dky, Nx, Ny);
                end
            end
        end

    otherwise
        error('mf12_from_linear_surface:UnknownComponent', ...
            'Unknown MF12 component: %s', component);
end

eta = real(ifft2(spec_eta)) * (Nx * Ny);
end

function spec = deposit(spec, kx_in, ky_in, values, dkx, dky, Nx, Ny)
ux = (kx_in(:) / dkx);
uy = (ky_in(:) / dky);
vals = values(:);

valid = isfinite(ux) & isfinite(uy) & isfinite(vals);
ux = ux(valid);
uy = uy(valid);
vals = vals(valid);
if isempty(vals)
    return;
end

ix0 = floor(ux);
iy0 = floor(uy);
fx = ux - ix0;
fy = uy - iy0;
tol = 1e-12;
fx(abs(fx) < tol) = 0;
fy(abs(fy) < tol) = 0;
fx(abs(fx-1) < tol) = 1;
fy(abs(fy-1) < tol) = 1;

ix1 = ix0 + 1;
iy1 = iy0 + 1;
idx_x00 = mod(ix0, Nx) + 1;
idx_y00 = mod(iy0, Ny) + 1;
idx_x10 = mod(ix1, Nx) + 1;
idx_y10 = idx_y00;
idx_x01 = idx_x00;
idx_y01 = mod(iy1, Ny) + 1;
idx_x11 = idx_x10;
idx_y11 = idx_y01;

w00 = (1-fx).*(1-fy);
w10 = fx.*(1-fy);
w01 = (1-fx).*fy;
w11 = fx.*fy;

spec = accum_add(spec, idx_y00, idx_x00, vals.*w00);
spec = accum_add(spec, idx_y10, idx_x10, vals.*w10);
spec = accum_add(spec, idx_y01, idx_x01, vals.*w01);
spec = accum_add(spec, idx_y11, idx_x11, vals.*w11);
end

function spec = accum_add(spec, iy, ix, vals)
lin = sub2ind(size(spec), iy, ix);
spec = spec + reshape(accumarray(lin, vals, [numel(spec), 1], @sum, 0), size(spec));
end
