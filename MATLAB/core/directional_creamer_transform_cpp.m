function [eta_nl, diagnostics, phi_nl, spectra] = directional_creamer_transform_cpp(eta_lin, x_vec, y_vec, cfg)
%DIRECTIONAL_CREAMER_TRANSFORM_CPP Use the standalone C++ Creamer flow CLI.
% MATLAB prepares the initial spectra; C++ selects active modes, builds the
% canonical-pair interaction plan, runs fixed RK4 lambda-flow, and returns
% final spectra.

if nargin < 4 || isempty(cfg)
    cfg = struct();
end
cfg = local_defaults(cfg);

[ny, nx] = size(eta_lin);
x_vec = x_vec(:).';
y_vec = y_vec(:).';
dx = mean(diff(x_vec));
if ny > 1
    dy = mean(diff(y_vec));
else
    dy = 1;
end
n_total = nx * ny;

[kx, ~] = local_fft_wavenumbers(nx, dx);
[ky, ~] = local_fft_wavenumbers(ny, dy);
[KX, KY] = meshgrid(kx, ky);
Kmag = hypot(KX, KY);

zeta_hat = fft2(eta_lin) / n_total;
phi_hat = local_linear_phi_from_eta(zeta_hat, Kmag, KX, KY, cfg.g, cfg.propagation_direction_deg);

if ~cfg.preserve_mean
    zeta_hat(1, 1) = 0;
    phi_hat(1, 1) = 0;
end

job_dir = char(cfg.cpp_job_dir);
if ~exist(job_dir, 'dir')
    mkdir(job_dir);
end
local_write_job(job_dir, zeta_hat, phi_hat, ny, nx, dx, dy, cfg);

exe = char(cfg.cpp_exe);
if ~isfile(exe)
    error('directional_creamer_transform_cpp:MissingExecutable', ...
        'C++ executable not found: %s', exe);
end

cmd = sprintf('"%s" "%s"', exe, job_dir);
t_cpp = tic;
[status, output] = system(cmd);
cpp_wall_s = toc(t_cpp);
if status ~= 0
    error('directional_creamer_transform_cpp:CppFailed', ...
        'C++ Creamer flow failed with status %d:\n%s', status, output);
end
cpp_summary = local_read_summary(fullfile(job_dir, 'summary.txt'));

eta_nl_hat = reshape(local_read_complex(fullfile(job_dir, 'zeta_final.bin'), n_total), ny, nx);
phi_nl_hat = reshape(local_read_complex(fullfile(job_dir, 'phi_final.bin'), n_total), ny, nx);

eta_nl = real(ifft2(eta_nl_hat * n_total));
phi_nl = real(ifft2(phi_nl_hat * n_total));

diagnostics = struct();
diagnostics.nx = nx;
diagnostics.ny = ny;
diagnostics.dx = dx;
diagnostics.dy = dy;
diagnostics.active_mode_count = cpp_summary.n_active;
diagnostics.n_lambda_steps = cfg.n_lambda_steps;
diagnostics.lambda_flow_model = 'canonical_pair';
diagnostics.lambda_stepper = 'cpp_fixed_rk4';
diagnostics.cpp_wall_s = cpp_wall_s;
diagnostics.cpp_kernel_runtime_s = cpp_summary.runtime_s;
diagnostics.cpp_openmp_threads = cpp_summary.openmp_threads;
diagnostics.max_eta_lin = max(abs(eta_lin(:)));
diagnostics.max_eta_nl = max(abs(eta_nl(:)));
diagnostics.max_imag_eta_nl = max(abs(imag(ifft2(eta_nl_hat * n_total))), [], 'all');
diagnostics.max_imag_phi_nl = max(abs(imag(ifft2(phi_nl_hat * n_total))), [], 'all');

spectra = struct('zeta_hat', eta_nl_hat, 'phi_hat', phi_nl_hat, ...
    'zeta0_hat', zeta_hat, 'phi0_hat', phi_hat);
end

function cfg = local_defaults(cfg)
defaults = struct( ...
    'g', 9.81, ...
    'energy_fraction', 0.99, ...
    'min_active_modes', 0, ...
    'max_active_modes', inf, ...
    'n_lambda_steps', 6, ...
    'propagation_direction_deg', 0, ...
    'preserve_mean', true, ...
    'cpp_exe', fullfile(pwd, '..', '..', 'cpp', 'creamer_flow', 'build', 'creamer_flow_plan.exe'), ...
    'cpp_job_dir', fullfile(tempdir, 'creamer_cpp_job'));
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(cfg, names{i}) || isempty(cfg.(names{i}))
        cfg.(names{i}) = defaults.(names{i});
    end
end
end

function local_write_job(job_dir, zeta_hat, phi_hat, ny, nx, dx, dy, cfg)
n_total = nx * ny;
local_write_int64(fullfile(job_dir, 'meta.bin'), ...
    int64([ny nx n_total cfg.n_lambda_steps double(cfg.preserve_mean)]));
local_write_double(fullfile(job_dir, 'params.bin'), ...
    [dx dy cfg.g cfg.energy_fraction cfg.min_active_modes cfg.max_active_modes]);
local_write_complex(fullfile(job_dir, 'zeta0.bin'), zeta_hat(:));
local_write_complex(fullfile(job_dir, 'phi0.bin'), phi_hat(:));
end

function local_write_int64(path, data)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, data, 'int64');
end

function local_write_double(path, data)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, data, 'double');
end

function local_write_complex(path, z)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, [real(z(:)).'; imag(z(:)).'], 'double');
end

function z = local_read_complex(path, n)
fid = fopen(path, 'r');
cleanup = onCleanup(@() fclose(fid));
raw = fread(fid, [2, n], 'double');
z = complex(raw(1, :).', raw(2, :).');
end

function [kvec, mvec] = local_fft_wavenumbers(n, d)
if n == 1
    mvec = 0;
    kvec = 0;
    return;
end
L = n * d;
mvec = local_fft_mode_numbers(n);
kvec = (2 * pi / L) * mvec;
end

function mvec = local_fft_mode_numbers(n)
if n == 1
    mvec = 0;
else
    mvec = [0:(n/2), (-n/2 + 1):-1];
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
phi_hat(positive_mask) = -1i * sqrt(g ./ Kmag(positive_mask)) .* zeta_hat(positive_mask);

[rows, cols] = find(positive_mask);
for n = 1:numel(rows)
    row_neg = mod(size(zeta_hat, 1) - rows(n) + 1, size(zeta_hat, 1)) + 1;
    col_neg = mod(size(zeta_hat, 2) - cols(n) + 1, size(zeta_hat, 2)) + 1;
    phi_hat(row_neg, col_neg) = conj(phi_hat(rows(n), cols(n)));
end
phi_hat(~nonzero_mask) = 0;
end

function summary = local_read_summary(path)
summary = struct('n_active', NaN, 'runtime_s', NaN, 'openmp_threads', NaN);
if ~isfile(path)
    return;
end
lines = splitlines(string(fileread(path)));
for i = 1:numel(lines)
    parts = split(strtrim(lines(i)));
    if numel(parts) ~= 2
        continue;
    end
    key = char(parts(1));
    val = str2double(parts(2));
    if isfield(summary, key)
        summary.(key) = val;
    end
end
end
