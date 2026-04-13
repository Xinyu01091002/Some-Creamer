function out = creamer_eta33_h3h4_triad_1d_cpp(eta_lin, x_vec, h, cfg)
%CREAMER_ETA33_H3H4_TRIAD_1D_CPP C++ backend for 1D finite-depth eta33 triads.

if nargin < 4 || isempty(cfg)
    cfg = struct();
end
cfg = local_defaults(cfg);

eta_lin = reshape(real(eta_lin), 1, []);
x_vec = x_vec(:).';
nx = numel(eta_lin);
if numel(x_vec) ~= nx
    error('creamer_eta33_h3h4_triad_1d_cpp:SizeMismatch', ...
        'Length of x_vec (%d) must match eta length (%d).', numel(x_vec), nx);
end
if mod(nx, 2) ~= 0
    error('creamer_eta33_h3h4_triad_1d_cpp:EvenGridRequired', ...
        'This helper assumes an even grid size. Got nx=%d.', nx);
end

dx = mean(diff(x_vec));
Lx = dx * nx;
mx = [0:(nx/2), (-nx/2 + 1):-1];
k = (2*pi/Lx) * mx;
eta_hat = fft(eta_lin) / nx;

if ~exist(cfg.cpp_job_dir, 'dir')
    mkdir(cfg.cpp_job_dir);
end
local_write_int64(fullfile(cfg.cpp_job_dir, 'meta.bin'), int64([nx cfg.max_triad_active_modes]));
local_write_double(fullfile(cfg.cpp_job_dir, 'params.bin'), [dx cfg.g h 0]);
local_write_complex(fullfile(cfg.cpp_job_dir, 'zeta0.bin'), eta_hat(:));

if ~isfile(cfg.cpp_exe)
    error('creamer_eta33_h3h4_triad_1d_cpp:MissingExecutable', ...
        'C++ eta33 triad executable not found: %s', cfg.cpp_exe);
end

cmd = sprintf('"%s" "%s"', cfg.cpp_exe, cfg.cpp_job_dir);
t_cpp = tic;
[status, output] = system(cmd);
cpp_wall_s = toc(t_cpp);
if status ~= 0
    error('creamer_eta33_h3h4_triad_1d_cpp:CppFailed', ...
        'C++ eta33 triad backend failed with status %d:\n%s', status, output);
end

eta33_h3_hat = reshape(local_read_complex(fullfile(cfg.cpp_job_dir, 'eta33_h3_hat.bin'), nx), 1, []);
eta33_h4_delta_hat = reshape(local_read_complex(fullfile(cfg.cpp_job_dir, 'eta33_h4_delta_hat.bin'), nx), 1, []);
eta33_h3h4_hat = eta33_h3_hat + eta33_h4_delta_hat;

summary = local_read_summary(fullfile(cfg.cpp_job_dir, 'summary.txt'));

out = struct();
out.eta33_h3 = real(ifft(eta33_h3_hat * nx));
out.eta33_h4_delta = real(ifft(eta33_h4_delta_hat * nx));
out.eta33_h3h4 = real(ifft(eta33_h3h4_hat * nx));
out.eta33_h3_hat = eta33_h3_hat;
out.eta33_h4_delta_hat = eta33_h4_delta_hat;
out.eta33_h3h4_hat = eta33_h3h4_hat;
out.k = k;
out.mx = mx;
out.diagnostics = summary;
out.diagnostics.cpp_wall_s = cpp_wall_s;
out.diagnostics.backend = 'cpp_eta33_triad_1d';
end

function cfg = local_defaults(cfg)
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
defaults = struct( ...
    'g', 9.81, ...
    'max_triad_active_modes', 32, ...
    'cpp_exe', fullfile(repo_root, 'cpp', 'creamer_flow', 'build', 'creamer_eta33_triad_1d.exe'), ...
    'cpp_job_dir', fullfile(tempdir, 'creamer_eta33_triad_1d_cpp_job'));
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(cfg, names{i}) || isempty(cfg.(names{i}))
        cfg.(names{i}) = defaults.(names{i});
    end
end
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

function summary = local_read_summary(path)
summary = struct();
if ~isfile(path)
    return;
end
lines = splitlines(string(fileread(path)));
for i = 1:numel(lines)
    parts = split(strtrim(lines(i)));
    if numel(parts) ~= 2
        continue;
    end
    summary.(matlab.lang.makeValidName(char(parts(1)))) = str2double(parts(2));
end
end
