function [field_filtered, diagnostics] = directional_frequency_filter(field_in, x_vec, y_vec, cfg)
%DIRECTIONAL_FREQUENCY_FILTER 2D directional band-pass filter in wavenumber space.
%
% The filter is designed for directional wave fields where we want to
% isolate content near a target radial wavenumber and, optionally, near a
% preferred propagation direction.
%
% Required cfg fields:
%   cfg.kp                - peak wavenumber used for normalization
%
% Either:
%   cfg.k_ratio_center    - target radius in units of kp, e.g. 2.0 or 0.5
% and optionally
%   cfg.k_ratio_halfwidth - radial half-width in units of kp
%
% Or:
%   cfg.k_ratio_min       - lower radial cutoff in units of kp
%   cfg.k_ratio_max       - upper radial cutoff in units of kp
%
% Optional cfg fields:
%   cfg.k_ratio_halfwidth - radial half-width in units of kp (default 0.25)
%   cfg.k_ratio_min       - lower radial cutoff in units of kp
%   cfg.k_ratio_max       - upper radial cutoff in units of kp
%   cfg.direction_deg     - preferred direction in degrees (default: [])
%   cfg.angle_halfwidth_deg - angular half-width in degrees (default 180)
%   cfg.mask_type         - 'hard' or 'gaussian' (default 'gaussian')
%   cfg.verbose           - print summary (default true)

if nargin < 4 || isempty(cfg)
    cfg = struct();
end
cfg = local_apply_defaults(cfg);

[ny, nx] = size(field_in);
x_vec = x_vec(:).';
y_vec = y_vec(:).';

if numel(x_vec) ~= nx || numel(y_vec) ~= ny
    error('directional_frequency_filter:SizeMismatch', ...
        'Input field size must match x_vec and y_vec.');
end

dx = mean(diff(x_vec));
dy = mean(diff(y_vec));

[kx, ~] = local_fft_wavenumbers(nx, dx);
[ky, ~] = local_fft_wavenumbers(ny, dy);
[KX, KY] = meshgrid(kx, ky);
K = hypot(KX, KY);
theta = atan2(KY, KX);

if ~isempty(cfg.k_ratio_min) && ~isempty(cfg.k_ratio_max)
    k_min = cfg.k_ratio_min * cfg.kp;
    k_max = cfg.k_ratio_max * cfg.kp;
    if strcmpi(cfg.mask_type, 'hard')
        radial_mask = (K >= k_min) & (K <= k_max);
    else
        k_center = 0.5 * (k_min + k_max);
        k_halfwidth = 0.5 * (k_max - k_min);
        sigma_k = k_halfwidth / sqrt(2 * log(2));
        radial_mask = exp(-0.5 * ((K - k_center) ./ max(sigma_k, eps)).^2);
        radial_mask((K < k_min) | (K > k_max)) = 0;
    end
else
    k_center = cfg.k_ratio_center * cfg.kp;
    k_halfwidth = cfg.k_ratio_halfwidth * cfg.kp;
    if strcmpi(cfg.mask_type, 'hard')
        radial_mask = abs(K - k_center) <= k_halfwidth;
    else
        sigma_k = k_halfwidth / sqrt(2 * log(2));
        radial_mask = exp(-0.5 * ((K - k_center) ./ max(sigma_k, eps)).^2);
    end
    k_min = k_center - k_halfwidth;
    k_max = k_center + k_halfwidth;
end

if isempty(cfg.direction_deg)
    angular_mask = ones(size(K));
else
    theta0 = deg2rad(cfg.direction_deg);
    dtheta = angle(exp(1i * (theta - theta0)));
    dtheta_mirror = angle(exp(1i * (theta - theta0 - pi)));
    dtheta_eff = min(abs(dtheta), abs(dtheta_mirror));
    if strcmpi(cfg.mask_type, 'hard')
        angular_mask = double(dtheta_eff <= deg2rad(cfg.angle_halfwidth_deg));
    else
        sigma_theta = deg2rad(cfg.angle_halfwidth_deg) / sqrt(2 * log(2));
        angular_mask = exp(-0.5 * (dtheta_eff ./ max(sigma_theta, eps)).^2);
    end
end

mask = radial_mask .* angular_mask;
spec_in = fft2(field_in);
spec_filtered = spec_in .* mask;
field_filtered = real(ifft2(spec_filtered));

diagnostics = struct();
diagnostics.kx = kx;
diagnostics.ky = ky;
diagnostics.mask = mask;
diagnostics.k_min = k_min;
diagnostics.k_max = k_max;
diagnostics.direction_deg = cfg.direction_deg;
diagnostics.angle_halfwidth_deg = cfg.angle_halfwidth_deg;
diagnostics.mask_type = cfg.mask_type;
diagnostics.input_max = max(abs(field_in(:)));
diagnostics.output_max = max(abs(field_filtered(:)));

if cfg.verbose
    fprintf('Directional frequency filter:\n');
    fprintf('  radial band         = [%.6g, %.6g] kp\n', k_min / cfg.kp, k_max / cfg.kp);
    if isempty(cfg.direction_deg)
        fprintf('  direction           = all directions\n');
    else
        fprintf('  direction           = %.6g deg +- %.6g deg\n', ...
            cfg.direction_deg, cfg.angle_halfwidth_deg);
    end
    fprintf('  mask type           = %s\n', cfg.mask_type);
    fprintf('  max |input|         = %.6g\n', diagnostics.input_max);
    fprintf('  max |filtered|      = %.6g\n', diagnostics.output_max);
end
end

function cfg = local_apply_defaults(cfg)
defaults = struct( ...
    'kp', [], ...
    'k_ratio_center', [], ...
    'k_ratio_halfwidth', 0.25, ...
    'k_ratio_min', [], ...
    'k_ratio_max', [], ...
    'direction_deg', [], ...
    'angle_halfwidth_deg', 180, ...
    'mask_type', 'gaussian', ...
    'verbose', true);

names = fieldnames(defaults);
for n = 1:numel(names)
    name = names{n};
    if ~isfield(cfg, name) || isempty(cfg.(name))
        cfg.(name) = defaults.(name);
    end
end

if isempty(cfg.kp)
    error('directional_frequency_filter:MissingConfig', ...
        'cfg.kp is required.');
end

has_band = ~isempty(cfg.k_ratio_min) && ~isempty(cfg.k_ratio_max);
has_center = ~isempty(cfg.k_ratio_center);

if ~(has_band || has_center)
    error('directional_frequency_filter:MissingConfig', ...
        'Provide either (k_ratio_min,k_ratio_max) or k_ratio_center.');
end
end

function [kvec, mvec] = local_fft_wavenumbers(n, d)
L = n * d;
mvec = [0:(n/2), (-n/2 + 1):-1];
kvec = (2 * pi / L) * mvec;
end
