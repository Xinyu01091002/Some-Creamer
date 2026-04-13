function eta = linear_focus_gaussian_spectrum_1d(Akp, bandwidth_ratio, x_vec, cfg)
%LINEAR_FOCUS_GAUSSIAN_SPECTRUM_1D Symmetric narrow-band focused 1D spectrum.
%
% bandwidth_ratio is sigma_k / k_p.  The positive-frequency amplitudes are
% normalized so that sum(a_m) is approximately Akp/kp.  The phase is chosen
% so the packet focuses at cfg.focus_x in the physical x_vec coordinates.

if nargin < 4 || isempty(cfg)
    cfg = struct();
end
cfg = local_defaults(cfg);

x_vec = x_vec(:).';
nx = numel(x_vec);
if mod(nx, 2) ~= 0
    error('linear_focus_gaussian_spectrum_1d:EvenGridRequired', ...
        'This helper assumes an even grid size. Got nx=%d.', nx);
end

dx = mean(diff(x_vec));
Lx = dx * nx;
mx = [0:(nx/2), (-nx/2 + 1):-1];
k = (2*pi/Lx) * mx;
pos = 2:(nx/2);
kp = cfg.kp;
A = Akp / kp;
sigma_k = bandwidth_ratio * kp;

shape = exp(-0.5 * ((k(pos) - kp) / sigma_k).^2);
shape(k(pos) <= 0) = 0;
dk = 2*pi / Lx;
amp = A * shape * dk / sum(shape * dk);

eta_hat = complex(zeros(1, nx));
eta_hat(pos) = 0.5 * amp .* exp(1i * (k(pos) * (x_vec(1) - cfg.focus_x) + cfg.phase_shift));
for idx = pos
    neg_idx = mod(nx - idx + 2 - 1, nx) + 1;
    eta_hat(neg_idx) = conj(eta_hat(idx));
end
eta = real(ifft(eta_hat * nx));
end

function cfg = local_defaults(cfg)
defaults = struct( ...
    'kp', 0.0279, ...
    'focus_x', 0, ...
    'phase_shift', 0);
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(cfg, names{i}) || isempty(cfg.(names{i}))
        cfg.(names{i}) = defaults.(names{i});
    end
end
end
