function eta_shifted = apply_global_phase_shift_1d(eta, x_vec, phase_shift)
%APPLY_GLOBAL_PHASE_SHIFT_1D Apply a uniform carrier-phase shift to a real
%1D wave field by rotating the positive Fourier half-plane.

eta = reshape(real(eta), 1, []);
x_vec = x_vec(:).';
nx = numel(eta);

if numel(x_vec) ~= nx
    error('apply_global_phase_shift_1d:SizeMismatch', ...
        'Length of x_vec (%d) must match numel(eta) (%d).', numel(x_vec), nx);
end

dx = mean(diff(x_vec));
Lx = dx * nx;

if mod(nx, 2) ~= 0
    error('apply_global_phase_shift_1d:EvenGridRequired', ...
        'Current helper assumes an even grid size. Got nx=%d.', nx);
end

mx = [0:(nx/2), (-nx/2 + 1):-1];
kx = (2 * pi / Lx) * mx; %#ok<NASGU>

eta_hat = fft(eta) / nx;
shifted_hat = eta_hat;

positive_mask = mx > 0;
negative_mask = mx < 0;

shifted_hat(positive_mask) = eta_hat(positive_mask) .* exp(1i * phase_shift);
shifted_hat(negative_mask) = eta_hat(negative_mask) .* exp(-1i * phase_shift);

eta_shifted = real(ifft(shifted_hat * nx));
end
