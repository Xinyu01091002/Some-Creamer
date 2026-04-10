function eta_shifted = apply_global_phase_shift_2d(eta, x_vec, y_vec, phase_shift)
%APPLY_GLOBAL_PHASE_SHIFT_2D Apply a uniform carrier-phase shift to a real
%2D wave field by rotating one representative half-plane of Fourier modes.

[ny, nx] = size(eta);
x_vec = x_vec(:).';
y_vec = y_vec(:).';
dx = mean(diff(x_vec));
dy = mean(diff(y_vec));
Lx = dx * nx;
Ly = dy * ny;
n_total = nx * ny;

mx = [0:(nx/2), (-nx/2 + 1):-1];
my = [0:(ny/2), (-ny/2 + 1):-1];
kx = (2 * pi / Lx) * mx;
ky = (2 * pi / Ly) * my;
[KX, KY] = meshgrid(kx, ky);

eta_hat = fft2(eta) / n_total;
shifted_hat = eta_hat;

positive_mask = (KX > 0) | (KX == 0 & KY > 0);
negative_mask = (KX < 0) | (KX == 0 & KY < 0);

shifted_hat(positive_mask) = eta_hat(positive_mask) .* exp(1i * phase_shift);
shifted_hat(negative_mask) = eta_hat(negative_mask) .* exp(-1i * phase_shift);

eta_shifted = real(ifft2(shifted_hat * n_total));
end
