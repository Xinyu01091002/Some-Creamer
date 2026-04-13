function out = creamer_eta33_h3h4_single_frequency_1d(eta_lin, x_vec, h)
%CREAMER_ETA33_H3H4_SINGLE_FREQUENCY_1D Self-interaction eta33 check.
%
% Builds the one-frequency third superharmonic contribution implied by the
% finite-depth single-frequency limit:
%
%   eta_33(k -> 3k) = 4 C33(kh) k^2 eta_hat(k)^3.
%
% The H3-only coefficient is the finite-depth Creamer lambda-flow result.
% The H3+H4 coefficient is the Stokes/VWA coefficient recovered after
% absorbing the non-resonant H4 part in the Mathematica check.
%
% This is exact for a regular single-frequency wave. For a broadband focused
% packet it is a self-interaction diagnostic only; it does not include all
% cross-triad third-order MF12 interactions.

eta_lin = reshape(real(eta_lin), 1, []);
x_vec = x_vec(:).';
nx = numel(eta_lin);

if numel(x_vec) ~= nx
    error('creamer_eta33_h3h4_single_frequency_1d:SizeMismatch', ...
        'Length of x_vec (%d) must match eta length (%d).', numel(x_vec), nx);
end
if mod(nx, 2) ~= 0
    error('creamer_eta33_h3h4_single_frequency_1d:EvenGridRequired', ...
        'This helper assumes an even grid size. Got nx=%d.', nx);
end

dx = mean(diff(x_vec));
Lx = dx * nx;
mx = [0:(nx/2), (-nx/2 + 1):-1];
k = (2*pi/Lx) * mx;
eta_hat = fft(eta_lin) / nx;

eta33_h3_hat = complex(zeros(1, nx));
eta33_h3h4_hat = complex(zeros(1, nx));

positive_idx = find(mx > 0);
for ii = 1:numel(positive_idx)
    idx = positive_idx(ii);
    dest_mode = 3 * mx(idx);
    if dest_mode > nx/2
        continue;
    end

    dest_idx = find(mx == dest_mode, 1);
    if isempty(dest_idx)
        continue;
    end

    kk = abs(k(idx));
    if kk == 0
        continue;
    end

    sigma = tanh(kk * h);
    c33_h3 = local_creamer_h3_c33(sigma);
    c33_h3h4 = local_stokes_c33(sigma);

    z3 = 4 * kk^2 * eta_hat(idx)^3;
    eta33_h3_hat(dest_idx) = eta33_h3_hat(dest_idx) + c33_h3 * z3;
    eta33_h3h4_hat(dest_idx) = eta33_h3h4_hat(dest_idx) + c33_h3h4 * z3;
end

eta33_h3_hat = local_enforce_hermitian(eta33_h3_hat, mx);
eta33_h3h4_hat = local_enforce_hermitian(eta33_h3h4_hat, mx);

out = struct();
out.eta33_h3 = real(ifft(eta33_h3_hat * nx));
out.eta33_h3h4 = real(ifft(eta33_h3h4_hat * nx));
out.eta33_h3_hat = eta33_h3_hat;
out.eta33_h3h4_hat = eta33_h3h4_hat;
out.k = k;
out.mx = mx;
end

function c = local_creamer_h3_c33(sigma)
c = (3 * (1 + sigma.^2) .* ...
    (54 + 39*sigma.^2 + 70*sigma.^4 - 28*sigma.^6 - 4*sigma.^8 - 3*sigma.^10)) ./ ...
    (64 * sigma.^6 .* (9 + 14*sigma.^2 + 9*sigma.^4));
end

function c = local_stokes_c33(sigma)
c = (27 - 9*sigma.^2 + 9*sigma.^4 - 3*sigma.^6) ./ (64 * sigma.^6);
end

function spec = local_enforce_hermitian(spec, mx)
for idx = 1:numel(mx)
    m = mx(idx);
    if m <= 0
        continue;
    end
    neg_idx = find(mx == -m, 1);
    if ~isempty(neg_idx)
        spec(neg_idx) = conj(spec(idx));
    end
end
zero_idx = find(mx == 0, 1);
if ~isempty(zero_idx)
    spec(zero_idx) = real(spec(zero_idx));
end
nyquist_idx = find(mx == max(mx), 1);
if ~isempty(nyquist_idx)
    spec(nyquist_idx) = real(spec(nyquist_idx));
end
end

