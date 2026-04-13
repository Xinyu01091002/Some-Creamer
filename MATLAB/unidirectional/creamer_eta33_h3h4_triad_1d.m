function out = creamer_eta33_h3h4_triad_1d(eta_lin, x_vec, h, cfg)
%CREAMER_ETA33_H3H4_TRIAD_1D Prototype 1D broadband H3+H4 eta33.
%
% This correctness-first prototype keeps a small set of positive active
% parent modes, builds the quartic polynomial K4 from direct H4 convolution
% and -1/2{H3,W3}, absorbs only non-resonant normal-variable monomials, and
% evaluates {zeta,W4} on the linear wave state.
%
% It is intended for 1D finite-depth diagnostics, not production-size 2D runs.
% Current status: experimental broadband prototype. The H3 and H3+H4 pieces
% are calibrated against the single-frequency finite-depth Stokes/MF12 limit,
% but the sparse-dictionary implementation is still slow for broad active sets.

if nargin < 4 || isempty(cfg)
    cfg = struct();
end
cfg = local_defaults(cfg);

eta_lin = reshape(real(eta_lin), 1, []);
x_vec = x_vec(:).';
nx = numel(eta_lin);
if numel(x_vec) ~= nx
    error('creamer_eta33_h3h4_triad_1d:SizeMismatch', ...
        'Length of x_vec (%d) must match eta length (%d).', numel(x_vec), nx);
end
if mod(nx, 2) ~= 0
    error('creamer_eta33_h3h4_triad_1d:EvenGridRequired', ...
        'This helper assumes an even grid size. Got nx=%d.', nx);
end

dx = mean(diff(x_vec));
Lx = dx * nx;
mx_grid = [0:(nx/2), (-nx/2 + 1):-1];
k_grid = (2*pi/Lx) * mx_grid;
eta_hat = fft(eta_lin) / nx;
theta_grid = abs(k_grid) .* tanh(abs(k_grid) * h);
theta_grid(mx_grid == 0) = 0;

positive_idx = find(mx_grid > 0);
[~, order] = sort(abs(eta_hat(positive_idx)).^2, 'descend');
keep_count = min(cfg.max_triad_active_modes, numel(order));
positive_active_idx = positive_idx(order(1:keep_count));
positive_modes = sort(mx_grid(positive_active_idx));
linear_modes = unique([positive_modes, -positive_modes]);

max_abs_linear_mode = max(abs(linear_modes));
mode_set = unique([linear_modes, ...
    local_pair_sums(linear_modes), local_triple_sums(linear_modes)]);
mode_set = mode_set(mode_set ~= 0);
mode_set = mode_set(abs(mode_set) <= nx/2);
mode_set = sort(mode_set(:).');

zlin = containers.Map('KeyType', 'double', 'ValueType', 'any');
plin = containers.Map('KeyType', 'double', 'ValueType', 'any');
for m = mode_set
    idx = local_mode_to_index(m, nx);
    zlin(m) = eta_hat(idx);
    if m == 0 || theta_grid(idx) == 0
        plin(m) = 0;
    elseif m > 0
        plin(m) = -1i * sqrt(cfg.g / theta_grid(idx)) * eta_hat(idx);
    else
        plin(m) = 1i * sqrt(cfg.g / theta_grid(idx)) * eta_hat(idx);
    end
end

[z2, p2] = local_v2(zlin, plin, mode_set, nx, k_grid, theta_grid, cfg.g);
z3_h3 = local_v2_cross_z(zlin, plin, z2, p2, mode_set, nx, k_grid, theta_grid, cfg.g);
z3_h3 = local_scale_map(z3_h3, 0.5);

target_outputs = local_third_superharmonic_targets(positive_modes, nx);
if strcmpi(cfg.evaluation_model, 'dictionary')
    H3 = local_build_h3_poly(mode_set, nx, k_grid, theta_grid);
    W3 = local_build_w3_poly(mode_set, nx, k_grid, theta_grid, cfg.g);
    H4 = local_build_h4_direct_poly(mode_set, nx, k_grid, theta_grid);
    br = local_poisson_bracket(H3, W3, mode_set);
    K4 = local_poly_add(local_poly_scale(H4, cfg.h4_scale), br, cfg.bracket_scale);

    [z3_h4, resonance_stats] = local_eval_z_w4(K4, zlin, plin, target_outputs, nx, k_grid, theta_grid, cfg);
    target_normal_coeffs = local_target_normal_coeffs(K4, nx, theta_grid, cfg);
    poly_counts = struct('n_h3_terms', local_poly_count(H3), ...
        'n_w3_terms', local_poly_count(W3), ...
        'n_h4_terms', local_poly_count(H4), ...
        'n_k4_terms', local_poly_count(K4));
else
    H3 = local_build_h3_terms(mode_set, nx, k_grid, theta_grid);
    W3 = local_build_w3_terms(mode_set, nx, k_grid, theta_grid, cfg.g);
    support_modes = local_target_support_modes(cfg, linear_modes, mode_set);
    [z3_h4, resonance_stats, poly_counts] = local_eval_z_w4_targeted( ...
        H3, W3, zlin, plin, mode_set, linear_modes, support_modes, target_outputs, nx, k_grid, theta_grid, cfg);
    target_normal_coeffs = struct('available', false);
end

eta33_h3_hat = local_map_to_spec(z3_h3, nx);
eta33_h4_delta_hat = local_map_to_spec(z3_h4, nx);
eta33_h3h4_hat = eta33_h3_hat + eta33_h4_delta_hat;

eta33_h3_hat = local_keep_third_superharmonic(eta33_h3_hat, positive_modes, nx);
eta33_h4_delta_hat = local_keep_third_superharmonic(eta33_h4_delta_hat, positive_modes, nx);
eta33_h3h4_hat = local_keep_third_superharmonic(eta33_h3h4_hat, positive_modes, nx);

out = struct();
out.eta33_h3 = real(ifft(eta33_h3_hat * nx));
out.eta33_h4_delta = real(ifft(eta33_h4_delta_hat * nx));
out.eta33_h3h4 = real(ifft(eta33_h3h4_hat * nx));
out.eta33_h3_hat = eta33_h3_hat;
out.eta33_h4_delta_hat = eta33_h4_delta_hat;
out.eta33_h3h4_hat = eta33_h3h4_hat;
out.k = k_grid;
out.mx = mx_grid;
out.diagnostics = struct( ...
    'positive_modes', positive_modes, ...
    'linear_modes', linear_modes, ...
    'mode_set', mode_set, ...
    'target_outputs', target_outputs, ...
    'max_abs_linear_mode', max_abs_linear_mode, ...
    'n_h3_terms', poly_counts.n_h3_terms, ...
    'n_w3_terms', poly_counts.n_w3_terms, ...
    'n_h4_terms', poly_counts.n_h4_terms, ...
    'n_k4_terms', poly_counts.n_k4_terms, ...
    'resonant_terms', resonance_stats.resonant_terms, ...
    'nonresonant_terms', resonance_stats.nonresonant_terms, ...
    'candidate_terms', resonance_stats.candidate_terms, ...
    'used_terms', resonance_stats.used_terms, ...
    'evaluation_model', cfg.evaluation_model, ...
    'support_modes_kind', cfg.support_modes_kind, ...
    'target_normal_coeffs', target_normal_coeffs, ...
    'pruning_stats', resonance_stats.pruning_stats, ...
    'status', 'experimental_broadband_single_frequency_calibrated');
end

function cfg = local_defaults(cfg)
defaults = struct( ...
    'g', 9.81, ...
    'max_triad_active_modes', 3, ...
    'resonance_tol', 1e-10, ...
    'homological_inverse_sign', -1, ...
    'h4_scale', 1, ...
    'bracket_scale', -0.5, ...
    'evaluation_model', 'targeted', ...
    'support_modes_kind', '');
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(cfg, names{i}) || isempty(cfg.(names{i}))
        cfg.(names{i}) = defaults.(names{i});
    end
end
if isempty(cfg.support_modes_kind)
    if strcmpi(cfg.evaluation_model, 'less_pruned') || local_is_semi_pruned(cfg)
        cfg.support_modes_kind = 'mode_set';
    else
        cfg.support_modes_kind = 'linear_modes';
    end
end
end

function support_modes = local_target_support_modes(cfg, linear_modes, mode_set)
if strcmpi(cfg.support_modes_kind, 'mode_set')
    support_modes = mode_set;
else
    support_modes = linear_modes;
end
support_modes = sort(unique(support_modes(:).'));
end

function tf = local_is_semi_pruned(cfg)
tf = strcmpi(cfg.evaluation_model, 'semi_pruned') || strcmpi(cfg.evaluation_model, 'semi_pruned_v2');
end

function tf = local_is_semi_pruned_v2(cfg)
tf = strcmpi(cfg.evaluation_model, 'semi_pruned_v2');
end

function sums = local_pair_sums(modes)
[a, b] = ndgrid(modes, modes);
sums = unique(a(:) + b(:)).';
end

function sums = local_triple_sums(modes)
[a, b, c] = ndgrid(modes, modes, modes);
sums = unique(a(:) + b(:) + c(:)).';
end

function idx = local_mode_to_index(m, nx)
idx = mod(m, nx) + 1;
end

function val = local_z(map_obj, m)
if isKey(map_obj, m)
    val = map_obj(m);
else
    val = 0;
end
end

function [zout, pout] = local_v2(zmap, pmap, mode_set, nx, k_grid, theta_grid, g)
zout = containers.Map('KeyType', 'double', 'ValueType', 'any');
pout = containers.Map('KeyType', 'double', 'ValueType', 'any');
for dest = mode_set
    zsum = 0;
    psum = 0;
    for a = mode_set
        b = dest - a;
        if ~isKey(zmap, b) && ~isKey(pmap, b)
            continue;
        end
        [D, Bz, Bphi] = local_b_d(dest, a, b, nx, k_grid, theta_grid, g);
        za = local_z(zmap, a); zb = local_z(zmap, b);
        pa = local_z(pmap, a); pb = local_z(pmap, b);
        zsum = zsum + 3 * D * pa * pb + Bz * za * zb;
        psum = psum - 2 * Bphi * za * pb;
    end
    if abs(zsum) > 0
        zout(dest) = zsum;
    end
    if abs(psum) > 0
        pout(dest) = psum;
    end
end
end

function zout = local_v2_cross_z(z1, p1, z2, p2, mode_set, nx, k_grid, theta_grid, g)
zout = containers.Map('KeyType', 'double', 'ValueType', 'any');
for dest = mode_set
    zsum = 0;
    for a = mode_set
        b = dest - a;
        [D, Bz, ~] = local_b_d(dest, a, b, nx, k_grid, theta_grid, g);
        zsum = zsum + 3 * D * (local_z(p2, a) * local_z(p1, b) + local_z(p1, a) * local_z(p2, b)) ...
            + Bz * (local_z(z2, a) * local_z(z1, b) + local_z(z1, a) * local_z(z2, b));
    end
    if abs(zsum) > 0
        zout(dest) = zsum;
    end
end
end

function out = local_scale_map(in, scale)
out = containers.Map('KeyType', 'double', 'ValueType', 'any');
keys_in = cell2mat(keys(in));
for m = keys_in
    out(m) = scale * in(m);
end
end

function spec = local_map_to_spec(map_obj, nx)
spec = complex(zeros(1, nx));
keys_in = cell2mat(keys(map_obj));
for m = keys_in
    spec(local_mode_to_index(m, nx)) = spec(local_mode_to_index(m, nx)) + map_obj(m);
end
spec = local_enforce_hermitian(spec, nx);
end

function spec = local_keep_third_superharmonic(spec, positive_modes, nx)
keep = false(1, nx);
targets = local_third_superharmonic_targets(positive_modes, nx);
for m = targets(:).'
    keep(local_mode_to_index(m, nx)) = true;
    keep(local_mode_to_index(-m, nx)) = true;
end
spec(~keep) = 0;
spec = local_enforce_hermitian(spec, nx);
end

function targets = local_third_superharmonic_targets(positive_modes, nx)
[a, b, c] = ndgrid(positive_modes, positive_modes, positive_modes);
targets = unique(a(:) + b(:) + c(:));
targets = targets(abs(targets) <= nx/2);
targets = sort(targets(:).');
end

function spec = local_enforce_hermitian(spec, nx)
mx = [0:(nx/2), (-nx/2 + 1):-1];
for idx = 1:numel(mx)
    m = mx(idx);
    if m <= 0
        continue;
    end
    neg_idx = local_mode_to_index(-m, nx);
    spec(neg_idx) = conj(spec(idx));
end
spec(local_mode_to_index(0, nx)) = real(spec(local_mode_to_index(0, nx)));
spec(local_mode_to_index(nx/2, nx)) = real(spec(local_mode_to_index(nx/2, nx)));
end

function H3 = local_build_h3_poly(mode_set, nx, k_grid, theta_grid)
H3 = containers.Map('KeyType', 'char', 'ValueType', 'any');
for a = mode_set
    for b = mode_set
        c = -(a + b);
        if ~ismember(c, mode_set)
            continue;
        end
        thb = local_theta(b, nx, theta_grid);
        thc = local_theta(c, nx, theta_grid);
        kb = local_k(b, nx, k_grid);
        kc = local_k(c, nx, k_grid);
        ka = local_k(a, nx, k_grid);
        tha = local_theta(a, nx, theta_grid);
        coeff = 0.25 * (ka^2 + kb^2 - kc^2 - 2 * tha * thb);
        H3 = local_add_poly_term(H3, {local_var('z', c), local_var('p', a), local_var('p', b)}, coeff);
    end
end
end

function W3 = local_build_w3_poly(mode_set, nx, k_grid, theta_grid, g)
W3 = containers.Map('KeyType', 'char', 'ValueType', 'any');
for a = mode_set
    for b = mode_set
        c = -(a + b);
        if ~ismember(c, mode_set)
            continue;
        end
        [D, ~, Bphi] = local_b_d(a, b, c, nx, k_grid, theta_grid, g);
        W3 = local_add_poly_term(W3, {local_var('p', a), local_var('p', b), local_var('p', c)}, D);
        W3 = local_add_poly_term(W3, {local_var('z', a), local_var('z', b), local_var('p', c)}, Bphi);
    end
end
end

function H3 = local_build_h3_terms(mode_set, nx, k_grid, theta_grid)
terms = repmat(struct('vars', {{}}, 'coeff', 0), 1, numel(mode_set)^2);
n = 0;
for a = mode_set
    for b = mode_set
        c = -(a + b);
        if ~ismember(c, mode_set)
            continue;
        end
        thb = local_theta(b, nx, theta_grid);
        kb = local_k(b, nx, k_grid);
        ka = local_k(a, nx, k_grid);
        tha = local_theta(a, nx, theta_grid);
        kc = local_k(c, nx, k_grid);
        coeff = 0.25 * (ka^2 + kb^2 - kc^2 - 2 * tha * thb);
        if coeff == 0
            continue;
        end
        n = n + 1;
        terms(n).vars = {local_var('z', c), local_var('p', a), local_var('p', b)};
        terms(n).coeff = coeff;
    end
end
H3 = terms(1:n);
end

function W3 = local_build_w3_terms(mode_set, nx, k_grid, theta_grid, g)
terms = repmat(struct('vars', {{}}, 'coeff', 0), 1, 2 * numel(mode_set)^2);
n = 0;
for a = mode_set
    for b = mode_set
        c = -(a + b);
        if ~ismember(c, mode_set)
            continue;
        end
        [D, ~, Bphi] = local_b_d(a, b, c, nx, k_grid, theta_grid, g);
        if D ~= 0
            n = n + 1;
            terms(n).vars = {local_var('p', a), local_var('p', b), local_var('p', c)};
            terms(n).coeff = D;
        end
        if Bphi ~= 0
            n = n + 1;
            terms(n).vars = {local_var('z', a), local_var('z', b), local_var('p', c)};
            terms(n).coeff = Bphi;
        end
    end
end
W3 = terms(1:n);
end

function H4 = local_build_h4_direct_poly(mode_set, nx, k_grid, theta_grid)
H4 = containers.Map('KeyType', 'char', 'ValueType', 'any');
for k = mode_set
    for b = mode_set
        for d = mode_set
            c = k - b - d;
            if ~ismember(c, mode_set)
                continue;
            end
            coeff = 0.5 * (local_theta(k, nx, theta_grid) * local_theta(k - b, nx, theta_grid) * local_theta(d, nx, theta_grid) ...
                - 0.5 * local_theta(k, nx, theta_grid) * local_k(d, nx, k_grid)^2 ...
                - 0.5 * local_k(k, nx, k_grid)^2 * local_theta(d, nx, theta_grid));
            H4 = local_add_poly_term(H4, {local_var('p', -k), local_var('z', b), local_var('z', c), local_var('p', d)}, coeff);
        end
    end
end
end

function [z3_h4, stats] = local_eval_z_w4(K4, zlin, plin, target_outputs, nx, k_grid, theta_grid, cfg)
z3_h4 = containers.Map('KeyType', 'double', 'ValueType', 'any');
terms = local_poly_terms(K4);
linear_modes = sort(cell2mat(keys(zlin)));
stats = local_empty_w4_stats();
for i = 1:numel(terms)
    [z3_h4, stats] = local_accumulate_w4_term(z3_h4, stats, terms(i).vars, ...
        terms(i).coeff, zlin, plin, linear_modes, target_outputs, nx, theta_grid, cfg);
end
end

function [z3_h4, stats, counts] = local_eval_z_w4_targeted(H3, W3, zlin, plin, mode_set, linear_modes, support_modes, target_outputs, nx, k_grid, theta_grid, cfg)
z3_h4 = containers.Map('KeyType', 'double', 'ValueType', 'any');
stats = local_empty_w4_stats();

n_h4_terms = 0;
[h4_mode_rows, h4_coeffs, h4_diag] = local_targeted_h4_terms(linear_modes, support_modes, target_outputs, nx, k_grid, theta_grid, cfg);
stats.pruning_stats.targeted_h4_combo_total = h4_diag.combo_total;
stats.pruning_stats.targeted_h4_consistency_pass = h4_diag.consistency_pass;
stats.pruning_stats.targeted_h4_nonzero_pass = h4_diag.nonzero_pass;
stats.pruning_stats.targeted_h4_unique_rows = size(h4_mode_rows, 1);
stats.pruning_stats.targeted_h4_support_pruned = h4_diag.support_pruned;
for row = 1:size(h4_mode_rows, 1)
    n_h4_terms = n_h4_terms + 1;
    [z3_h4, stats] = local_accumulate_w4_term(z3_h4, stats, ...
        {local_var('p', h4_mode_rows(row, 1)), local_var('z', h4_mode_rows(row, 2)), ...
         local_var('z', h4_mode_rows(row, 3)), local_var('p', h4_mode_rows(row, 4))}, ...
        h4_coeffs(row), zlin, plin, support_modes, target_outputs, nx, theta_grid, cfg);
end

n_bracket_terms = 0;
for k = mode_set
    f_z = local_terms_derivative(H3, local_var('z', k));
    if ~isempty(f_z)
        g_p = local_terms_derivative(W3, local_var('p', -k));
        [z3_h4, stats, n_added] = local_accumulate_w4_product(z3_h4, stats, f_z, g_p, ...
            cfg.bracket_scale, zlin, plin, support_modes, target_outputs, nx, theta_grid, cfg);
        n_bracket_terms = n_bracket_terms + n_added;
    end
    f_p = local_terms_derivative(H3, local_var('p', k));
    if ~isempty(f_p)
        g_z = local_terms_derivative(W3, local_var('z', -k));
        [z3_h4, stats, n_added] = local_accumulate_w4_product(z3_h4, stats, f_p, g_z, ...
            -cfg.bracket_scale, zlin, plin, support_modes, target_outputs, nx, theta_grid, cfg);
        n_bracket_terms = n_bracket_terms + n_added;
    end
end

counts = struct('n_h3_terms', numel(H3), ...
    'n_w3_terms', numel(W3), ...
    'n_h4_terms', n_h4_terms, ...
    'n_k4_terms', n_h4_terms + n_bracket_terms);
end

function [mode_rows, coeffs, diag] = local_targeted_h4_terms(linear_modes, support_modes, target_outputs, nx, k_grid, theta_grid, cfg)
rows = zeros(0, 4);
coeffs = zeros(0, 1);
diag = struct('combo_total', 0, 'consistency_pass', 0, 'nonzero_pass', 0, 'support_pruned', 0);
for target = target_outputs
    for forced_pos = 1:4
        for a = support_modes
            for b = support_modes
                for c = support_modes
                    diag.combo_total = diag.combo_total + 1;
                    modes = [a, b, c];
                    if local_h4_combo_support_pruned(modes, linear_modes, cfg)
                        diag.support_pruned = diag.support_pruned + 1;
                        continue;
                    end
                    term_modes = zeros(1, 4);
                    term_modes(forced_pos) = -target;
                    term_modes(setdiff(1:4, forced_pos, 'stable')) = modes;
                    k = -term_modes(1);
                    zb = term_modes(2);
                    zc = term_modes(3);
                    pd = term_modes(4);
                    if zc ~= k - zb - pd
                        continue;
                    end
                    diag.consistency_pass = diag.consistency_pass + 1;
                    coeff = cfg.h4_scale * 0.5 * (local_theta(k, nx, theta_grid) * local_theta(k - zb, nx, theta_grid) * local_theta(pd, nx, theta_grid) ...
                        - 0.5 * local_theta(k, nx, theta_grid) * local_k(pd, nx, k_grid)^2 ...
                        - 0.5 * local_k(k, nx, k_grid)^2 * local_theta(pd, nx, theta_grid));
                    if coeff == 0
                        continue;
                    end
                    diag.nonzero_pass = diag.nonzero_pass + 1;
                    rows(end+1, :) = term_modes; %#ok<AGROW>
                    coeffs(end+1, 1) = coeff; %#ok<AGROW>
                end
            end
        end
    end
end
[mode_rows, ~, idx] = unique(rows, 'rows', 'stable');
merged = accumarray(idx, coeffs, [], @sum);
keep = (merged ~= 0);
mode_rows = mode_rows(keep, :);
coeffs = merged(keep);
end

function tf = local_h4_combo_support_pruned(modes, linear_modes, cfg)
if strcmpi(cfg.support_modes_kind, 'linear_modes')
    tf = any(~ismember(modes, linear_modes));
    return;
end
if local_is_semi_pruned(cfg)
    tf = sum(~ismember(modes, linear_modes)) > 1;
    return;
end
tf = false;
end

function stats = local_empty_w4_stats()
stats = struct('resonant_terms', 0, 'nonresonant_terms', 0, ...
    'candidate_terms', 0, 'used_terms', 0, ...
    'pruning_stats', local_empty_pruning_stats());
end

function stats = local_empty_pruning_stats()
stats = struct( ...
    'terms_seen', 0, ...
    'terms_rejected_target_gate', 0, ...
    'derivative_slots_seen', 0, ...
    'derivative_slots_rejected_nonlinear_support', 0, ...
    'normal_forms_seen', 0, ...
    'normal_forms_zero_coeff', 0, ...
    'normal_forms_resonant', 0, ...
    'normal_forms_used', 0, ...
    'targeted_h4_combo_total', 0, ...
    'targeted_h4_consistency_pass', 0, ...
    'targeted_h4_nonzero_pass', 0, ...
    'targeted_h4_unique_rows', 0, ...
    'targeted_h4_support_pruned', 0, ...
    'targeted_bracket_terms_a_seen', 0, ...
    'targeted_bracket_needed_keys', 0, ...
    'targeted_bracket_matched_keys', 0, ...
    'targeted_bracket_products_formed', 0, ...
    'targeted_bracket_support_pruned', 0);
end

function [z3_h4, stats, n_added] = local_accumulate_w4_product(z3_h4, stats, A, B, scale, zlin, plin, linear_modes, target_outputs, nx, theta_grid, cfg)
n_added = 0;
if isempty(A) || isempty(B) || scale == 0
    return;
end
terms_a = A;
terms_b = B;
b_by_modes = local_terms_by_mode_pair(terms_b);
for i = 1:numel(terms_a)
    stats.pruning_stats.targeted_bracket_terms_a_seen = stats.pruning_stats.targeted_bracket_terms_a_seen + 1;
    modes_a = local_term_modes(terms_a(i).vars);
    [needed_keys, n_support_pruned] = local_needed_pair_keys(modes_a, linear_modes, target_outputs, cfg);
    stats.pruning_stats.targeted_bracket_needed_keys = stats.pruning_stats.targeted_bracket_needed_keys + numel(needed_keys);
    stats.pruning_stats.targeted_bracket_support_pruned = stats.pruning_stats.targeted_bracket_support_pruned + n_support_pruned;
    for key_idx = 1:numel(needed_keys)
        if ~isKey(b_by_modes, needed_keys{key_idx})
            continue;
        end
        stats.pruning_stats.targeted_bracket_matched_keys = stats.pruning_stats.targeted_bracket_matched_keys + 1;
        b_indices = b_by_modes(needed_keys{key_idx});
        for b_pos = 1:numel(b_indices)
            j = b_indices(b_pos);
            n_added = n_added + 1;
            stats.pruning_stats.targeted_bracket_products_formed = stats.pruning_stats.targeted_bracket_products_formed + 1;
            [z3_h4, stats] = local_accumulate_w4_term(z3_h4, stats, ...
                [terms_a(i).vars, terms_b(j).vars], scale * terms_a(i).coeff * terms_b(j).coeff, ...
                zlin, plin, linear_modes, target_outputs, nx, theta_grid, cfg);
        end
    end
end
end

function out = local_terms_derivative(terms, target)
out = repmat(struct('vars', {{}}, 'coeff', 0), 1, numel(terms));
n = 0;
for i = 1:numel(terms)
    [new_vars, count] = local_derivative_term(terms(i).vars, target);
    if count ~= 0
        n = n + 1;
        out(n).vars = new_vars;
        out(n).coeff = terms(i).coeff * count;
    end
end
out = out(1:n);
end

function by_modes = local_terms_by_mode_pair(terms)
by_modes = containers.Map('KeyType', 'char', 'ValueType', 'any');
for i = 1:numel(terms)
    key = local_mode_pair_key(local_term_modes(terms(i).vars));
    if isKey(by_modes, key)
        by_modes(key) = [by_modes(key), i];
    else
        by_modes(key) = i;
    end
end
end

function [keys_needed, n_support_pruned] = local_needed_pair_keys(modes_a, linear_modes, target_outputs, cfg)
keys_needed = {};
n_support_pruned = 0;
for target = target_outputs
    target_mode = -target;
    target_count = sum(modes_a == target_mode);
    if target_count > 1
        continue;
    end
    if target_count == 1
        rest = modes_a(modes_a ~= target_mode);
        if local_rest_support_allowed(rest, linear_modes, cfg)
            partner_pool = local_partner_mode_pool(linear_modes, target_mode, cfg);
            for i = 1:numel(partner_pool)
                for j = 1:numel(partner_pool)
                    pair = [partner_pool(i), partner_pool(j)];
                    if local_pair_support_allowed(pair, linear_modes, target_mode, cfg)
                        keys_needed{end+1} = local_mode_pair_key(pair); %#ok<AGROW>
                    else
                        n_support_pruned = n_support_pruned + 1;
                    end
                end
            end
        end
    else
        if local_rest_support_allowed(modes_a, linear_modes, cfg)
            partner_pool = local_partner_mode_pool(linear_modes, target_mode, cfg);
            for i = 1:numel(partner_pool)
                pair = [target_mode, partner_pool(i)];
                if local_pair_support_allowed(pair, linear_modes, target_mode, cfg)
                    keys_needed{end+1} = local_mode_pair_key(pair); %#ok<AGROW>
                else
                    n_support_pruned = n_support_pruned + 1;
                end
            end
        end
    end
end
if ~isempty(keys_needed)
    keys_needed = unique(keys_needed);
end
end

function tf = local_rest_support_allowed(modes, linear_modes, cfg)
if strcmpi(cfg.support_modes_kind, 'linear_modes')
    tf = all(ismember(modes, linear_modes));
    return;
end
if local_is_semi_pruned(cfg)
    tf = sum(~ismember(modes, linear_modes)) <= 1;
    return;
end
tf = true;
end

function pool = local_partner_mode_pool(linear_modes, target_mode, cfg)
if strcmpi(cfg.support_modes_kind, 'linear_modes')
    pool = linear_modes;
else
    pool = unique([linear_modes, target_mode]);
end
pool = sort(unique(pool(:).'));
end

function tf = local_pair_support_allowed(pair, linear_modes, target_mode, cfg)
if strcmpi(cfg.support_modes_kind, 'linear_modes')
    tf = all(ismember(pair, linear_modes));
    return;
end
extra_count = sum(~ismember(pair, linear_modes));
if local_is_semi_pruned_v2(cfg)
    tf = extra_count <= 1 && any(pair == target_mode);
    return;
end
if local_is_semi_pruned(cfg)
    tf = extra_count <= 1;
    return;
end
tf = true;
end

function modes = local_term_modes(vars)
modes = cellfun(@(v) v.mode, vars);
end

function key = local_mode_pair_key(modes)
modes = sort(modes);
key = sprintf('%d,%d', modes(1), modes(2));
end

function [z3_h4, stats] = local_accumulate_w4_term(z3_h4, stats, vars, coeff, zlin, plin, linear_modes, target_outputs, nx, theta_grid, cfg)
if coeff == 0
    return;
end
stats.pruning_stats.terms_seen = stats.pruning_stats.terms_seen + 1;
physical_modes = cellfun(@(v) v.mode, vars);
if ~local_can_contribute_to_targets(physical_modes, linear_modes, target_outputs)
    stats.pruning_stats.terms_rejected_target_gate = stats.pruning_stats.terms_rejected_target_gate + 1;
    return;
end
for out = target_outputs
    deriv_indices = find(physical_modes == -out);
    for deriv_pos = deriv_indices
        stats.pruning_stats.derivative_slots_seen = stats.pruning_stats.derivative_slots_seen + 1;
        other_modes = physical_modes;
        other_modes(deriv_pos) = [];
        if any(~ismember(other_modes, linear_modes))
            stats.pruning_stats.derivative_slots_rejected_nonlinear_support = ...
                stats.pruning_stats.derivative_slots_rejected_nonlinear_support + 1;
            continue;
        end
        for deriv_kind = {'a', 'b'}
            stats.pruning_stats.normal_forms_seen = stats.pruning_stats.normal_forms_seen + 1;
            [vars_n, normal_coeff] = local_linear_state_normal_term(vars, coeff, deriv_pos, deriv_kind{1}, nx, theta_grid, cfg.g);
            if normal_coeff == 0
                stats.pruning_stats.normal_forms_zero_coeff = stats.pruning_stats.normal_forms_zero_coeff + 1;
                continue;
            end
            stats.candidate_terms = stats.candidate_terms + 1;
            weight = local_normal_weight(vars_n, nx, theta_grid, cfg.g);
            if abs(weight) < cfg.resonance_tol
                stats.pruning_stats.normal_forms_resonant = stats.pruning_stats.normal_forms_resonant + 1;
                stats.resonant_terms = stats.resonant_terms + 1;
                continue;
            end
            stats.pruning_stats.normal_forms_used = stats.pruning_stats.normal_forms_used + 1;
            stats.nonresonant_terms = stats.nonresonant_terms + 1;
            stats.used_terms = stats.used_terms + 1;
            wcoeff = cfg.homological_inverse_sign * normal_coeff / (1i * weight);
            deriv_val = local_eval_normal_derivative_at_index(vars_n, wcoeff, deriv_pos, zlin, plin, nx, theta_grid, cfg.g);
            if deriv_val ~= 0
                if isKey(z3_h4, out)
                    z3_h4(out) = z3_h4(out) + deriv_val;
                else
                    z3_h4(out) = deriv_val;
                end
            end
        end
    end
end
end

function [vars_n, coeff_n] = local_linear_state_normal_term(vars, coeff, deriv_pos, deriv_kind, nx, theta_grid, g)
vars_n = cell(1, numel(vars));
coeff_n = coeff;
for i = 1:numel(vars)
    v = vars{i};
    if i == deriv_pos
        kind = deriv_kind;
    elseif v.mode > 0
        kind = 'a';
    elseif v.mode < 0
        kind = 'b';
    else
        coeff_n = 0;
        return;
    end
    gamma = sqrt(local_theta(v.mode, nx, theta_grid) / g);
    if strcmp(v.kind, 'z')
        piece_coeff = 0.5;
    elseif gamma == 0
        piece_coeff = 0;
    elseif strcmp(kind, 'a')
        piece_coeff = 1/(2i * gamma);
    else
        piece_coeff = -1/(2i * gamma);
    end
    coeff_n = coeff_n * piece_coeff;
    if coeff_n == 0
        return;
    end
    vars_n{i} = local_nvar(kind, v.mode);
end
end

function [z3_h4, stats] = local_accumulate_w4_term_full_expand(z3_h4, stats, vars, coeff, zlin, plin, linear_modes, target_outputs, nx, theta_grid, cfg)
if coeff == 0
    return;
end
[normal_vars, normal_coeffs] = local_expand_to_normal(vars, coeff, nx, theta_grid, cfg.g);
for j = 1:numel(normal_coeffs)
    vars_n = normal_vars{j};
    stats.candidate_terms = stats.candidate_terms + 1;
    normal_modes = cellfun(@(v) v.mode, vars_n);
    for i = 1:numel(vars_n)
        out = -vars_n{i}.mode;
        if ~ismember(out, target_outputs)
            continue;
        end
        other_modes = normal_modes;
        other_modes(i) = [];
        if any(~ismember(other_modes, linear_modes))
            continue;
        end
        weight = local_normal_weight(vars_n, nx, theta_grid, cfg.g);
        if abs(weight) < cfg.resonance_tol
            stats.resonant_terms = stats.resonant_terms + 1;
            continue;
        end
        stats.nonresonant_terms = stats.nonresonant_terms + 1;
        stats.used_terms = stats.used_terms + 1;
        wcoeff = cfg.homological_inverse_sign * normal_coeffs(j) / (1i * weight);
        deriv_val = local_eval_normal_derivative_at_index(vars_n, wcoeff, i, zlin, plin, nx, theta_grid, cfg.g);
        if deriv_val ~= 0
            if isKey(z3_h4, out)
                z3_h4(out) = z3_h4(out) + deriv_val;
            else
                z3_h4(out) = deriv_val;
            end
        end
    end
end
end

function tf = local_can_contribute_to_targets(term_modes, linear_modes, target_outputs)
tf = false;
for target = target_outputs
    hit_idx = find(term_modes == -target);
    for j = 1:numel(hit_idx)
        other_modes = term_modes;
        other_modes(hit_idx(j)) = [];
        if all(ismember(other_modes, linear_modes))
            tf = true;
            return;
        end
    end
end
end

function out = local_target_normal_coeffs(K4, nx, theta_grid, cfg)
out = struct('available', false);
if ~isfield(cfg, 'debug_target_mode') || isempty(cfg.debug_target_mode)
    return;
end
target = cfg.debug_target_mode;
target_a = local_normal_key({local_nvar('a', -3 * target), ...
    local_nvar('a', target), local_nvar('a', target), local_nvar('a', target)});
target_b = local_normal_key({local_nvar('b', -3 * target), ...
    local_nvar('a', target), local_nvar('a', target), local_nvar('a', target)});
k4_a = 0;
k4_b = 0;
w4_a = 0;
w4_b = 0;
terms = local_poly_terms(K4);
for i = 1:numel(terms)
    [normal_vars, normal_coeffs] = local_expand_to_normal(terms(i).vars, terms(i).coeff, nx, theta_grid, cfg.g);
    for j = 1:numel(normal_coeffs)
        key = local_normal_key(normal_vars{j});
        weight = local_normal_weight(normal_vars{j}, nx, theta_grid, cfg.g);
        wcoeff = 0;
        if abs(weight) >= cfg.resonance_tol
            wcoeff = cfg.homological_inverse_sign * normal_coeffs(j) / (1i * weight);
        end
        if strcmp(key, target_a)
            k4_a = k4_a + normal_coeffs(j);
            w4_a = w4_a + wcoeff;
        elseif strcmp(key, target_b)
            k4_b = k4_b + normal_coeffs(j);
            w4_b = w4_b + wcoeff;
        end
    end
end
out = struct( ...
    'available', true, ...
    'target_mode', target, ...
    'k4_a_minus3_a1_cubed', k4_a, ...
    'k4_b_minus3_a1_cubed', k4_b, ...
    'w4_a_minus3_a1_cubed', w4_a, ...
    'w4_b_minus3_a1_cubed', w4_b);
end

function key = local_normal_key(vars)
parts = cell(1, numel(vars));
for i = 1:numel(vars)
    parts{i} = sprintf('%s%d', vars{i}.kind, vars{i}.mode);
end
parts = sort(parts);
key = strjoin(parts, '|');
end

function [normal_vars, normal_coeffs] = local_expand_to_normal(vars, coeff, nx, theta_grid, g)
normal_vars = {{}};
normal_coeffs = coeff;
for i = 1:numel(vars)
    v = vars{i};
    gamma = sqrt(local_theta(v.mode, nx, theta_grid) / g);
    if strcmp(v.kind, 'z')
        pieces = {local_nvar('a', v.mode), local_nvar('b', v.mode)};
        coeffs = [0.5, 0.5];
    else
        if gamma == 0
            pieces = {local_nvar('a', v.mode), local_nvar('b', v.mode)};
            coeffs = [0, 0];
        else
            pieces = {local_nvar('a', v.mode), local_nvar('b', v.mode)};
            coeffs = [1/(2i*gamma), -1/(2i*gamma)];
        end
    end
    new_vars = {};
    new_coeffs = [];
    for j = 1:numel(normal_coeffs)
        for p = 1:2
            new_vars{end+1} = [normal_vars{j}, pieces(p)]; %#ok<AGROW>
            new_coeffs(end+1) = normal_coeffs(j) * coeffs(p); %#ok<AGROW>
        end
    end
    normal_vars = new_vars;
    normal_coeffs = new_coeffs;
end
end

function weight = local_normal_weight(vars, nx, theta_grid, g)
weight = 0;
for i = 1:numel(vars)
    om = sqrt(g * local_theta(vars{i}.mode, nx, theta_grid));
    if strcmp(vars{i}.kind, 'a')
        weight = weight + om;
    else
        weight = weight - om;
    end
end
end

function val = local_eval_normal_derivative(vars, coeff, deriv_mode, zlin, plin, nx, theta_grid, g)
val = 0;
for i = 1:numel(vars)
    if vars{i}.mode ~= deriv_mode
        continue;
    end
    val = val + local_eval_normal_derivative_at_index(vars, coeff, i, zlin, plin, nx, theta_grid, g);
end
end

function val = local_eval_normal_derivative_at_index(vars, coeff, deriv_index, zlin, plin, nx, theta_grid, g)
deriv_mode = vars{deriv_index}.mode;
gamma = sqrt(local_theta(deriv_mode, nx, theta_grid) / g);
if strcmp(vars{deriv_index}.kind, 'a')
    dval = 1i * gamma;
else
    dval = -1i * gamma;
end
val = coeff * dval;
for j = 1:numel(vars)
    if j == deriv_index
        continue;
    end
    val = val * local_eval_normal_var(vars{j}, zlin, plin, nx, theta_grid, g);
    if val == 0
        break;
    end
end
end

function val = local_eval_normal_var(v, zlin, plin, nx, theta_grid, g)
zv = local_z(zlin, v.mode);
pv = local_z(plin, v.mode);
gamma = sqrt(local_theta(v.mode, nx, theta_grid) / g);
if strcmp(v.kind, 'a')
    val = zv + 1i * gamma * pv;
else
    val = zv - 1i * gamma * pv;
end
end

function br = local_poisson_bracket(F, G, mode_set)
br = containers.Map('KeyType', 'char', 'ValueType', 'any');
for k = mode_set
    f_z = local_poly_derivative(F, local_var('z', k));
    if f_z.Count > 0
        g_p = local_poly_derivative(G, local_var('p', -k));
        br = local_poly_add(br, local_poly_multiply(f_z, g_p), 1);
    end
    f_p = local_poly_derivative(F, local_var('p', k));
    if f_p.Count > 0
        g_z = local_poly_derivative(G, local_var('z', -k));
        br = local_poly_add(br, local_poly_multiply(f_p, g_z), -1);
    end
end
end

function out = local_poly_derivative(poly, target)
out = containers.Map('KeyType', 'char', 'ValueType', 'any');
terms = local_poly_terms(poly);
for i = 1:numel(terms)
    [new_vars, count] = local_derivative_term(terms(i).vars, target);
    if count ~= 0
        out = local_add_poly_term(out, new_vars, terms(i).coeff * count);
    end
end
end

function out = local_poly_multiply(A, B)
out = containers.Map('KeyType', 'char', 'ValueType', 'any');
if A.Count == 0 || B.Count == 0
    return;
end
terms_a = local_poly_terms(A);
terms_b = local_poly_terms(B);
for i = 1:numel(terms_a)
    for j = 1:numel(terms_b)
        out = local_add_poly_term(out, [terms_a(i).vars, terms_b(j).vars], ...
            terms_a(i).coeff * terms_b(j).coeff);
    end
end
end

function [new_vars, count] = local_derivative_term(vars, target)
matches = false(1, numel(vars));
for i = 1:numel(vars)
    matches(i) = strcmp(vars{i}.kind, target.kind) && vars{i}.mode == target.mode;
end
count = sum(matches);
if count == 0
    new_vars = {};
    return;
end
first = find(matches, 1);
new_vars = vars;
new_vars(first) = [];
end

function out = local_poly_add(A, B, scaleB)
out = containers.Map('KeyType', 'char', 'ValueType', 'any');
keysA = keys(A);
for i = 1:numel(keysA)
    out(keysA{i}) = A(keysA{i});
end
keysB = keys(B);
for i = 1:numel(keysB)
    key = keysB{i};
    if isKey(out, key)
        out(key) = out(key) + scaleB * B(key);
    else
        out(key) = scaleB * B(key);
    end
end
end

function out = local_poly_scale(A, scaleA)
out = containers.Map('KeyType', 'char', 'ValueType', 'any');
keysA = keys(A);
for i = 1:numel(keysA)
    out(keysA{i}) = scaleA * A(keysA{i});
end
end

function n = local_poly_count(poly)
n = poly.Count;
end

function terms = local_poly_terms(poly)
keys_poly = keys(poly);
terms = repmat(struct('vars', {{}}, 'coeff', 0), 1, numel(keys_poly));
for i = 1:numel(keys_poly)
    terms(i).vars = local_key_to_vars(keys_poly{i});
    terms(i).coeff = poly(keys_poly{i});
end
end

function poly = local_add_poly_term(poly, vars, coeff)
if coeff == 0
    return;
end
key = local_vars_to_key(vars);
if isKey(poly, key)
    poly(key) = poly(key) + coeff;
else
    poly(key) = coeff;
end
end

function key = local_vars_to_key(vars)
parts = cell(1, numel(vars));
for i = 1:numel(vars)
    parts{i} = sprintf('%s%d', vars{i}.kind, vars{i}.mode);
end
parts = sort(parts);
key = strjoin(parts, '|');
end

function vars = local_key_to_vars(key)
if isempty(key)
    vars = {};
    return;
end
parts = strsplit(key, '|');
vars = cell(1, numel(parts));
for i = 1:numel(parts)
    vars{i} = local_parse_var(parts{i});
end
end

function v = local_parse_var(part)
v = struct('kind', part(1), 'mode', str2double(part(2:end)));
end

function v = local_var(kind, mode)
v = struct('kind', kind, 'mode', mode);
end

function v = local_nvar(kind, mode)
v = struct('kind', kind, 'mode', mode);
end

function [D, Bz, Bphi] = local_b_d(dest, a, b, nx, k_grid, theta_grid, g)
theta_dest = local_theta(dest, nx, theta_grid);
theta_a = local_theta(a, nx, theta_grid);
theta_b = local_theta(b, nx, theta_grid);
k_dest_sq = local_k(dest, nx, k_grid)^2;
k_a_sq = local_k(a, nx, k_grid)^2;
k_b_sq = local_k(b, nx, k_grid)^2;
D = local_kernel_D_finite(theta_dest, theta_a, theta_b, k_dest_sq, k_a_sq, k_b_sq, g);
Bz = local_kernel_B_finite(theta_a, theta_b, theta_dest, k_a_sq, k_b_sq, k_dest_sq);
Bphi = local_kernel_B_finite(theta_dest, theta_a, theta_b, k_dest_sq, k_a_sq, k_b_sq);
end

function val = local_theta(m, nx, theta_grid)
if abs(m) > nx/2
    val = 0;
else
    val = theta_grid(local_mode_to_index(m, nx));
end
end

function val = local_k(m, nx, k_grid)
if abs(m) > nx/2
    val = 0;
else
    val = k_grid(local_mode_to_index(m, nx));
end
end

function d_val = local_kernel_D_finite(theta1, theta2, theta3, k1sq, k2sq, k3sq, g)
den_base = theta1 * (theta2 + theta3 - theta1) ...
         + theta2 * (theta1 + theta3 - theta2) ...
         + theta3 * (theta1 + theta2 - theta3);
if abs(den_base) < 1e-14 || theta1 == 0 || theta2 == 0 || theta3 == 0
    d_val = 0;
    return;
end
num = k1sq * (theta1^2 - (theta2 - theta3)^2) ...
    + k2sq * (theta2^2 - (theta1 - theta3)^2) ...
    + k3sq * (theta3^2 - (theta1 - theta2)^2) ...
    - 2 * theta1 * theta2 * theta3 * (theta1 + theta2 + theta3);
d_val = num / (12 * g * den_base);
end

function b_val = local_kernel_B_finite(theta1, theta2, theta3, k1sq, k2sq, k3sq)
den_base = theta1 * (theta2 + theta3 - theta1) ...
         + theta2 * (theta1 + theta3 - theta2) ...
         + theta3 * (theta1 + theta2 - theta3);
if abs(den_base) < 1e-14
    b_val = 0;
    return;
end
num = theta3 * (theta3 * (theta1 + theta2) - theta1^2 - theta2^2) ...
    + theta3 * (k1sq + k2sq - 2 * k3sq) ...
    + (theta1 - theta2) * (k1sq - k2sq);
b_val = num / (2 * den_base);
end
