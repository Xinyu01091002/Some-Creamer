function [best_case, match_summary] = match_validation_second_order_case(validation_file, alpha_candidates, chi_deg, spread_deg, kd_target, eta_ref, x_ref, y_ref, akp_target)
%MATCH_VALIDATION_SECOND_ORDER_CASE Choose the MF12 validation case whose
%stored linear input best matches the current Creamer parent surface.

if nargin < 2 || isempty(alpha_candidates)
    error('match_validation_second_order_case:MissingAlphaCandidates', ...
        'alpha_candidates must be provided.');
end

x_ref = x_ref(:).';
y_ref = y_ref(:).';
[X_ref, Y_ref] = meshgrid(x_ref, y_ref);

best_idx = 0;
best_score = inf;
cases = cell(1, numel(alpha_candidates));
scores = nan(1, numel(alpha_candidates));

for idx = 1:numel(alpha_candidates)
    if nargin >= 9
        this_case = load_validation_second_order_case( ...
            validation_file, alpha_candidates(idx), chi_deg, spread_deg, kd_target, akp_target);
    else
        this_case = load_validation_second_order_case( ...
            validation_file, alpha_candidates(idx), chi_deg, spread_deg, kd_target);
    end

    [X_val, Y_val] = meshgrid(this_case.x(:).', this_case.y(:));
    eta_interp = interp2(X_val, Y_val, this_case.eta_linear_input, X_ref, Y_ref, 'linear', 0);

    rel_l2 = norm(eta_interp(:) - eta_ref(:)) / max(norm(eta_ref(:)), eps);
    amp_ratio = max(abs(eta_interp), [], 'all') / max(max(abs(eta_ref), [], 'all'), eps);
    corr_mat = corrcoef(eta_interp(:), eta_ref(:));
    corr_val = corr_mat(1, 2);

    this_case.eta_linear_input_interp = eta_interp;
    this_case.match_rel_l2 = rel_l2;
    this_case.match_amp_ratio = amp_ratio;
    this_case.match_corr = corr_val;

    cases{idx} = this_case;
    scores(idx) = rel_l2;

    if rel_l2 < best_score
        best_score = rel_l2;
        best_idx = idx;
    end
end

best_case = cases{best_idx};
match_summary = struct();
match_summary.alpha_candidates = alpha_candidates(:).';
match_summary.scores_rel_l2 = scores;
match_summary.selected_alpha = best_case.alpha;
match_summary.selected_case_path = best_case.case_path;
match_summary.selected_rel_l2 = best_case.match_rel_l2;
match_summary.selected_amp_ratio = best_case.match_amp_ratio;
match_summary.selected_corr = best_case.match_corr;
if nargin >= 9
    match_summary.selected_akp = akp_target;
end
match_summary.cases = cases;
end
