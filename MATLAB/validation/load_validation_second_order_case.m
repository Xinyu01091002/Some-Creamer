function val = load_validation_second_order_case(validation_file, alpha, chi_deg, spread_deg, kd_target, akp_target)
%LOAD_VALIDATION_SECOND_ORDER_CASE Read one validation case from the HDF5 MAT file.

if nargin >= 6 && ~isempty(akp_target)
    case_path = sprintf('/out/second_order/data/akp_%s_alpha_%s_chi_%s_spread_%s_kd_%s', ...
        scalar_key(akp_target), scalar_key(alpha), scalar_key(chi_deg), ...
        scalar_key(spread_deg), scalar_key(kd_target));
else
    case_path = sprintf('/out/second_order/data/alpha_%g_chi_%g_spread_%g_kd_%g', ...
        alpha, chi_deg, spread_deg, kd_target);
end

val = struct();
val.case_path = case_path;
val.x = double(h5read(validation_file, [case_path '/x']));
val.y = double(h5read(validation_file, [case_path '/y']));
val.kp = double(h5read(validation_file, [case_path '/kp']));
val.kd = double(h5read(validation_file, [case_path '/kd']));
val.alpha = double(h5read(validation_file, [case_path '/alpha']));
val.direction_deg = double(h5read(validation_file, [case_path '/direction_deg']));
val.spread_deg = double(h5read(validation_file, [case_path '/spread_deg']));
val.eta_linear_input = double(h5read(validation_file, [case_path '/eta_linear_input']));
val.eta20 = double(h5read(validation_file, [case_path '/eta20_mf12']));
val.eta22 = double(h5read(validation_file, [case_path '/eta22_mf12']));
val.eta2_total = double(h5read(validation_file, [case_path '/eta2_total_mf12']));
val.phi20 = double(h5read(validation_file, [case_path '/phi20_mf12']));
val.phi22 = double(h5read(validation_file, [case_path '/phi22_mf12']));
val.source_file = char(h5read(validation_file, [case_path '/source_file'])).';
end

function text = scalar_key(value)
text = strrep(sprintf('%.6g', value), '.', 'p');
text = strrep(text, '-', 'm');
end
