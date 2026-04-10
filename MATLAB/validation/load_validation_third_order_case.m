function val = load_validation_third_order_case(validation_file, alpha, chi_deg, spread_deg, kd_target, akp_target)
%LOAD_VALIDATION_THIRD_ORDER_CASE Read one third-order validation case from
%the HDF5 MAT file.

case_path = sprintf('/out/third_order/data/akp_%s_alpha_%s_chi_%s_spread_%s_kd_%s', ...
    scalar_key(akp_target), scalar_key(alpha), scalar_key(chi_deg), ...
    scalar_key(spread_deg), scalar_key(kd_target));

val = struct();
val.case_path = case_path;
val.x = double(h5read(validation_file, [case_path '/x']));
val.y = double(h5read(validation_file, [case_path '/y']));
val.kp = double(h5read(validation_file, [case_path '/kp']));
val.kd = double(h5read(validation_file, [case_path '/kd']));
val.akp = double(h5read(validation_file, [case_path '/akp']));
val.alpha = double(h5read(validation_file, [case_path '/alpha']));
val.direction_deg = double(h5read(validation_file, [case_path '/direction_deg']));
val.spread_deg = double(h5read(validation_file, [case_path '/spread_deg']));
val.eta_linear_input = double(h5read(validation_file, [case_path '/eta_linear_input']));
val.eta31 = double(h5read(validation_file, [case_path '/eta31_mf12']));
val.eta33 = double(h5read(validation_file, [case_path '/eta33_mf12']));
val.source_file = char(h5read(validation_file, [case_path '/source_file'])).';
end

function text = scalar_key(value)
text = strrep(sprintf('%.6g', value), '.', 'p');
text = strrep(text, '-', 'm');
end
