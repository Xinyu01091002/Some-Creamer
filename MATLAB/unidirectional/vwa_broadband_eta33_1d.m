function out = vwa_broadband_eta33_1d(eta11, x_vec, h, g, cfg)
%VWA_BROADBAND_ETA33_1D Bridge to the external broadband VWA eta33 helper.
%
% Uses the existing VWA project at:
%   C:\Research\VWA\VWA Unidirectinal\Unidirectional
%
% and calls its opensource wrapper:
%   vwa_eval_from_opensource(eta11, x_vec, h, g, 3)

if nargin < 5 || isempty(cfg)
    cfg = struct();
end
cfg = local_defaults(cfg);

eta11 = reshape(real(eta11), [], 1);
x_vec = x_vec(:);
nx = numel(eta11);

if numel(x_vec) ~= nx
    error('vwa_broadband_eta33_1d:SizeMismatch', ...
        'Length of x_vec (%d) must match eta length (%d).', numel(x_vec), nx);
end

persistent vwa_path_added;
if isempty(vwa_path_added) || ~vwa_path_added
    if exist(cfg.vwa_root_dir, 'dir') ~= 7
        error('vwa_broadband_eta33_1d:MissingVwaRoot', ...
            'VWA root directory not found: %s', cfg.vwa_root_dir);
    end

    vwa_wrapper = fullfile(cfg.vwa_root_dir, 'vwa_eval_from_opensource.m');
    if exist(vwa_wrapper, 'file') ~= 2
        error('vwa_broadband_eta33_1d:MissingWrapper', ...
            'VWA wrapper not found: %s', vwa_wrapper);
    end

    addpath(cfg.vwa_root_dir, '-begin');
    clear('vwa_eval_from_opensource');
    vwa_path_added = true;
end

[eta_vwa, phi_vwa] = vwa_eval_from_opensource(eta11, x_vec, h, g, 3);

if ~iscell(eta_vwa) || numel(eta_vwa) < 3 || isempty(eta_vwa{3})
    error('vwa_broadband_eta33_1d:MissingEta33', ...
        'External VWA wrapper did not return eta_vwa{3}.');
end

eta33 = reshape(real(eta_vwa{3}), 1, []);

out = struct();
out.eta33 = eta33;
out.eta33_hat = fft(eta33) / nx;
out.phi33 = [];
if iscell(phi_vwa) && numel(phi_vwa) >= 3 && ~isempty(phi_vwa{3})
    out.phi33 = reshape(real(phi_vwa{3}), 1, []);
end
out.diagnostics = struct( ...
    'backend', 'external_vwa_opensource', ...
    'vwa_root_dir', cfg.vwa_root_dir);
end

function cfg = local_defaults(cfg)
defaults = struct( ...
    'vwa_root_dir', 'C:\Research\VWA\VWA Unidirectinal\Unidirectional');

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(cfg, names{i}) || isempty(cfg.(names{i}))
        cfg.(names{i}) = defaults.(names{i});
    end
end
end
