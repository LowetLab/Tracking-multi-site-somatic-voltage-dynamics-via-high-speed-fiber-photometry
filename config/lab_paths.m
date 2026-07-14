function p = lab_paths()
%LAB_PATHS  Canonical paths for this pipeline (single source of truth).
%
%   p = LAB_PATHS() returns a struct of all project paths, derived by
%   self-location from this file. Nothing is hardcoded to a specific drive,
%   so moving / copying the project tree (preserving its internal structure)
%   just works.
%
%   Two anchors:
%     * project_root  -- SELF-LOCATED from this file (works from a clone
%                        anywhere on disk). Code and figure outputs live
%                        under here.
%     * lab_root      -- an EXPLICIT default you must edit (see below). Your
%                        data and external MATLAB toolboxes live OUTSIDE the
%                        repo, so they cannot be derived from project_root.
%                        Change it per machine via paths_local.m.
%
%   PER-MACHINE OVERRIDES:
%   To override any field on a specific machine (e.g. data on a local C:/D:
%   drive, or a different share), copy config/paths_local.example.m to
%   config/paths_local.m and edit it. paths_local.m is gitignored and, if
%   present, is applied last.
%
%   See also SETUP_LAB_PATHS.

    this_dir = fileparts(mfilename('fullpath'));   % <project_root>/config
    p.project_root = up(this_dir, 1);              % <project_root>  (self-located)

    % Lab/data share root -- the one explicit anchor for everything outside the
    % repo. EDIT THIS to point at your own data location, or (preferred) leave
    % it and set p.data_root in a local config/paths_local.m override instead.
    p.lab_root  = 'C:\PATH\TO\YOUR\DATA_SHARE';

    % --- Data & outputs ------------------------------------------------------
    p.data_root      = fullfile(p.lab_root, 'Data');   % raw imaging + ephys data root
    p.figures_root   = fullfile(p.project_root, 'Figures');
    p.spectral_output_root = fullfile(p.figures_root, 'Spectral_data_outputs');

    % --- Internal code roots -------------------------------------------------
    p.matlab_codes   = p.project_root;
    p.spectral_code  = fullfile(p.matlab_codes, 'spectral_analysis');

    % --- External MATLAB toolboxes / dependency script folders ---------------
    % These are NOT bundled with this repo -- install them separately and point
    % at your local copies here (or override in paths_local.m). Missing entries
    % are skipped (with an optional warning), so a machine without all of them
    % still runs whatever does not depend on them. See environment/SETUP.md.
    p.toolboxes = {
        'C:\Toolboxes\fieldtrip';        % FieldTrip (github.com/fieldtrip/fieldtrip) -- coherence
        'C:\Toolboxes\NoRMCorre';        % NoRMCorre (optional motion correction)
        'C:\Toolboxes\open-ephys-matlab-tools';  % provides load_open_ephys_data.m
        'C:\Toolboxes\spike_detection';  % your spike-detection scripts, see preprocessing/cellular_imaging/README.md
    };

    % --- Per-machine overrides (optional, gitignored) ------------------------
    if exist('paths_local', 'file') == 2
        p = paths_local(p);
    end
end

function d = up(d, n)
%UP  Go up N directory levels.
    for i = 1:n
        d = fileparts(d);
    end
end
