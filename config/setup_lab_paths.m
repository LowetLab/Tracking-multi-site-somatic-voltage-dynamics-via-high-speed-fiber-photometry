function p = setup_lab_paths(verbose)
%SETUP_LAB_PATHS  Add all external toolbox folders to the MATLAB path.
%
%   setup_lab_paths()         adds the dependency folders defined in lab_paths()
%                             to the MATLAB path (silently skips missing ones).
%   setup_lab_paths(true)     also warns about any dependency folder not found.
%   p = setup_lab_paths(...)  additionally returns the lab_paths() struct.
%
%   This REPLACES per-script blocks of hardcoded
%       addpath(genpath('C:\some\absolute\path'))
%       ...
%   Call it once at the top of any entry-point script:
%       setup_lab_paths();
%
%   See also LAB_PATHS.

    if nargin < 1 || isempty(verbose)
        verbose = false;
    end

    p = lab_paths();

    for i = 1:numel(p.toolboxes)
        tb = p.toolboxes{i};
        if isfolder(tb)
            addpath(genpath(tb));
        elseif verbose
            warning('lab_paths:missingToolbox', ...
                    'Dependency folder not found (skipped): %s', tb);
        end
    end
end
