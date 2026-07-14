function p = paths_local(p)
%PATHS_LOCAL  Per-machine path overrides (MATLAB) -- EXAMPLE TEMPLATE.
%
%   Copy this file to config/paths_local.m (same folder) and edit. paths_local.m
%   is gitignored, so each machine can point at local data / toolboxes without
%   touching tracked code. lab_paths() calls this LAST, so anything you set here
%   wins. Only override what differs on this machine.
%
%   The function receives the fully-derived struct `p` and must return it.

    % Example: data on a fast local drive on the acquisition PC.
    % p.data_root = 'D:\Imaging_Data';

    % Example: prepend a machine-local toolbox copy so it is found first.
    % p.toolboxes = [ {'C:\Codes'}; p.toolboxes ];
end
