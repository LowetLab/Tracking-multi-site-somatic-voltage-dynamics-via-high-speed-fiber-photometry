function [ephys_path, LOAD_EPHYS_DATA] = prompt_manual_ephys_file_selection()
%PROMPT_MANUAL_EPHYS_FILE_SELECTION  Fallback UI prompt for locating Open Ephys data.
%   [ephys_path, LOAD_EPHYS_DATA] = PROMPT_MANUAL_EPHYS_FILE_SELECTION() asks
%   the user to manually pick an Open Ephys .continuous file when automatic
%   path discovery has failed. Returns LOAD_EPHYS_DATA = false if the user
%   cancels the dialog (caller should also set its own EPHYS_LOADED = false
%   in that case, since local functions cannot modify caller-workspace
%   variables that were not passed in).
%
%   Shared by the single- and multi-trial preprocessing scripts, which each
%   called an identical "prompt, uigetfile, bail out on cancel" block at
%   every automatic-discovery failure point (missing base folder, missing
%   recording folder, missing Record Node folder, missing .continuous files).

fprintf('Prompting user to select file manually...\n');
[file, ephys_path] = uigetfile('.continuous', 'Select any Open Ephys .continuous file');
LOAD_EPHYS_DATA = true;
if file == 0
    warning('No Open Ephys file selected. Skipping ephys analysis.');
    LOAD_EPHYS_DATA = false;
end
end
