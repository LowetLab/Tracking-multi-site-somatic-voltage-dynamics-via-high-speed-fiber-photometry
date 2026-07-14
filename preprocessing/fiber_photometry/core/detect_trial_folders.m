function [trial_folders, trial_numbers, num_trials] = detect_trial_folders(base_folder)
%DETECT_TRIAL_FOLDERS  Find and sort the per-trial subfolders of a recording.
%   [trial_folders, trial_numbers, num_trials] = DETECT_TRIAL_FOLDERS(base_folder)
%   scans base_folder for subfolders whose name ends in "_<number>" (e.g.
%   recording_1, recording_2, recording_10), and returns them sorted by that
%   trailing number.
%
%   Inputs
%   ------
%   base_folder : path to the recording folder containing the trial subfolders.
%
%   Outputs
%   -------
%   trial_folders : 1xN cell array of folder NAMES (not full paths), ascending
%                   by trial number.
%   trial_numbers : 1xN double array of the trailing trial numbers (sorted).
%   num_trials    : N (number of trial folders found).
%
%   Errors if base_folder does not exist, is empty, or contains no "_<number>"
%   subfolders -- with the same messages the inline version used.
%
%   Extracted verbatim (behaviour-preserving) from the multi-trial fiber
%   preprocessing script. Unit-tested in core/tests/test_detect_trial_folders.m.

% Verify base folder exists
if ~exist(base_folder, 'dir')
    error('Base folder does not exist: %s\nPlease check MOUSE_NAME, RECORDING_DATE, and RECORDING_ID parameters.', base_folder);
end

% Find all trial folders with _1, _2, _3 suffixes
fprintf('Searching for multi-trial folders in: %s\n', base_folder);

folder_contents = dir(base_folder);
folder_contents = folder_contents([folder_contents.isdir]);  % Keep only directories
folder_contents = folder_contents(~ismember({folder_contents.name}, {'.', '..'}));  % Remove . and ..

if isempty(folder_contents)
    error('No trial folders found in: %s', base_folder);
end

% Extract trial numbers from folder names (look for _1, _2, _3, etc.)
trial_folders = {};
trial_numbers = [];

for i = 1:length(folder_contents)
    folder_name = folder_contents(i).name;
    % Look for pattern ending with _ followed by digits
    match = regexp(folder_name, '_(\d+)$', 'tokens');
    if ~isempty(match)
        trial_num = str2double(match{1}{1});
        trial_folders{end+1} = folder_name; %#ok<AGROW>
        trial_numbers(end+1) = trial_num; %#ok<AGROW>
    end
end

if isempty(trial_folders)
    error('No trial folders with _1, _2, _3 suffixes found in: %s\nPlease ensure trial folders are named with _1, _2, _3, etc. suffixes.', base_folder);
end

% Sort by trial number
[trial_numbers, sort_idx] = sort(trial_numbers);
trial_folders = trial_folders(sort_idx);
num_trials = length(trial_folders);

fprintf('Found %d trial folders:\n', num_trials);
for i = 1:num_trials
    fprintf('  Trial %d: %s\n', trial_numbers(i), trial_folders{i});
end
end
