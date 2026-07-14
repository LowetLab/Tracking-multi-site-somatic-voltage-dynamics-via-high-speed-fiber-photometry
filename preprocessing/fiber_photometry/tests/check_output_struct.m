function ok = check_output_struct(mat_file)
%CHECK_OUTPUT_STRUCT  Sanity-check a *_FiberPhotometry_Analysis.mat output.
%   ok = CHECK_OUTPUT_STRUCT(mat_file) loads a saved FiberPhotometryAnalysis
%   struct and verifies that the fields fed by the Tier-B step-4 extractions
%   (photobleaching, ΔF/F, stim detection) are present, non-empty and finite.
%
%   This is the INTEGRATION check the unit tests can't do: the unit tests prove
%   the extracted functions are numerically identical to the old inline code,
%   but only a real run confirms the script WIRING still populates the saved
%   struct (e.g. the photobleaching fit_diag values that used to be loop
%   leftovers and are restored with isfield/exist guards). A clean run that
%   passes this check means no field silently went missing or empty.
%
%   Run after a real-recording run:
%       addpath('tests')
%       check_output_struct('...\<MOUSE_NAME>-..._Trial1_FiberPhotometry_Analysis.mat')
%
%   Prints a per-field PASS/FAIL report and returns true iff all required
%   checks pass. Optional fields (only saved for stim trials / when a fit
%   succeeded) are reported as SKIP when absent rather than failing.

S = load(mat_file);
assert(isfield(S, 'FiberPhotometryAnalysis'), ...
    'check_output_struct:noStruct', 'No FiberPhotometryAnalysis variable in %s', mat_file);
A = S.FiberPhotometryAnalysis;

fprintf('\n=== check_output_struct ===\n');
fprintf('File: %s\n\n', mat_file);
nfail = 0;

% --- stage 6: ΔF/F + F0 (always expected) -------------------------------
nfail = nfail + req(A, {'signals','raw_traces'},        'raw ROI traces');
nfail = nfail + req(A, {'signals','corrected_traces'},  'photobleach-corrected traces (stage 3/4)');
nfail = nfail + req(A, {'signals','deltaF_F_traces'},   'ΔF/F traces (stage 6)');
nfail = nfail + req(A, {'signals','zscored_traces'},    'z-scored ΔF/F (stage 6)');
nfail = nfail + req(A, {'signals','final_processed_traces'}, 'final processed traces');
nfail = nfail + req(A, {'signals','F0_values'},         'per-fiber F0 (stage 6)');

% Cross-check shapes line up across the trace fields.
nfail = nfail + same_size(A, {'signals','corrected_traces'}, {'signals','deltaF_F_traces'}, ...
    'corrected vs ΔF/F trace size');
nfail = nfail + same_size(A, {'signals','deltaF_F_traces'}, {'signals','zscored_traces'}, ...
    'ΔF/F vs z-scored trace size');

% --- stage 3/4: photobleaching fit_diag (only if a fit succeeded) -------
if has(A, {'photobleaching','fit_success'}) && get(A, {'photobleaching','fit_success'})
    fprintf('  (photobleaching fit_success = true -> checking fit_diag fields)\n');
    nfail = nfail + req(A, {'photobleaching','fit_params'},   'fit params (fit_diag, stage 3/4)');
    nfail = nfail + req(A, {'photobleaching','time_full'},    'time_full (fit_diag, stage 3/4)');
    nfail = nfail + req(A, {'photobleaching','fitted_curve'}, 'fitted_curve (fit_diag, stage 3/4)');
    nfail = nfail + opt(A, {'photobleaching','time_pre_stim'},  'time_pre_stim (fit_diag)');
    nfail = nfail + opt(A, {'photobleaching','pre_stim_trace'}, 'pre_stim_trace (fit_diag)');
else
    fprintf('  SKIP photobleaching fit_diag (no successful double-exp fit recorded)\n');
end

% --- stage 5/5b: stim detection (only on stim trials) -------------------
% The onset can live under a couple of names across versions; accept any.
stim_found = opt_any(A, { {'time_periods','stim_onset_frame'}, ...
                          {'stimulation','onset_sample'}, ...
                          {'stimulation','stim_onset_sample'}, ...
                          {'metadata','stim_onset_ephys_sample'} }, ...
                          'stim onset (stage 5/5b)');
if ~stim_found
    fprintf('  SKIP stim onset (baseline trial or stored under another name)\n');
end

fprintf('\n');
if nfail == 0
    fprintf('RESULT: OUTPUT STRUCT OK  (all required fields present, non-empty, finite)\n');
    ok = true;
else
    fprintf('RESULT: %d REQUIRED CHECK(S) FAILED -- send me the list above\n', nfail);
    ok = false;
end
end

% =========================================================================
function nf = req(A, path, label)
% Required: must exist, be non-empty, and contain no NaN/Inf in numeric data.
if ~has(A, path)
    fprintf('  FAIL %-45s : MISSING (%s)\n', label, dotpath(path)); nf = 1; return;
end
v = get(A, path);
if isempty(v)
    fprintf('  FAIL %-45s : EMPTY\n', label); nf = 1; return;
end
if isnumeric(v) && ~all(isfinite(v(:)))
    fprintf('  FAIL %-45s : contains NaN/Inf\n', label); nf = 1; return;
end
fprintf('  PASS %-45s : %s\n', label, describe(v)); nf = 0;
end

function nf = opt(A, path, label)
% Optional: report SKIP when absent, but validate when present.
if ~has(A, path)
    fprintf('  SKIP %-45s : not present\n', label); nf = 0; return;
end
nf = req(A, path, label);
end

function found = opt_any(A, paths, label)
% Optional, accept the first of several possible field names.
for k = 1:numel(paths)
    if has(A, paths{k})
        opt(A, paths{k}, label);
        found = true; return;
    end
end
found = false;
end

function nf = same_size(A, p1, p2, label)
if ~has(A, p1) || ~has(A, p2)
    fprintf('  SKIP %-45s : a field is missing\n', label); nf = 0; return;
end
if isequal(size(get(A, p1)), size(get(A, p2)))
    fprintf('  PASS %-45s : %s\n', label, mat2str(size(get(A, p1)))); nf = 0;
else
    fprintf('  FAIL %-45s : %s vs %s\n', label, mat2str(size(get(A, p1))), mat2str(size(get(A, p2)))); nf = 1;
end
end

% --- small struct-path helpers ------------------------------------------
function tf = has(A, path)
tf = true;
for k = 1:numel(path)
    if isstruct(A) && isfield(A, path{k})
        A = A.(path{k});
    else
        tf = false; return;
    end
end
end

function v = get(A, path)
for k = 1:numel(path), A = A.(path{k}); end
v = A;
end

function s = dotpath(path)
s = strjoin(path, '.');
end

function s = describe(v)
if isnumeric(v)
    s = sprintf('numeric %s', mat2str(size(v)));
elseif ischar(v)
    s = sprintf('char ''%s''', v);
elseif islogical(v)
    s = sprintf('logical %s', mat2str(v));
else
    s = class(v);
end
end
