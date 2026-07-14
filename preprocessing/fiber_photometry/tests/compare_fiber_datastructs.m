function report = compare_fiber_datastructs(a, b, varargin)
%COMPARE_FIBER_DATASTRUCTS  Recursively compare two preprocessing outputs.
%   Tells you -- field by field -- whether a refactor changed the saved output,
%   so validating a refactor is "run this, read PASS/FAIL" instead of eyeballing
%   structs by hand.
%
%   report = COMPARE_FIBER_DATASTRUCTS(a, b)
%   report = COMPARE_FIBER_DATASTRUCTS(file_a, file_b)
%   report = COMPARE_FIBER_DATASTRUCTS(..., 'Tolerance', 1e-9, 'MaxReport', 40)
%
%   Inputs
%   ------
%   a, b : either two structs, OR two paths to .mat files. When a path is given,
%          the WHOLE file is loaded (load(file) -> struct of all variables) so it
%          works whatever the saved variable is called.
%
%   Name-value options
%   ------------------
%   'Tolerance' : max allowed abs difference for numeric arrays (default 0 =
%                 bit-identical). Use e.g. 1e-9 to allow tiny rounding.
%   'MaxReport' : max number of mismatch lines to print (default 40).
%   'Verbose'   : true to also print each field as it passes (default false).
%
%   Output
%   ------
%   report : struct with fields
%       .pass        - true if no mismatches
%       .n_compared  - number of leaf values compared
%       .n_mismatch  - number of mismatches
%       .messages    - cellstr of human-readable mismatch descriptions
%
%   Typical use (see VALIDATION.md): save a reference output BEFORE a refactor,
%   re-run AFTER, then
%       r = compare_fiber_datastructs('before.mat', 'after.mat');
%   A green ALL MATCH means the refactor preserved the output exactly.

p = inputParser;
p.addParameter('Tolerance', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('MaxReport', 40, @(x) isnumeric(x) && isscalar(x));
p.addParameter('Verbose', false, @(x) islogical(x) || ismember(x, [0 1]));
p.parse(varargin{:});
opt = p.Results;

a = load_if_path(a);
b = load_if_path(b);

state.tol = opt.Tolerance;
state.verbose = logical(opt.Verbose);
state.n_compared = 0;
state.messages = {};

state = compare_node('', a, b, state);

report.pass       = isempty(state.messages);
report.n_compared = state.n_compared;
report.n_mismatch = numel(state.messages);
report.messages   = state.messages;

% --- Print a friendly summary --------------------------------------------
fprintf('\n=== compare_fiber_datastructs ===\n');
fprintf('  leaf values compared : %d\n', report.n_compared);
fprintf('  mismatches           : %d  (tolerance = %g)\n', report.n_mismatch, opt.Tolerance);
nshow = min(report.n_mismatch, opt.MaxReport);
for i = 1:nshow
    fprintf('    - %s\n', report.messages{i});
end
if report.n_mismatch > nshow
    fprintf('    ... and %d more (raise MaxReport to see them)\n', report.n_mismatch - nshow);
end
if report.pass
    fprintf('  RESULT: ALL MATCH  ✓\n\n');
else
    fprintf('  RESULT: MISMATCH  ✗\n\n');
end
end

% =========================================================================
function s = load_if_path(s)
if ischar(s) || (isstring(s) && isscalar(s))
    s = load(char(s));   % whole-file load -> struct of all saved variables
end
end

% =========================================================================
function state = compare_node(path, a, b, state)
% Recursively walk two values in lock-step.

% --- type mismatch -------------------------------------------------------
if ~strcmp(class(a), class(b))
    state = add(state, path, sprintf('class differs: %s vs %s', class(a), class(b)));
    return;
end

if isstruct(a)
    state = compare_struct(path, a, b, state);
elseif iscell(a)
    state = compare_cell(path, a, b, state);
else
    state = compare_leaf(path, a, b, state);
end
end

% =========================================================================
function state = compare_struct(path, a, b, state)
if ~isequal(size(a), size(b))
    state = add(state, path, sprintf('struct size differs: %s vs %s', sz(a), sz(b)));
    return;
end
fa = fieldnames(a); fb = fieldnames(b);
missing = setdiff(fa, fb);
extra   = setdiff(fb, fa);
for i = 1:numel(missing), state = add(state, [path '.' missing{i}], 'field only in A'); end
for i = 1:numel(extra),   state = add(state, [path '.' extra{i}],   'field only in B'); end
common = intersect(fa, fb);
% Handle struct arrays element-by-element
for e = 1:numel(a)
    if numel(a) > 1, epath = sprintf('%s(%d)', path, e); else, epath = path; end
    for i = 1:numel(common)
        f = common{i};
        state = compare_node([epath '.' f], a(e).(f), b(e).(f), state);
    end
end
end

% =========================================================================
function state = compare_cell(path, a, b, state)
if ~isequal(size(a), size(b))
    state = add(state, path, sprintf('cell size differs: %s vs %s', sz(a), sz(b)));
    return;
end
for i = 1:numel(a)
    state = compare_node(sprintf('%s{%d}', path, i), a{i}, b{i}, state);
end
end

% =========================================================================
function state = compare_leaf(path, a, b, state)
state.n_compared = state.n_compared + 1;

if ischar(a)
    if ~strcmp(a, b)
        state = add(state, path, sprintf('char differs: "%s" vs "%s"', a, b));
    elseif state.verbose, fprintf('    ok  %s\n', path);
    end
    return;
end

if ~isequal(size(a), size(b))
    state = add(state, path, sprintf('size differs: %s vs %s', sz(a), sz(b)));
    return;
end

if isnumeric(a) || islogical(a)
    da = double(a); db = double(b);
    % NaN positions must match
    if ~isequal(isnan(da), isnan(db))
        state = add(state, path, 'NaN pattern differs');
        return;
    end
    mask = ~isnan(da);
    if any(mask(:))
        d = max(abs(da(mask) - db(mask)));
    else
        d = 0;
    end
    if d > state.tol
        state = add(state, path, sprintf('numeric differs: max|diff|=%.3g', d));
    elseif state.verbose, fprintf('    ok  %s (max|diff|=%.3g)\n', path, d);
    end
    return;
end

% Fallback for any other type
if ~isequaln(a, b)
    state = add(state, path, sprintf('values differ (class %s)', class(a)));
elseif state.verbose, fprintf('    ok  %s\n', path);
end
end

% =========================================================================
function state = add(state, path, msg)
if isempty(path), path = '<root>'; end
state.messages{end+1} = sprintf('%s : %s', path, msg);
end

function s = sz(x)
s = ['[' strjoin(arrayfun(@(d) num2str(d), size(x), 'uni', 0), 'x') ']'];
end
