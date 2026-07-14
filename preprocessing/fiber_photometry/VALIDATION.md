# Validating the preprocessing pipeline

Because the full pipeline can't run automatically (it needs raw imaging + Open Ephys data +
interactive ROI selection), use the two checks below to confirm changes to `core/` or the entry
scripts haven't broken anything. Tooling lives in `preprocessing/fiber_photometry/tests/`.

There are **two levels**. Always do Level 1 after any change. Do Level 2 as well if the change
touches signal processing, output-struct assembly, or anything that could affect the saved
`.mat` contents.

---

## Level 1 — Offline check (~30 seconds, no data needed)

This confirms the shared helpers are correct and resolve to `core/utils` (not an
external or built-in version), and that the comparison tool itself works.

```matlab
cd 'C:\path\to\your_project_root\preprocessing\fiber_photometry'
addpath('tests')
validate_refactor
```

**Expect:** a block ending in `OFFLINE VALIDATION: PASS`.

---

## Level 2 — Real-recording check (~5–10 min)

The idea: produce the output **before** and **after** the change for the **same
recording with the same ROIs**, then let the tool compare them automatically.

### 2a. Make a "before" reference (do this ONCE, on the current code)

The cleanest reference is an output file you already trust. Either:

- **Use an existing processed file** you produced earlier (an
  `*_FiberPhotometry_Analysis.mat`), **or**
- Make a fresh one now: open
  `run_fiber_preprocessing_singletrial.m`, set the config
  at the top to a recording, run it, and note where it saved the
  `*_FiberPhotometry_Analysis.mat`.

Copy that file somewhere stable and rename it, e.g. `before.mat`.

### 2b. Re-run AFTER the stage, same recording, same ROIs

Pull the new commit, run the **same** script on the **same** recording, and
**select the same ROIs in the same order** when prompted. Note the new output
file; call it `after.mat`.

### 2c. Compare automatically

```matlab
addpath('tests')
r = compare_fiber_datastructs('before.mat', 'after.mat');
```

**Expect:** `RESULT: ALL MATCH ✓` (0 mismatches).

---

## Reading the result

`compare_fiber_datastructs` prints every difference with its exact location, e.g.

```
.trials{2}.signals.fluorescence_corrected : numeric differs: max|diff|=3.1e-04
```

- **0 mismatches** → the refactor preserved the output exactly. 
- **Only ROI-derived fields differ by tiny amounts** → almost always means your
  hand-drawn ROIs were a few pixels different between runs, not a code bug. Re-run
  2b being careful to match ROIs, or allow a small tolerance to confirm:
  ```matlab
  r = compare_fiber_datastructs('before.mat', 'after.mat', 'Tolerance', 1e-6);
  ```
- **Whole fields missing / structural differences / large diffs** → a real
  regression. Send me the printed list.

---

## If something fails

Collect:
1. the full `validate_refactor` output (Level 1), or the
   `compare_fiber_datastructs` mismatch list (Level 2), and
2. any red MATLAB error text.

Compare against the last known-good commit to narrow down which change introduced the issue.
