# MATLAB Preprocessing — Workflow Schematics

Visual + textual walkthrough of the two MATLAB preprocessing pipelines, traced
directly from the source scripts. Line references point at the canonical entry
scripts so you can jump to the exact code.

- **Fiber photometry + Open Ephys** →
  [`preprocessing/fiber_photometry/run_fiber_preprocessing_singletrial.m`](preprocessing/fiber_photometry/run_fiber_preprocessing_singletrial.m)
  (the `*_multitrial_*` script is the same flow wrapped in a per-trial loop).
  Config: [`config/fiber_preprocessing_singletrial_config.m`](preprocessing/fiber_photometry/config/fiber_preprocessing_singletrial_config.m).
- **Single-cell (cellular) voltage imaging** →
  [`preprocessing/cellular_imaging/cellular_processing_multitrial_DBS.m`](preprocessing/cellular_imaging/cellular_processing_multitrial_DBS.m).
  Batch wrapper: [`run_cellular_batch_processing.m`](preprocessing/cellular_imaging/run_cellular_batch_processing.m).

Both share the same recipe — *imaging → ROI traces → align to Open Ephys via
camera triggers → photobleaching/ΔF/F → spectral/spike analysis → `.mat`
datastruct* — but differ in ROI granularity (fiber bundles vs single neurons),
the downstream analysis (LFP phase-locking vs spike detection), and trial
handling.

---

## 1 · Fiber Photometry + Open Ephys Preprocessing

### 1.1 Inputs

```
RAW DATA (external, under  …\Data\FiberVoltageImaging\<Mouse>\…)
┌────────────────────────────────┐     ┌──────────────────────────────────────────┐
│ IMAGING                        │     │ OPEN EPHYS  (Record Node 103/104, 30 kHz)  │
│ <…>\ImagingData\<Exp>\<Date>\  │     │ <…>\OpenEphys\<Exp>\<Date>\…<R##>…\         │
│   <R##>\<trial>\ *.ome.tif      │     │   100_RhythmData_Ch11.continuous  → HP LFP │
│   (multi-part OME-TIFF stacks,  │     │   100_RhythmData_Ch1 / Ch3        → mPFC ± │
│    GEVI voltage imaging through │     │   100_RhythmData_ADC1   → camera triggers  │
│    optic fibers; HP + mPFC)     │     │   100_RhythmData_ADC7   → trial markers    │
│                                 │     │   100_RhythmData_ADC5/6 → stim pulse trains│
│                                 │     │   100_RhythmData_ADC4   → running wheel     │
└────────────────────────────────┘     └──────────────────────────────────────────┘
```
Paths are auto-constructed from `MOUSE_NAME / RECORDING_DATE / RECORDING_ID`
(config §"BASE PATH"). The recording folder is matched by a digit-safe regex so
`R1` never matches `R10`; the Record Node is probed 103 → 104; on any miss the
script falls back to a `uigetfile` prompt.

### 1.2 The pipeline (10 stages)

```
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │ STAGE 0 · CONFIG & SETUP                              §1, lines 16–113         │
 │  run() config include  →  all UPPER_CASE params into workspace                 │
 │  setup_lab_paths() (FieldTrip, NoRMCorre, load_open_ephys_data, core/)         │
 │  ANALYSIS_MODE = 'single_trial' | 'multi_fov'  → auto-find trial folder        │
 └──────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │ STAGE 1 · LOAD IMAGING                               §2, lines 115–325         │
 │  • find & name-sort *.ome.tif parts → concatenate into one uint16 stack        │
 │  • [optional] MOTION CORRECTION (NoRMCorre, rigid):                            │
 │      single-trial → draw rect ROI, normcorre on window                        │
 │      multi-FOV    → boxfilt3 + DoG high-pass → estimate shifts → circshift     │
 │  • permute [2 1 3] to fix orientation                                          │
 └──────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │ STAGE 2 · EXTRACT FIBER TRACES                       §3, lines 327–388         │
 │  PROCESS_FULL_FIELD ?                                                          │
 │    true  → single whole-field mean trace                                       │
 │    false → INTERACTIVE: draw N polygon ROIs (one per fiber) on a              │
 │            clean_display_frame (99.9-pct clip); trace = mean over ROI pixels   │
 │  → all_traces  [frames × num_fibers]                                           │
 └──────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │ STAGE 3 · ARTIFACT REMOVAL (per trace)              §4, lines 390–433         │
 │  • clamp last ~80 samples to a median (end artifact)                           │
 │  • notch-out 120–124 Hz and 130–132 Hz (Butter bandpass + subtract)           │
 │    → removes stim/line harmonics from the optical signal                       │
 │  (uses provisional IMAGING_FS=500; final Fs set in Stage 4)                    │
 └──────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │ STAGE 4 · LOAD OPEN EPHYS + ALIGN  (the heart)      §5, lines 436–917         │
 │  Load Ch11 (HP LFP), Ch1−Ch3 (mPFC differential), ADC1/4/5/6/7                 │
 │                                                                               │
 │  ① STIM CHANNEL ID  — count pulses on ADC5 vs ADC6 → pulse channel            │
 │  ② BASELINE vs STIM TRIAL  — IS_BASELINE_TRIAL=[] auto: <5 pulses ⇒ baseline   │
 │  ③ STIM ONSET  — detect_stim_onset() on pulse train (first onset in trial),    │
 │       offset = onset + fixed STIMULATION_DURATION (10 s)  → STIM_DETECTION     │
 │  ④ TRIAL WINDOW — ADC7 > 0.5 → trial_starts/stops                             │
 │  ⑤ CAMERA TRIGGERS — ADC1 rising edges, min-spacing>8 samp, keep in-trial     │
 │  ⑥ IMAGING_FS  = median(1/Δtrigger)   ← true frame rate from hardware          │
 │  ⑦ ALIGN: clip traces & triggers to min length; LFP_aligned = lfp(trigIdx)     │
 │       (≈30 kHz → ~frame rate downsample by indexing at trigger samples)        │
 │  ⑧ CLEAN: replace_outliers_with_median(LFP,10); smooth running velocity        │
 │  ⑨ map STIM samples → imaging frames → STIM_PERIOD (seconds)                   │
 └──────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │ STAGE 5 · PHOTOBLEACHING + ΔF/F                     §6/6B, lines 935–1067      │
 │  correct_photobleaching():  linear detrend  +  double-exponential fit          │
 │     fit on PRE-stim baseline window only  → traces_detrended,                  │
 │                                              traces_exp_corrected              │
 │  compute_deltaF_F(): F0 = mean over baseline window                            │
 │     stim trial → [−60 s … stim onset];   baseline trial → whole 60 s trace     │
 │     → traces_deltaF_F  +  z-scored version                                     │
 │  processed_traces := ΔF/F   (canonical signal downstream)                      │
 └──────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │ STAGE 6 · SPECTRAL / TIME-FREQUENCY  (per fiber)    §7+§8, lines 1090–1850     │
 │  • FieldTrip ft_freqanalysis (mtmconvol, hanning, foi 2–70 Hz, 0.5 s win)     │
 │      on [Fiber, LFP-HP, (mPFC)] → complex spectra (Fourier) → PHASE per chan   │
 │  • MATLAB spectrogram() (480/449 win, 1–150 Hz) per signal                     │
 │  • compute_band_power() over BAND_RANGES (δθ/α/β/lowγ/highγ)                    │
 │      → zscore_smooth_bands()                                                    │
 └──────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │ STAGE 7 · COUPLING ANALYSIS  (per fiber)            §8, lines 2560–2750        │
 │  • PLV(f) = |mean_t exp(i·Δφ)|,  Δφ = circ_dist(phase_fiber, phase_LFP)        │
 │      computed overall AND split by RUN vs REST (running_velocity threshold)     │
 │  • Theta (5–10 Hz) envelope: |hilbert(LFP θ)| vs fiber ΔF/F → corr (r,p)        │
 │  • HP↔mPFC LFP phase difference (when mPFC loaded)                              │
 └──────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │ STAGE 8 · FIGURES                                   §8, lines 1281–2756        │
 │  per-fiber: raw-vs-corrected traces, spectrograms (fiber/LFP/mPFC),            │
 │  band-power stacks, PSD, photobleach-fit quality, PLV(f), run-vs-rest PLV,      │
 │  θ-envelope correlation, multi-fiber overlay, stim-detection diagnostic        │
 │  → saved as .fig + .png in the session output folder                           │
 └──────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │ STAGE 9 · SAVE DATASTRUCT                           §9/§10, lines 2757–3160    │
 │  FiberPhotometryAnalysis struct (saved -v7.3) →                               │
 │   …\<Mouse>\Fiber_Voltage_Processed\<Date>-<R##>\                              │
 └──────────────────────────────────────────────────────────────────────────────┘
```

### 1.3 Output datastruct — `FiberPhotometryAnalysis` (frozen schema)

```
FiberPhotometryAnalysis
├── metadata        mouse / date / recording_id / experimenter
├── time            time_vector_seconds, sampling_rate (= IMAGING_FS)
├── signals         raw_traces → filtered_traces → corrected/exp → deltaF_F
│                   → zscored → final_processed; F0_values, baseline window
├── ephys           lfp_HP_full_original, lfp_raw_aligned_HP, lfp_z_HP,
│                   lfp_mPFC_* , running_velocity(_smooth), stim_pulses,
│                   stim_onset, camera_trigger_indices, theta_env_correlation
├── photobleaching  methods_used, fit_params, fitted_curve, per-fiber rate
├── rois            type (manual/full_field) + ROI masks
└── spectral        band_names/ranges/colors + per-fiber spectrogram & bands
```

> **Sidecar step — LFP artifact masks.**
> [`artifact_removal_lfp.m`](preprocessing/fiber_photometry/artifact_removal_lfp.m)
> / `_multisession.m` run separately to detect stim/movement artifacts in the
> LFP and write `*_artifact_removal.mat` masks. Downstream
> `spectral_analysis/` consumes these in its `exclude` / `clean` modes.

### 1.4 Key parameters (edit in the config file, not the script)

| Param | Default | Meaning |
|---|---|---|
| `EPHYS_FS` | 30000 | Open Ephys hardware sample rate |
| `IMAGING_FS` | `[]` → measured | Frame rate from camera-trigger intervals |
| `MOTION_CORRECTION` | 0 | NoRMCorre rigid motion correction on/off |
| `APPLY_PHOTOBLEACHING_CORRECTION` | true | linear + double-exp detrend |
| `LOAD_mPFC_LFP` | true | also build mPFC differential (Ch1−Ch3) |
| `IS_BASELINE_TRIAL` | `[]` | auto (<5 stim pulses ⇒ baseline) / force |
| `BAND_RANGES` | `[1 8;8 12;13 30;31 70;71 185]` | δθ / α / β / lowγ / highγ |
| `STIMULATION_DURATION` | 10 s | fixed stim offset = onset + duration |

---

## 2 · Single-Cell (Cellular) Voltage Imaging Pipeline

Same backbone, but ROIs are **individual neurons**, ROIs are **shared across
trials** of a session, the headline analysis is **spike detection** (not LFP
phase-locking), and it loops over **multiple trials** per session. Built for DBS
(40 Hz vs 130 Hz, amplitude- vs energy-balanced).

### 2.1 Inputs

```
RAW DATA (…\DBSData\<Mouse>\<Exp>\Voltage_Imaging\<Date>\<R##>\<trial_1,2,3…>\*.tif)
┌────────────────────────────────┐     ┌──────────────────────────────────────────┐
│ IMAGING (per trial)            │     │ OPEN EPHYS (30 kHz)                         │
│   *.tif (multi-file, name-sort │     │   Ch11           → LFP                      │
│    & concatenated per trial)   │     │   ADC1           → camera triggers         │
│   CamKII GEVI, single neurons  │     │   ADC7           → trial illumination mark │
│                                 │     │   ADC6           → stimulus pulse train    │
│                                 │     │   ADC2/3/4       → motion X/Y/Z (defined)  │
└────────────────────────────────┘     └──────────────────────────────────────────┘
Session metadata (freq, #trials, #neurons, comparison) in ANIMAL_SESSIONS table.
```

### 2.2 The pipeline

```
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │ STAGE 0 · CONFIG                                     §1, lines 34–143          │
 │  Mouse/date/rec id, INDICATOR (CamKII, polarity), DBS params, CHANNELS map,    │
 │  SPIKE_PARAMS, session DB → build DATA_FOLDER + OUTPUT_FOLDER                  │
 └──────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │ STAGE 1 · DETECT TRIALS + EPHYS PATH                §2/§3, lines 145–281       │
 │  auto-detect trial_1/_2/_3… folders;  build Open Ephys Record-Node path        │
 └──────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
        ╔═══════════════ FOR EACH TRIAL  (§4, lines 282–1243) ═══════════════╗
        ║                                                                     ║
        ║  A. LOAD IMAGING        concat all *.tif → uint16 stack             ║
        ║                                                                     ║
        ║  B. MOTION ROI  ─ trial 1 only ─ draw rect (imrect)  ──┐            ║
        ║     (shared_motion_roi_position reused for all trials) │            ║
        ║  C. MOTION CORRECT      NoRMCorre rigid (boxfilt3 +    │  SHARED    ║
        ║     DoG high-pass → shifts → circshift) → permute      │  across    ║
        ║                                                        │  trials    ║
        ║  D. NEURON ROIs ─ trial 1 only ─ draw N polygons      ─┘            ║
        ║     → shared_neuron_roi_masks + centroids (reused → same neurons    ║
        ║       tracked across every trial of the session)                    ║
        ║                                                                     ║
        ║  E. EXTRACT TRACES      per-neuron mean over ROI pixels             ║
        ║                         + background trace (all non-ROI pixels)     ║
        ║                                                                     ║
        ║  F. LOAD + ALIGN EPHYS                                              ║
        ║     • camera triggers (ADC1, edges, spacing>8) → IMAGING_FS=median  ║
        ║       (fallback: 650 Hz + illumination-window resampling)           ║
        ║     • trial window from ADC7 illumination on/off                    ║
        ║     • stim onset from ADC6 pulse train inside trial window,         ║
        ║       offset = onset + DBS_DURATION_SEC → frame indices             ║
        ║     • LFP_aligned = lfp(trigger samples);  time_vector zeroed       ║
        ║       at stimulus onset                                             ║
        ║                                                                     ║
        ║  G. PHOTOBLEACHING (per neuron, fit on PRE-stim only)              ║
        ║     • linear detrend → normalize to baseline                        ║
        ║     • double-exponential A1e^-t/τ1+A2e^-t/τ2+C (lsqcurvefit),        ║
        ║       corrected = trace./fit → ΔF/F = (corr−F0)/F0; z-scored too    ║
        ║     • on fit failure → fall back to linear-detrended trace          ║
        ║                                                                     ║
        ║  H. SPIKE DETECTION (per neuron)        SPIKE_DETECTION_TYPE        ║
        ║     • 'baselineSD' → spike_detect_SNR_sim3_SC                       ║
        ║     • 'thresholdSD' → spike_detect_baseline_threshold_SC           ║
        ║       baseline window = pre-stim period (else first 30%)            ║
        ║     → rasters, SNR, denoised/highpass/subthreshold traces,          ║
        ║       firing rate, amplitude, noise level                           ║
        ║                                                                     ║
        ║  I. STORE trial_result struct (metadata, params, time, signals,     ║
        ║     spikes, photobleaching) → all_trial_results{trial_idx}          ║
        ╚═════════════════════════════════════════════════════════════════════╝
                                      │
                                      ▼
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │ STAGE 2 · SAVE + FIGURES                            §5/§6, lines 1244–end      │
 │  CellularAnalysis = {metadata, trials[], shared_rois}                          │
 │  → <Mouse>\<Exp>\CellularDataProcessed\<Date>-<R##>\                           │
 │       <Mouse>_<Date>-<R##>_CellularAnalysis.mat   (-v7.3)                      │
 │  + per-trial / per-neuron figures (traces, rasters, ROI map)                   │
 └──────────────────────────────────────────────────────────────────────────────┘
```

### 2.3 Output datastruct — `CellularAnalysis`

```
CellularAnalysis
├── metadata          mouse/date/rec, indicator, DBS freq + comparison,
│                     is_baseline, num_trials, num_neurons
├── shared_rois       masks{}, centroids   ← SAME neurons across all trials
└── trials{ }         one entry per trial:
    ├── metadata      trial number / folder / index
    ├── parameters    imaging_fs, ephys_fs, num_frames/neurons, duration,
    │                 motion+photobleach flags, spike params
    ├── time(_periods) time_vector (0 = stim onset), stim onset/offset frames
    ├── signals       raw / detrended / exp-corrected (ΔF/F) / z-scored,
    │                 F0_values, denoised / highpass / subthreshold, background
    ├── spikes        rasters, firing_rate, SNR, amplitude, noise level
    └── ephys         LFP aligned, stimulus aligned
```

### 2.4 Fiber vs Cellular — the differences at a glance

| | Fiber photometry | Single-cell |
|---|---|---|
| ROI unit | fiber bundle / whole field | individual neuron (polygon) |
| ROIs across trials | per-trial (single-trial script) | **shared** (drawn on trial 1) |
| Trial handling | one trial / FOV loop | **multi-trial loop**, one session |
| Stim channel | ADC5 **or** ADC6 (auto-pick) | ADC6 (fixed) |
| 2nd LFP site | mPFC differential (Ch1−Ch3) | — |
| Headline analysis | LFP–fiber **PLV**, θ-envelope corr | **spike detection** (SNR/SD) |
| Spectral | FieldTrip TFR + band power | (not in preprocessing) |
| Output | `FiberPhotometryAnalysis` | `CellularAnalysis` (+ `trials{}`) |
| Output root | `Fiber_Voltage_Processed\` | `CellularDataProcessed\` |

---

### Shared conventions (both pipelines)
- **Alignment is trigger-driven:** the true imaging sample rate is *measured*
  from camera-trigger spacing on ADC1, and every ephys signal is brought to the
  imaging timebase by indexing at the trigger samples (`lfp(trigger_indices)`).
- **Stim detection is data-driven** from the pulse train, with a fixed-duration
  offset; baseline trials are auto-recognised by pulse count.
- **Photobleaching** = linear detrend + double-exponential fit on the pre-stim
  baseline; ΔF/F uses a static baseline window.
- **Outputs** are `-v7.3` `.mat` structs (HDF5) → read by the Python/spectral
  downstream stages via `h5py` / `scipy.io`.
