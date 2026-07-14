"""
Quick stats-only script: print unsplit theta cross-correlation statistics
for a given animal from the multisite analysis without generating figures.

USAGE:  python print_xcorr_stats.py [animal_name]
        python print_xcorr_stats.py              # defaults to Animal01
"""

import sys
from pathlib import Path
import numpy as np

# multisite_fiber_analysis.py (and its RECORDINGS cohort config) lives in the
# Fig6_bilateral_ca1/ figure folder, not here -- point sys.path at it.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Fig6_bilateral_ca1"))
from multisite_fiber_analysis import (
    collect_trialwise_xcorr_unsplit,
    peak_lag_index_restricted,
    RECORDINGS,
    THETA_BAND,
    XCORR_UNSPLIT_N_SURROGATES,
)

ANIMAL = sys.argv[1] if len(sys.argv) > 1 else "Animal01"

if ANIMAL not in RECORDINGS:
    print(f"ERROR: '{ANIMAL}' not in RECORDINGS. "
          f"Available: {list(RECORDINGS.keys())}")
    sys.exit(1)

sessions = RECORDINGS[ANIMAL]
n_sessions = len(sessions)
n_trials_total = sum(v["n_trials"] for v in sessions.values())

print("=" * 70)
print(f"UNSPLIT THETA CROSS-CORRELATION STATS: {ANIMAL}")
print(f"Sessions: {n_sessions} ({', '.join(sessions.keys())}), "
      f"max trials: {n_trials_total}")
print(f"Band: {THETA_BAND[0]}-{THETA_BAND[1]} Hz, "
      f"Surrogates: {XCORR_UNSPLIT_N_SURROGATES}/trial")
print("=" * 70)

result = collect_trialwise_xcorr_unsplit(ANIMAL, session=None, band=THETA_BAND)

if result is None:
    print("ERROR: No data returned.")
    sys.exit(1)

n_loaded = result["n_loaded"]
n_sig = result["n_sig_trials"]
pooled_p = result["pooled_p"]
obs = result["observed_peaks"]
median_r = float(np.median(obs))
mean_r = float(np.mean(obs))
mean_lag = float(np.mean([t["peak_lag"] for t in result["trials"]]))
null_p95 = float(np.percentile(result["pooled_null"], 95))

ir_gm, _ = peak_lag_index_restricted(
    result["grand_mean"], result["lags_ms"], THETA_BAND)
gm_peak_r = float(result["grand_mean"][ir_gm])
gm_peak_lag = float(result["lags_ms"][ir_gm])

p_str = (f"p = {pooled_p:.4f}" if np.isfinite(pooled_p) and pooled_p >= 0.001
         else f"p = {pooled_p:.2e}" if np.isfinite(pooled_p) else "p = n/a")

print()
print("=" * 70)
print(f"SURROGATE TEST RESULTS — {ANIMAL}")
print("=" * 70)
print(f"  Sessions:              {n_sessions}")
print(f"  Trials loaded:         {n_loaded}")
print(f"  Band:                  {THETA_BAND[0]}-{THETA_BAND[1]} Hz (theta)")
print(f"  Surrogates/trial:      {XCORR_UNSPLIT_N_SURROGATES}")
print(f"  Grand-mean peak r:     {gm_peak_r:.4f} at {gm_peak_lag:.1f} ms")
print(f"  Median per-trial |r|:  {median_r:.4f}")
print(f"  Mean per-trial |r|:    {mean_r:.4f}")
print(f"  Mean peak lag:         {mean_lag:.1f} ms")
print(f"  Trials significant:    {n_sig}/{n_loaded} (p<0.05, per-trial)")
print(f"  Pooled {p_str}")
print(f"  Null 95th pctl:        {null_p95:.4f}")

per_trial_lines = []
for t in result["trials"]:
    per_trial_lines.append(
        f"    {t['session']} T{t['trial_num']}: "
        f"|r|={t['peak_r']:.4f}, lag={t['peak_lag']:.1f} ms, "
        f"p={t['p_surrogate']:.4f}")
print(f"\n  Per-trial breakdown:")
for line in per_trial_lines:
    print(line)

print("-" * 70)
print(f"RESULT: Fiber-fiber theta-band cross-correlation "
      f"({ANIMAL}, n={n_loaded} trials across {n_sessions} sessions): "
      f"grand-mean peak r = {gm_peak_r:.4f} at {gm_peak_lag:.1f} ms, "
      f"median per-trial peak |r| = {median_r:.4f}, "
      f"mean peak lag = {mean_lag:.1f} ms, "
      f"{n_sig}/{n_loaded} trials individually significant, "
      f"pooled {p_str} "
      f"(circular-shift surrogate test, "
      f"{XCORR_UNSPLIT_N_SURROGATES} surrogates/trial)")
print("=" * 70)
