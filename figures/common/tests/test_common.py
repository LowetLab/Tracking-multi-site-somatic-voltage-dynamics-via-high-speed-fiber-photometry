"""Unit tests for figures/common/common.py shared helpers (no data needed).

Run:  python figures/common/tests/test_common.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import common  # noqa: E402


def _check(name, got, want, fails):
    ok = got == want
    print(f"  {'PASS' if ok else 'FAIL'} {name}" + ("" if ok else f"  got={got!r} want={want!r}"))
    return fails + (0 if ok else 1)


def test_to_long_path():
    f = 0
    # UNC network share -> \\?\UNC\...
    f = _check("UNC share", common.to_long_path(r"\\store\dept\file.mat"),
               r"\\?\UNC\store\dept\file.mat", f)
    # Local drive path -> \\?\C:\...
    f = _check("local drive", common.to_long_path(r"C:\data\file.mat"),
               r"\\?\C:\data\file.mat", f)
    # Already extended -> unchanged
    f = _check("already extended UNC", common.to_long_path(r"\\?\UNC\x\y"),
               r"\\?\UNC\x\y", f)
    f = _check("already extended drive", common.to_long_path(r"\\?\C:\x"),
               r"\\?\C:\x", f)
    # Relative path -> unchanged
    f = _check("relative path", common.to_long_path(r"folder\file.mat"),
               r"folder\file.mat", f)
    # Accepts non-str (e.g. pathlib) via str()
    f = _check("non-str coerced", common.to_long_path(123), "123", f)
    return f


def test_next_pow2():
    f = 0
    f = _check("_next_pow2(1)", common._next_pow2(1), 1, f)
    f = _check("_next_pow2(5)", common._next_pow2(5), 8, f)
    f = _check("_next_pow2(16)", common._next_pow2(16), 16, f)
    f = _check("_next_pow2(17)", common._next_pow2(17), 32, f)
    return f


def test_infer_trial_from_name():
    f = 0
    f = _check("trial match", common._infer_trial_from_name("M_Trial7_FiberPhotometry_Analysis.mat"), 7, f)
    f = _check("trial match 2-digit", common._infer_trial_from_name("x_Trial12_FiberPhotometry_Analysis.mat"), 12, f)
    f = _check("case-insensitive", common._infer_trial_from_name("x_trial3_fiberphotometry_analysis.MAT"), 3, f)
    f = _check("no match -> None", common._infer_trial_from_name("Summary.txt"), None, f)
    return f


def test_cmaps():
    f = 0
    om = common.create_monochromatic_orange_cmap()
    pa = common.create_parula_like_cmap()
    f = _check("orange cmap N", om.N, 256, f)
    f = _check("orange cmap name", om.name, "mono_orange", f)
    f = _check("parula cmap N", pa.N, 256, f)
    f = _check("parula cmap name", pa.name, "parula_like", f)
    return f


if __name__ == "__main__":
    print("=== test_common ===")
    fails = test_to_long_path()
    fails += test_next_pow2()
    fails += test_infer_trial_from_name()
    fails += test_cmaps()
    if fails == 0:
        print("\nALL test_common TESTS PASSED")
        sys.exit(0)
    print(f"\n{fails} test(s) FAILED")
    sys.exit(1)
