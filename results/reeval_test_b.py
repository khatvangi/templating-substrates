"""
re-evaluate test b criterion 3 with corrected periodicity logic.

original criterion 3 demanded autocorrelation_peak_lag == N for every cell,
which fails because Y is also periodic at 2N, 3N, ...; sampling noise
decides which multiple wins the argmax. corrected logic:
  3a) peak_lag is a multiple of N
  3b) the autocorrelation value at the peak (which equals the value at lag N
      up to sampling noise, when peak_lag is a multiple of N) matches the
      predicted match prob (1-eps)^2 + eps^2/(N-1)
cells with L == N (only one cycle) are skipped: lag-N autocorr has no pairs.
"""

import csv
from pathlib import Path

CSV_PATH = Path("/storage/kiran-stuff/templating_framework/results/test_b_results.csv")
OUT_PATH = Path("/storage/kiran-stuff/templating_framework/results/test_b_verdict_corrected.md")
TOL = 0.05

rows = []
with CSV_PATH.open() as f:
    reader = csv.DictReader(f)
    for r in reader:
        rows.append({
            "N": int(r["N"]),
            "epsilon": float(r["epsilon"]),
            "L": int(r["L"]),
            "peak_lag": int(r["autocorrelation_peak_lag"]),
            "peak_val": float(r["autocorrelation_peak_value"]),
            "I_emp": float(r["I_empirical"]),
            "log2_N": float(r["log2_N_theoretical"]),
        })


def predicted_match_prob(N: int, eps: float) -> float:
    return (1 - eps) ** 2 + eps ** 2 / (N - 1)


# evaluate criterion 3 per row
n_skipped = 0
n_3a_pass = 0
n_3b_pass = 0
n_combined_pass = 0
n_eligible = 0
failures = []

for r in rows:
    N, L, eps = r["N"], r["L"], r["epsilon"]
    pred = predicted_match_prob(N, eps)
    r["predicted_match_prob"] = pred

    if L <= N:
        # L == N: only one cycle, lag-N autocorr undefined. L < N: should not occur.
        r["status"] = "L=N: skipped" if L == N else "L<N: skipped (defensive)"
        r["test_3a"] = None
        r["test_3b"] = None
        n_skipped += 1
        continue

    n_eligible += 1
    test_3a = (r["peak_lag"] % N == 0)
    test_3b = abs(r["peak_val"] - pred) < TOL
    r["test_3a"] = test_3a
    r["test_3b"] = test_3b
    if test_3a:
        n_3a_pass += 1
    if test_3b:
        n_3b_pass += 1
    if test_3a and test_3b:
        n_combined_pass += 1
        r["status"] = "PASS"
    else:
        r["status"] = "FAIL"
        failures.append(r)

# criterion 1 and 2 are taken as established PASS from prior analysis
# (per the task: only criterion 3 is being re-evaluated here).
crit1_pass = True
crit2_pass = True
crit3_pass = (n_combined_pass == n_eligible)
overall = crit1_pass and crit2_pass and crit3_pass

# build markdown verdict
lines = []
lines.append("# Test B verdict (corrected criterion 3)")
lines.append("")
lines.append("## Why this re-evaluation exists")
lines.append("")
lines.append(
    "The original criterion 3 required `autocorrelation_peak_lag == N` for "
    "every cell. This is mathematically wrong: a sequence Y with period N is "
    "also periodic with period 2N, 3N, ..., so the autocorrelation function "
    "has equal-magnitude peaks at every integer multiple of N. Sampling noise "
    "determines which multiple wins the argmax. The framework's actual "
    "prediction is *Y is periodic with period N* (i.e., y_{i+N} = y_i in the "
    "noiseless case)."
)
lines.append("")
lines.append("Corrected operationalization:")
lines.append("")
lines.append("- **Test 3a:** `peak_lag % N == 0` (peak at some integer multiple of N)")
lines.append(
    "- **Test 3b:** `|peak_val - predicted_match_prob| < 0.05`, with "
    "`predicted_match_prob = (1-eps)^2 + eps^2/(N-1)`. When peak_lag is a "
    "multiple of N, peak_val equals the lag-N autocorrelation up to sampling "
    "noise, so the peak_val column is a valid stand-in."
)
lines.append("- **Skip:** cells with `L == N` (only one polymer cycle, lag-N autocorr has zero pairs).")
lines.append("")
lines.append("## Summary of all three criteria")
lines.append("")
lines.append("- **Criterion 1** (low-noise saturation, I_empirical -> log_2(N) at eps=0.001): "
             "PASS (established by prior analysis; every cell within 0.00-0.01% relative error)")
lines.append("- **Criterion 2** (no L-growth of saturated information): "
             "PASS (established by prior analysis; every (N, eps) group has max/min ratio < 1.0001)")
lines.append(f"- **Criterion 3 (corrected)** (peak at multiple of N AND peak_val matches "
             f"predicted match prob within {TOL}): "
             f"{'PASS' if crit3_pass else 'FAIL'}")
lines.append("")
lines.append("## Criterion 3 breakdown")
lines.append("")
lines.append(f"- Total cells: {len(rows)}")
lines.append(f"- Skipped (L == N, no second cycle): {n_skipped}")
lines.append(f"- Eligible cells (L > N): {n_eligible}")
lines.append(f"- Test 3a passed (peak_lag is multiple of N): {n_3a_pass} / {n_eligible}")
lines.append(f"- Test 3b passed (peak_val within {TOL} of predicted): {n_3b_pass} / {n_eligible}")
lines.append(f"- Both tests passed: {n_combined_pass} / {n_eligible}")
lines.append("")

if failures:
    lines.append("## Remaining criterion-3 failures")
    lines.append("")
    lines.append("| N | eps | L | peak_lag | peak_val | predicted | 3a | 3b |")
    lines.append("|---|-----|---|----------|----------|-----------|----|----|")
    for r in failures:
        lines.append(
            f"| {r['N']} | {r['epsilon']} | {r['L']} | {r['peak_lag']} | "
            f"{r['peak_val']:.5f} | {r['predicted_match_prob']:.5f} | "
            f"{'PASS' if r['test_3a'] else 'FAIL'} | "
            f"{'PASS' if r['test_3b'] else 'FAIL'} |"
        )
    lines.append("")
else:
    lines.append("## Remaining criterion-3 failures")
    lines.append("")
    lines.append("None. All eligible cells pass both 3a and 3b.")
    lines.append("")

# tabulate skipped cells for transparency
skipped_rows = [r for r in rows if r.get("status", "").startswith("L=")]
if skipped_rows:
    lines.append("## Skipped cells (L == N)")
    lines.append("")
    lines.append("These have only one polymer cycle, so lag-N autocorrelation has zero")
    lines.append("pairs to average. The autocorrelation engine returned a small-lag peak")
    lines.append("with low value (consistent with no detectable periodicity from a single cycle).")
    lines.append("")
    lines.append("| N | eps | L | peak_lag | peak_val |")
    lines.append("|---|-----|---|----------|----------|")
    for r in skipped_rows:
        lines.append(f"| {r['N']} | {r['epsilon']} | {r['L']} | {r['peak_lag']} | {r['peak_val']:.5f} |")
    lines.append("")

lines.append("## Final verdict")
lines.append("")
lines.append(f"# **Test B: {'OVERALL PASS' if overall else 'OVERALL FAIL'}**")
lines.append("")
if overall:
    lines.append(
        "All three criteria pass under the corrected criterion-3 logic. The "
        "framework's prediction that Y is periodic with period N is "
        "borne out: peaks land at integer multiples of N (3a), and the lag-N "
        "autocorrelation matches the noise-model prediction "
        "`(1-eps)^2 + eps^2/(N-1)` to within sampling tolerance (3b)."
    )

OUT_PATH.write_text("\n".join(lines) + "\n")
print("\n".join(lines))
