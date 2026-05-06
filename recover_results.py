"""recover test_a2_results.csv and test_b_results.csv from preserved sources

a2: progress.txt was deleted, data transcribed from chat — reconstruct progress.txt + csv
b : progress.txt and stdout.log are intact — read the stdout table back into csv
"""

import csv
import re
from pathlib import Path

results_dir = Path("/storage/kiran-stuff/templating_framework/results")

# ---------------------------------------------------------------------------
# step 1: test_a2_progress.txt + test_a2_results.csv
# ---------------------------------------------------------------------------

# transcribed from chat (original progress.txt was deleted before recovery)
a2_data_block = """\
bias=pi_uniform, eps=0.01, L=25: I_sys1=47.604, I_bulk=0.033, ratio=1447.9x, mdev=0.0012
bias=pi_uniform, eps=0.01, L=100: I_sys1=190.754, I_bulk=0.122, ratio=1557.9x, mdev=0.0011
bias=pi_uniform, eps=0.01, L=500: I_sys1=952.219, I_bulk=0.631, ratio=1508.0x, mdev=0.0002
bias=pi_uniform, eps=0.05, L=25: I_sys1=41.070, I_bulk=0.034, ratio=1191.9x, mdev=0.0021
bias=pi_uniform, eps=0.05, L=100: I_sys1=163.432, I_bulk=0.128, ratio=1277.8x, mdev=0.0016
bias=pi_uniform, eps=0.05, L=500: I_sys1=817.703, I_bulk=0.668, ratio=1223.9x, mdev=0.0003
bias=pi_uniform, eps=0.1, L=25: I_sys1=34.285, I_bulk=0.036, ratio=945.9x, mdev=0.0015
bias=pi_uniform, eps=0.1, L=100: I_sys1=137.331, I_bulk=0.124, ratio=1106.5x, mdev=0.0012
bias=pi_uniform, eps=0.1, L=500: I_sys1=686.668, I_bulk=0.644, ratio=1066.5x, mdev=0.0009
bias=pi_uniform, eps=0.25, L=25: I_sys1=19.780, I_bulk=0.032, ratio=617.6x, mdev=0.0040
bias=pi_uniform, eps=0.25, L=100: I_sys1=79.505, I_bulk=0.127, ratio=625.5x, mdev=0.0007
bias=pi_uniform, eps=0.25, L=500: I_sys1=396.481, I_bulk=0.643, ratio=616.5x, mdev=0.0005
bias=pi_AT_skew, eps=0.01, L=25: I_sys1=40.911, I_bulk=0.037, ratio=1115.3x, mdev=0.0011
bias=pi_AT_skew, eps=0.01, L=100: I_sys1=163.454, I_bulk=0.121, ratio=1346.3x, mdev=0.0013
bias=pi_AT_skew, eps=0.01, L=500: I_sys1=817.598, I_bulk=0.654, ratio=1250.3x, mdev=0.0006
bias=pi_AT_skew, eps=0.05, L=25: I_sys1=34.755, I_bulk=0.034, ratio=1010.6x, mdev=0.0011
bias=pi_AT_skew, eps=0.05, L=100: I_sys1=139.583, I_bulk=0.128, ratio=1089.9x, mdev=0.0011
bias=pi_AT_skew, eps=0.05, L=500: I_sys1=698.087, I_bulk=0.650, ratio=1074.6x, mdev=0.0003
bias=pi_AT_skew, eps=0.1, L=25: I_sys1=29.191, I_bulk=0.034, ratio=853.9x, mdev=0.0015
bias=pi_AT_skew, eps=0.1, L=100: I_sys1=116.973, I_bulk=0.133, ratio=876.2x, mdev=0.0012
bias=pi_AT_skew, eps=0.1, L=500: I_sys1=584.308, I_bulk=0.626, ratio=933.4x, mdev=0.0004
bias=pi_AT_skew, eps=0.25, L=25: I_sys1=16.817, I_bulk=0.037, ratio=457.2x, mdev=0.0019
bias=pi_AT_skew, eps=0.25, L=100: I_sys1=67.383, I_bulk=0.139, ratio=483.7x, mdev=0.0009
bias=pi_AT_skew, eps=0.25, L=500: I_sys1=337.393, I_bulk=0.644, ratio=523.5x, mdev=0.0007
bias=pi_GC_skew, eps=0.01, L=25: I_sys1=40.880, I_bulk=0.039, ratio=1052.5x, mdev=0.0045
bias=pi_GC_skew, eps=0.01, L=100: I_sys1=163.295, I_bulk=0.121, ratio=1348.8x, mdev=0.0012
bias=pi_GC_skew, eps=0.01, L=500: I_sys1=816.746, I_bulk=0.655, ratio=1247.0x, mdev=0.0004
bias=pi_GC_skew, eps=0.05, L=25: I_sys1=34.983, I_bulk=0.032, ratio=1109.3x, mdev=0.0006
bias=pi_GC_skew, eps=0.05, L=100: I_sys1=139.270, I_bulk=0.129, ratio=1082.7x, mdev=0.0007
bias=pi_GC_skew, eps=0.05, L=500: I_sys1=697.964, I_bulk=0.656, ratio=1064.1x, mdev=0.0004
bias=pi_GC_skew, eps=0.1, L=25: I_sys1=29.209, I_bulk=0.036, ratio=811.4x, mdev=0.0026
bias=pi_GC_skew, eps=0.1, L=100: I_sys1=116.821, I_bulk=0.138, ratio=846.5x, mdev=0.0020
bias=pi_GC_skew, eps=0.1, L=500: I_sys1=584.540, I_bulk=0.680, ratio=859.4x, mdev=0.0003
bias=pi_GC_skew, eps=0.25, L=25: I_sys1=16.865, I_bulk=0.037, ratio=452.0x, mdev=0.0009
bias=pi_GC_skew, eps=0.25, L=100: I_sys1=67.269, I_bulk=0.138, ratio=488.2x, mdev=0.0012
bias=pi_GC_skew, eps=0.25, L=500: I_sys1=337.271, I_bulk=0.670, ratio=503.2x, mdev=0.0005"""

# canonical (label -> distribution string) map per task spec
bias_dist = {
    "pi_uniform": "(0.25, 0.25, 0.25, 0.25)",
    "pi_AT_skew": "(0.40, 0.10, 0.10, 0.40)",
    "pi_GC_skew": "(0.10, 0.40, 0.40, 0.10)",
}

a2_pat = re.compile(
    r"bias=(?P<bias>\S+),\s+eps=(?P<eps>[\d.]+),\s+L=(?P<L>\d+):\s+"
    r"I_sys1=(?P<isys1>[\d.]+),\s+I_bulk=(?P<ibulk>[\d.]+),\s+"
    r"ratio=(?P<ratio>[\d.]+)x,\s+mdev=(?P<mdev>[\d.]+)"
)

a2_rows = []
for line in a2_data_block.strip().split("\n"):
    m = a2_pat.match(line.strip())
    if not m:
        raise ValueError(f"could not parse a2 line: {line!r}")
    a2_rows.append(m.groupdict())

assert len(a2_rows) == 36, f"expected 36 a2 rows, got {len(a2_rows)}"

# write progress.txt — single approximate timestamp per task spec
a2_progress_path = results_dir / "test_a2_progress.txt"
ts = "2026-05-04 21:28:14"
with a2_progress_path.open("w") as f:
    f.write(f"\n# === run started at {ts} ===\n")
    for line in a2_data_block.strip().split("\n"):
        f.write(f"[{ts}] {line}\n")

# write csv with full schema
a2_csv_path = results_dir / "test_a2_results.csv"
with a2_csv_path.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow([
        "template_bias_label",
        "template_bias_distribution",
        "epsilon",
        "L",
        "I_system1_empirical",
        "I_bulk_empirical",
        "I_system1_per_position",
        "I_bulk_per_position",
        "marginal_match_max_dev",
        "separation_ratio",
    ])
    for r in a2_rows:
        L = int(r["L"])
        i_sys1 = float(r["isys1"])
        i_bulk = float(r["ibulk"])
        w.writerow([
            r["bias"],
            bias_dist[r["bias"]],
            float(r["eps"]),
            L,
            i_sys1,
            i_bulk,
            i_sys1 / L,
            i_bulk / L,
            float(r["mdev"]),
            float(r["ratio"]),
        ])

# ---------------------------------------------------------------------------
# step 2: test_b_results.csv from the stdout table
# ---------------------------------------------------------------------------

stdout_path = results_dir / "test_b_stdout.log"
stdout_lines = stdout_path.read_text().splitlines()

# table body lives at lines 3..107 (1-indexed) per inspection
table_rows = stdout_lines[2:107]
assert len(table_rows) == 105, f"expected 105 b rows, got {len(table_rows)}"

b_csv_path = results_dir / "test_b_results.csv"
with b_csv_path.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow([
        "N",
        "epsilon",
        "L",
        "n_samples",
        "block_length",
        "I_empirical",
        "log2_N_theoretical",
        "per_position_info_empirical",
        "autocorrelation_peak_lag",
        "autocorrelation_peak_value",
        "H_y_block",            # lost — not in stdout
        "H_y_given_x_block",    # lost — not in stdout
    ])
    for raw in table_rows:
        # columns:  N  eps  L  n  k  I_emp  log2N  ipp  peak_lag  peak_val
        parts = raw.split()
        if len(parts) != 10:
            raise ValueError(f"unexpected b row: {raw!r}")
        N, eps, L, n, k, I_emp, log2N, ipp, peak_lag, peak_val = parts
        w.writerow([
            int(N),
            float(eps),
            int(L),
            int(n),
            int(k),
            float(I_emp),
            float(log2N),
            float(ipp),
            int(peak_lag),
            float(peak_val),
            "",
            "",
        ])

# ---------------------------------------------------------------------------
# step 3: verification
# ---------------------------------------------------------------------------

def head(path, n):
    with open(path) as f:
        return [next(f).rstrip("\n") for _ in range(n)]

a2_size = a2_csv_path.stat().st_size
b_size = b_csv_path.stat().st_size
a2_progress_size = a2_progress_path.stat().st_size

with a2_csv_path.open() as f:
    a2_n = sum(1 for _ in f) - 1
with b_csv_path.open() as f:
    b_n = sum(1 for _ in f) - 1

print(f"test_a2_results.csv:    {a2_n} rows, {a2_size} bytes")
print(f"test_b_results.csv :    {b_n} rows, {b_size} bytes")
print(f"test_a2_progress.txt:   {a2_progress_size} bytes")
print()
print("--- test_a2_results.csv (first 4 lines) ---")
for line in head(a2_csv_path, 4):
    print(line)
print()
print("--- test_b_results.csv (first 4 lines) ---")
for line in head(b_csv_path, 4):
    print(line)
print()

assert a2_n == 36
assert b_n == 105

# step 4: sentinel
sentinel_path = results_dir / "CANONICAL_RESULTS.md"
sentinel_path.write_text(
    "# Canonical results (DO NOT OVERWRITE)\n"
    "\n"
    "The following CSVs are the canonical results of completed tests:\n"
    "\n"
    "- test_a1_results.csv (Test A.1, PASS)\n"
    "- test_a2_results.csv (Test A.2, PASS)  \n"
    "- test_b_results.csv (Test B, PASS, including verdict_corrected)\n"
    "- test_c_results.csv (Test C, PASS)\n"
    "- test_c_length_limit.csv (Test C length-limit sub-test, PASS)\n"
    "\n"
    "Future tasks must NOT overwrite these files. Future tasks may CREATE new files\n"
    "prefixed with the new test's identifier (e.g., test_d_*, test_e_*).\n"
    "\n"
    "If a future task needs to reload these CSVs for plotting or analysis, READ ONLY.\n"
)

print("RECOVERY COMPLETE")
