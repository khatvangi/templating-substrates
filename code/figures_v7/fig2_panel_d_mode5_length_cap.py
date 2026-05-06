"""fig 2 panel d: Mode 5 (modular conveyor) output length capped at module count N.

source: results/test_c_length_limit.csv

The CSV stores per-position information for two representative cases,
N=4 (A=4) and N=8 (A=20), each across 2N positions. Templated positions
(index < N) carry ~log_2(A) bits; positions index >= N carry ~0 bits because
no module exists past N.

We translate this into the dispatch's "requested L vs actual templated
length" view: for a Mode 5 system with N modules, requesting an output of
length L produces an actual templated length of min(L, N). The data confirms
the cap empirically -- positions past N drop to bias-floor information.

Note: the dispatch asks for N in {3, 5, 10}. The CSV only has N=4 and N=8
(both empirically validated). We use the available cases since the data is
load-bearing here -- per the dispatch's "trust the CSV, not the dispatch's
stated values" rule.

output: paper/figures/v7/fig2/fig2_panel_d.{pdf,png}

reader takeaway: Mode 5 (R2 fails -- no separable template) caps actual
output length at the module count N. The reference y = x line shows the
unbounded case (Mode 1); the staircase clamping at y = N is the Mode 5
length-limit signature.
"""
import csv
import sys
import subprocess
from pathlib import Path
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
from templating_style import apply_paper_style, OKABE_ITO  # noqa: E402

apply_paper_style()

# paths
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
OUT_DIR = PROJECT_ROOT / "paper" / "figures" / "v7" / "fig2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CSV_PATH = RESULTS_DIR / "test_c_length_limit.csv"
PDF_PATH = OUT_DIR / "fig2_panel_d.pdf"
PNG_PATH = OUT_DIR / "fig2_panel_d.png"

# bias-floor info threshold: a position with I_per_position above this is
# considered "templated" (carrying X-information). Templated values are
# ~2 bits (A=4) and ~4.2 bits (A=20); untemplated values are ~0.001-0.03.
TEMPLATE_THRESHOLD = 0.5


def loadLengthLimit():
    """parse test_c_length_limit.csv -> {case_label: dict(N, A, positions=[(i, I)])}"""
    cases = defaultdict(list)
    with open(CSV_PATH, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            cases[row["case_label"]].append((
                int(row["position_index"]),
                float(row["I_per_position_empirical"]),
                row["is_templated_region"] == "True",
            ))

    parsed = {}
    for label, rows in cases.items():
        # parse "N=<n>_A=<a>" label
        parts = label.split("_")
        N = int(parts[0].split("=")[1])
        A = int(parts[1].split("=")[1])
        rows.sort()
        parsed[label] = {
            "N": N,
            "A": A,
            "positions": rows,  # list of (i, I, is_templated_csv_label)
        }
    return parsed


def actualLengthCurve(case, L_max):
    """for a Mode 5 system with module count N, compute the empirical
    "actual templated length" given a requested output length L.

    The CSV gives per-position info up to 2N. The actual templated length at
    requested L is the count of positions in [0, L-1] whose info exceeds the
    bias floor. This naturally clamps at N for L >= N because positions past
    N have ~0 info."""
    N = case["N"]
    info_by_pos = {i: I for (i, I, _) in case["positions"]}

    Ls = np.arange(0, L_max + 1)
    actual = np.zeros_like(Ls, dtype=float)
    for j, L in enumerate(Ls):
        if L == 0:
            actual[j] = 0.0
            continue
        # for positions beyond what the CSV covers (i >= 2N), use the same
        # logic: untemplated tail. The empirical data confirms info ~ 0 for
        # i >= N regardless of how far we extend.
        count = 0
        for i in range(L):
            if i in info_by_pos:
                if info_by_pos[i] > TEMPLATE_THRESHOLD:
                    count += 1
            else:
                # past the CSV's range: by the panel D claim itself, info is
                # ~0 for i >= N; do not count
                pass
        actual[j] = float(count)
    return Ls, actual


def plotPanelD(parsed):
    fig, ax = plt.subplots(figsize=(85 / 25.4, 60 / 25.4))

    cividis = cm.get_cmap("cividis")
    sorted_cases = sorted(parsed.items(), key=lambda kv: kv[1]["N"])
    n_lines = len(sorted_cases)
    colors = [cividis(t) for t in np.linspace(0.15, 0.85, n_lines)]

    L_max_global = max(2 * c["N"] for _, c in sorted_cases)
    L_dense = np.linspace(0, L_max_global, 200)

    # reference y = x line for the unbounded case (Mode 1 / Mode 5 if no cap)
    ax.plot(L_dense, L_dense, color="black", linestyle="--", linewidth=1.0,
            label="$y=L$ (unbounded)", zorder=2)

    for i, (label, case) in enumerate(sorted_cases):
        N = case["N"]
        A = case["A"]
        Ls, actual = actualLengthCurve(case, 2 * N)
        legend_label = f"$N={N}$ ($A={A}$)"
        ax.plot(Ls, actual, color=colors[i], marker="o", markersize=3.0,
                linewidth=1.4, label=legend_label, zorder=3 + i)
        # horizontal cap line at y = N
        ax.axhline(N, color=colors[i], linestyle=":", linewidth=0.7,
                   alpha=0.7, zorder=1)

    ax.set_xlabel(r"requested output length $L$")
    ax.set_ylabel("actual templated length")
    ax.set_xlim(0, L_max_global)
    ax.set_ylim(0, L_max_global * 1.02)

    ax.legend(loc="lower right", frameon=False,
              handlelength=1.5, handletextpad=0.5, borderaxespad=0.3,
              labelspacing=0.3)

    fig.tight_layout(pad=0.3)
    return fig


def validate(pdf_path, png_path, fig):
    failures = []
    if not pdf_path.exists():
        failures.append(f"pdf missing: {pdf_path}")
    if not png_path.exists():
        failures.append(f"png missing: {png_path}")
    if failures:
        raise RuntimeError("; ".join(failures))

    try:
        subprocess.run(["pdfinfo", str(pdf_path)],
                       check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        failures.append(f"pdfinfo failed: {e.stderr}")

    try:
        with Image.open(png_path) as im:
            im.verify()
    except Exception as e:
        failures.append(f"png verify failed: {e}")

    for ax in fig.axes:
        if not ax.get_xlabel().strip():
            failures.append("empty xlabel")
        if not ax.get_ylabel().strip():
            failures.append("empty ylabel")

    min_font = min(
        plt.rcParams["font.size"],
        plt.rcParams["axes.labelsize"],
        plt.rcParams["xtick.labelsize"],
        plt.rcParams["ytick.labelsize"],
        plt.rcParams["legend.fontsize"],
    )
    if min_font < 7:
        failures.append(f"min font {min_font} < 7pt")

    for p in (pdf_path, png_path):
        size_mb = p.stat().st_size / 1e6
        if size_mb >= 5.0:
            failures.append(f"{p.name} is {size_mb:.2f} MB (>= 5 MB)")

    if failures:
        raise RuntimeError("validation failed: " + "; ".join(failures))


def main():
    parsed = loadLengthLimit()
    if not parsed:
        raise RuntimeError("no cases parsed from length-limit CSV")
    print(f"  loaded {len(parsed)} length-limit cases:")
    for label, case in parsed.items():
        print(f"    {label}: N={case['N']}, A={case['A']}, "
              f"{len(case['positions'])} positions")

    fig = plotPanelD(parsed)
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=600)
    validate(PDF_PATH, PNG_PATH, fig)
    plt.close(fig)
    print(f"OK panel D:")
    print(f"  pdf: {PDF_PATH} ({PDF_PATH.stat().st_size / 1024:.1f} kB)")
    print(f"  png: {PNG_PATH} ({PNG_PATH.stat().st_size / 1024:.1f} kB)")


if __name__ == "__main__":
    main()
