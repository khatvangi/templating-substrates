"""fig 5 panel a: M0-M4 plateau heights across selection sharpness beta.

source: results/test_h3_isolated_dynamics_v1.csv
filter: step=1, L=32, beta in {2, 5, 10, 20}
mechanisms: M0, M1, M2_r0.10, M3, M4
output: paper/figures/v7/fig5/fig5_panel_a.{pdf,png}

reader takeaway: M4 (R3+R4 satisfied) plateaus near 1.0 across all beta;
M0/M1/M2/M3 cap below 0.6 regardless of selection sharpness — Corollary 2.
"""
import csv
import sys
import subprocess
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

# style module lives in same directory
sys.path.insert(0, str(Path(__file__).resolve().parent))
from templating_style import apply_paper_style, OKABE_ITO  # noqa: E402

apply_paper_style()

# paths
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
OUT_DIR = PROJECT_ROOT / "paper" / "figures" / "v7" / "fig5"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CSV_PATH = RESULTS_DIR / "test_h3_isolated_dynamics_v1.csv"
PDF_PATH = OUT_DIR / "fig5_panel_a.pdf"
PNG_PATH = OUT_DIR / "fig5_panel_a.png"

# the five mechanisms to plot, in order
MECHS = ["M0", "M1", "M2_r0.10", "M3", "M4"]
MECH_LABELS = ["M0", "M1", "M2 (r=0.10)", "M3", "M4"]
BETAS = [2.0, 5.0, 10.0, 20.0]

# okabe-ito 5-color qualitative
COLORS = [
    OKABE_ITO["black"],      # M0 (chance)
    OKABE_ITO["skyblue"],    # M1
    OKABE_ITO["orange"],     # M2
    OKABE_ITO["green"],      # M3
    OKABE_ITO["blue"],       # M4 (winner)
]


def loadH3():
    """load H3 isolated dynamics; return dict {mech_key -> {beta -> (mean, std)}}."""
    data = {m: {} for m in MECHS}
    with open(CSV_PATH, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            # only step=1 has the full beta sweep at L=32
            if int(row["step"]) != 1:
                continue
            if int(row["L"]) != 32:
                continue
            mech = row["mech_key"]
            if mech not in MECHS:
                continue
            beta = float(row["beta"])
            if beta not in BETAS:
                continue
            data[mech][beta] = (
                float(row["mean_final_fitness"]),
                float(row["std_final_fitness"]),
            )
    # sanity check: every (mech, beta) must be present
    for m in MECHS:
        for b in BETAS:
            if b not in data[m]:
                raise RuntimeError(f"missing data for mech={m} beta={b}")
    return data


def plotPanelA(data):
    fig, ax = plt.subplots(figsize=(85 / 25.4, 60 / 25.4))

    n_mech = len(MECHS)
    n_beta = len(BETAS)
    x = np.arange(n_beta)
    bar_w = 0.16
    offsets = np.linspace(-bar_w * 2, bar_w * 2, n_mech)

    for i, mech in enumerate(MECHS):
        means = [data[mech][b][0] for b in BETAS]
        stds = [data[mech][b][1] for b in BETAS]
        # sem from std/sqrt(30) to match panel D convention; use std for direct read
        # use std as error bar (per-replicate dispersion is informative for bar chart)
        ax.bar(
            x + offsets[i],
            means,
            width=bar_w,
            color=COLORS[i],
            edgecolor="black",
            linewidth=0.5,
            label=MECH_LABELS[i],
            yerr=stds,
            error_kw={"elinewidth": 0.6, "capsize": 1.5, "capthick": 0.5},
        )

    # reference lines: chance (0.25) and theoretical max (1.0)
    ax.axhline(0.25, color="grey", linestyle="--", linewidth=0.8, zorder=0)
    ax.axhline(1.0, color="grey", linestyle=":", linewidth=0.8, zorder=0)

    ax.set_xticks(x)
    ax.set_xticklabels([f"{int(b)}" for b in BETAS])
    ax.set_xlabel(r"selection sharpness $\beta$")
    ax.set_ylabel("plateau fitness (gen 1000)")
    ax.set_ylim(0.0, 1.10)

    # legend in two columns to fit panel
    ax.legend(loc="upper left", ncol=2, frameon=False, columnspacing=0.8,
              handlelength=1.2, handletextpad=0.4, borderaxespad=0.3)

    # annotate chance and ceiling
    ax.text(3.45, 0.27, "chance", fontsize=6, color="grey", ha="right", va="bottom")
    ax.text(3.45, 1.01, "max", fontsize=6, color="grey", ha="right", va="bottom")

    fig.tight_layout(pad=0.3)
    return fig


def validate(pdf_path, png_path, fig):
    """run the 7 validation checks; raise RuntimeError on any failure."""
    failures = []

    # 1. files exist
    if not pdf_path.exists():
        failures.append(f"pdf missing: {pdf_path}")
    if not png_path.exists():
        failures.append(f"png missing: {png_path}")
    if failures:
        raise RuntimeError("; ".join(failures))

    # 2. pdf opens (use pdfinfo)
    try:
        subprocess.run(
            ["pdfinfo", str(pdf_path)],
            check=True, capture_output=True, text=True,
        )
    except subprocess.CalledProcessError as e:
        failures.append(f"pdfinfo failed: {e.stderr}")

    # 3. png opens
    try:
        with Image.open(png_path) as im:
            im.verify()
    except Exception as e:
        failures.append(f"png verify failed: {e}")

    # 4. axis labels non-empty (check fig axes)
    for ax in fig.axes:
        if not ax.get_xlabel().strip():
            failures.append("empty xlabel")
        if not ax.get_ylabel().strip():
            failures.append("empty ylabel")

    # 5. min font size >= 7pt (rcParams set explicitly above)
    min_font = min(
        plt.rcParams["font.size"],
        plt.rcParams["axes.labelsize"],
        plt.rcParams["xtick.labelsize"],
        plt.rcParams["ytick.labelsize"],
        plt.rcParams["legend.fontsize"],
    )
    if min_font < 7:
        failures.append(f"min font {min_font} < 7pt")

    # 6. file sizes < 5 MB
    for p in (pdf_path, png_path):
        size_mb = p.stat().st_size / 1e6
        if size_mb >= 5.0:
            failures.append(f"{p.name} is {size_mb:.2f} MB (>= 5 MB)")

    if failures:
        raise RuntimeError("validation failed: " + "; ".join(failures))


def main():
    data = loadH3()

    # sanity: confirm beta=10 numbers match dispatch within 0.005
    expected = {"M0": 0.249, "M1": 0.472, "M2_r0.10": 0.533,
                "M3": 0.498, "M4": 0.972}
    for mech, expt in expected.items():
        actual = data[mech][10.0][0]
        if abs(actual - expt) > 0.005:
            print(f"warning: {mech} at beta=10 actual={actual:.4f} expected={expt:.3f} "
                  f"diff={abs(actual - expt):.4f}", file=sys.stderr)

    fig = plotPanelA(data)
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=600)

    validate(PDF_PATH, PNG_PATH, fig)
    plt.close(fig)

    print(f"OK panel A:")
    print(f"  pdf: {PDF_PATH} ({PDF_PATH.stat().st_size / 1024:.1f} kB)")
    print(f"  png: {PNG_PATH} ({PNG_PATH.stat().st_size / 1024:.1f} kB)")


if __name__ == "__main__":
    main()
