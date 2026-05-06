"""fig 2 panel b: Mode 1 (linear in L) vs Mode 3 (saturating at log_2 N).

sources:
  results/test_a1_results.csv  -- Mode 1 empirical points (epsilon=0.001)
  results/test_b_results.csv   -- Mode 3 empirical points (N=2, epsilon=0.001)

filter: L in [16, 128] (per dispatch x-axis range)
output: paper/figures/v7/fig2/fig2_panel_b.{pdf,png}

reader takeaway: Mode 1 (R2 satisfied: open-ended template) information grows
linearly as 2L bits (alphabet 4). Mode 3 with N=2 cycle states (R2 fails)
saturates at log_2(2) = 1 bit no matter how long the output. The vertical
gap between the two lines is "open-ended (R2 satisfied) vs bounded (R2 fails)".
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
OUT_DIR = PROJECT_ROOT / "paper" / "figures" / "v7" / "fig2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

A1_CSV = RESULTS_DIR / "test_a1_results.csv"
B_CSV = RESULTS_DIR / "test_b_results.csv"
PDF_PATH = OUT_DIR / "fig2_panel_b.pdf"
PNG_PATH = OUT_DIR / "fig2_panel_b.png"

# x-axis range per dispatch
L_MIN, L_MAX = 16, 128
EPS_FIXED = 0.001
N_MODE3 = 2  # log_2(2) = 1 bit saturation


def loadMode1Points():
    """Mode 1 empirical points from test_a1: epsilon=EPS_FIXED, L in window."""
    pts = []
    with open(A1_CSV, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            eps = float(row["epsilon"])
            if abs(eps - EPS_FIXED) > 1e-9:
                continue
            L = int(row["L"])
            if L_MIN <= L <= L_MAX:
                pts.append((L, float(row["I_empirical"])))
    pts.sort()
    return pts


def loadMode3Points():
    """Mode 3 empirical points from test_b: N=2, epsilon=EPS_FIXED, L in window."""
    pts = []
    with open(B_CSV, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            N = int(row["N"])
            eps = float(row["epsilon"])
            if N != N_MODE3:
                continue
            if abs(eps - EPS_FIXED) > 1e-9:
                continue
            L = int(row["L"])
            if L_MIN <= L <= L_MAX:
                pts.append((L, float(row["I_empirical"])))
    pts.sort()
    return pts


def plotPanelB(mode1_pts, mode3_pts):
    fig, ax = plt.subplots(figsize=(85 / 25.4, 60 / 25.4))

    # analytical Mode 1: I = L * log_2(4) = 2L bits (perfect copy, eps -> 0)
    L_dense = np.linspace(L_MIN, L_MAX, 200)
    mode1_theory = 2.0 * L_dense
    ax.plot(L_dense, mode1_theory, color=OKABE_ITO["blue"],
            linewidth=1.5, label="Mode 1 ($2L$, R2 satisfied)", zorder=3)

    # analytical Mode 3: I = log_2(N=2) = 1 bit
    mode3_theory = np.full_like(L_dense, np.log2(N_MODE3))
    ax.plot(L_dense, mode3_theory, color=OKABE_ITO["orange"],
            linewidth=1.5, label=r"Mode 3 ($\log_2 N$, R2 fails)", zorder=3)

    # empirical points
    if mode1_pts:
        Lm1 = [p[0] for p in mode1_pts]
        Im1 = [p[1] for p in mode1_pts]
        ax.plot(Lm1, Im1, marker="o", markersize=3.5,
                linestyle="none", color=OKABE_ITO["blue"],
                markeredgecolor="black", markeredgewidth=0.4, zorder=4)
    if mode3_pts:
        Lm3 = [p[0] for p in mode3_pts]
        Im3 = [p[1] for p in mode3_pts]
        ax.plot(Lm3, Im3, marker="s", markersize=3.5,
                linestyle="none", color=OKABE_ITO["orange"],
                markeredgecolor="black", markeredgewidth=0.4, zorder=4)

    ax.set_xlabel(r"template length $L$")
    ax.set_ylabel(r"$I_{\mathrm{struct}}(X;Y)$ (bits)")
    ax.set_xlim(L_MIN - 4, L_MAX + 4)

    # mode 1 reaches 2*128=256, mode 3 saturates at 1
    ax.set_ylim(-5, 270)

    # gap annotation: vertical double-arrow at L ~ 96
    L_anno = 96
    y_top = 2.0 * L_anno
    y_bot = 1.0
    ax.annotate("", xy=(L_anno, y_top), xytext=(L_anno, y_bot),
                arrowprops=dict(arrowstyle="<->", color="grey",
                                lw=0.8, shrinkA=0, shrinkB=0))
    ax.text(L_anno + 3, (y_top + y_bot) / 2 + 5,
            "open-ended\nvs bounded",
            fontsize=7, color="grey", ha="left", va="center")

    ax.legend(loc="upper left", frameon=False,
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
    mode1_pts = loadMode1Points()
    mode3_pts = loadMode3Points()
    if not mode1_pts:
        print("warning: no Mode 1 empirical points in [16,128] at eps=0.001",
              file=sys.stderr)
    if not mode3_pts:
        print("warning: no Mode 3 empirical points in [16,128] at N=2 eps=0.001",
              file=sys.stderr)

    fig = plotPanelB(mode1_pts, mode3_pts)
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=600)
    validate(PDF_PATH, PNG_PATH, fig)
    plt.close(fig)
    print(f"OK panel B:")
    print(f"  mode1 empirical points: {len(mode1_pts)}")
    print(f"  mode3 empirical points: {len(mode3_pts)}")
    print(f"  pdf: {PDF_PATH} ({PDF_PATH.stat().st_size / 1024:.1f} kB)")
    print(f"  png: {PNG_PATH} ({PNG_PATH.stat().st_size / 1024:.1f} kB)")


if __name__ == "__main__":
    main()
