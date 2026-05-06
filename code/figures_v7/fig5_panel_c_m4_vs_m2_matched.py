"""fig 5 panel c: M4 vs M2 head-to-head at matched per-offspring drift.

source: results/test_h5_M2_vs_M4_matched_v1.csv
takes the n_m2_init=200, n_m4_init=200 row (the matched-population condition):
  - plateau heights: mean_final_fitness_m2 = 0.542, mean_final_fitness_m4 = 0.997
  - direct competition: p_m4_wins = 0.667, p_m2_wins = 0.333
output: paper/figures/v7/fig5/fig5_panel_c.{pdf,png}

annotation: "M2 drift rate is 100x M4's, but M4 still wins."
M2 at r*=0.10 introduces ~2.40 new positions/offspring; M4 at mu*=0.001
introduces ~0.02. So M2 has 100x more random search yet M4 wins.
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from templating_style import apply_paper_style, OKABE_ITO  # noqa: E402

apply_paper_style()

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
OUT_DIR = PROJECT_ROOT / "paper" / "figures" / "v7" / "fig5"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CSV_PATH = RESULTS_DIR / "test_h5_M2_vs_M4_matched_v1.csv"
PDF_PATH = OUT_DIR / "fig5_panel_c.pdf"
PNG_PATH = OUT_DIR / "fig5_panel_c.png"


def loadMatched():
    """find the 200/200 balanced row."""
    with open(CSV_PATH, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if int(row["n_m2_init"]) == 200 and int(row["n_m4_init"]) == 200:
                return {
                    "fit_m2": float(row["mean_final_fitness_m2"]),
                    "fit_m4": float(row["mean_final_fitness_m4"]),
                    "win_m2": float(row["p_m2_wins"]),
                    "win_m4": float(row["p_m4_wins"]),
                    "n": int(row["n_replicates"]),
                }
    raise RuntimeError("200/200 row not found in matched csv")


def plotPanelC(d):
    fig, ax = plt.subplots(figsize=(85 / 25.4, 60 / 25.4))

    # two metric groups: plateau height and win rate
    group_centers = np.array([0.0, 1.4])
    bar_w = 0.45
    # within each group, M4 then M2
    x_m4 = group_centers - bar_w / 2
    x_m2 = group_centers + bar_w / 2

    vals_m4 = np.array([d["fit_m4"], d["win_m4"]])
    vals_m2 = np.array([d["fit_m2"], d["win_m2"]])

    c_m4 = OKABE_ITO["blue"]
    c_m2 = OKABE_ITO["orange"]

    ax.bar(x_m4, vals_m4, width=bar_w, color=c_m4, edgecolor="black",
           linewidth=0.6, label=r"M4 ($\mu^*$=0.001)")
    ax.bar(x_m2, vals_m2, width=bar_w, color=c_m2, edgecolor="black",
           linewidth=0.6, label=r"M2 ($r^*$=0.10)")

    # value labels above bars
    for x, v in zip(x_m4, vals_m4):
        ax.text(x, v + 0.02, f"{v:.2f}", ha="center", va="bottom",
                fontsize=7, color=c_m4)
    for x, v in zip(x_m2, vals_m2):
        ax.text(x, v + 0.02, f"{v:.2f}", ha="center", va="bottom",
                fontsize=7, color=c_m2)

    ax.set_xticks(group_centers)
    ax.set_xticklabels(["plateau\nfitness", "direct-competition\nwin rate"])
    ax.set_ylabel("value")
    ax.set_ylim(0.0, 1.20)
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])

    # annotation: M2 drift 100x M4 yet M4 wins
    ax.text(
        0.5, 1.13,
        r"M2 drift rate is 100$\times$ M4's, but M4 still wins",
        transform=ax.transAxes, ha="center", va="center",
        fontsize=7, color="black",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor="grey", linewidth=0.6),
    )

    ax.legend(loc="upper right", frameon=False, handlelength=1.2,
              handletextpad=0.4, borderaxespad=0.3, fontsize=7)

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
        subprocess.run(["pdfinfo", str(pdf_path)], check=True,
                       capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        failures.append(f"pdfinfo failed: {e.stderr}")
    try:
        with Image.open(png_path) as im:
            im.verify()
    except Exception as e:
        failures.append(f"png verify: {e}")
    for ax in fig.axes:
        if not ax.get_xlabel().strip() and not ax.get_xticklabels():
            failures.append("empty xlabel and no xticklabels")
        if not ax.get_ylabel().strip():
            failures.append("empty ylabel")
    min_font = min(plt.rcParams["font.size"], plt.rcParams["axes.labelsize"],
                   plt.rcParams["xtick.labelsize"], plt.rcParams["ytick.labelsize"],
                   plt.rcParams["legend.fontsize"])
    if min_font < 7:
        failures.append(f"min font {min_font} < 7pt")
    for p in (pdf_path, png_path):
        size_mb = p.stat().st_size / 1e6
        if size_mb >= 5.0:
            failures.append(f"{p.name} {size_mb:.2f} MB >= 5")
    if failures:
        raise RuntimeError("validation failed: " + "; ".join(failures))


def main():
    d = loadMatched()
    print(f"  data: M4 plateau={d['fit_m4']:.4f}, M2 plateau={d['fit_m2']:.4f}, "
          f"M4 wins {d['win_m4']:.2%}, M2 wins {d['win_m2']:.2%} "
          f"(n_rep={d['n']})")

    fig = plotPanelC(d)
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=600)
    validate(PDF_PATH, PNG_PATH, fig)
    plt.close(fig)
    print(f"OK panel C:")
    print(f"  pdf: {PDF_PATH} ({PDF_PATH.stat().st_size / 1024:.1f} kB)")
    print(f"  png: {PNG_PATH} ({PNG_PATH.stat().st_size / 1024:.1f} kB)")


if __name__ == "__main__":
    main()
