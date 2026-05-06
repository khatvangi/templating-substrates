"""fig 1 panel A: empirical vs analytical I_struct^pop calibration (test A1).

source: results/test_a1_results.csv
columns: epsilon, L, I_empirical, I_theoretical, I_per_position_*

scatter of empirical I_struct^pop vs analytical prediction; one color per
epsilon (sequential viridis) and one marker per L. y=x reference line.
annotate with Pearson r and OLS slope from a fit on all (analytical, empirical)
points.

NOTE on dispatch deviation: dispatch states epsilon in {0.001, 0.01, 0.1} and
L in {16, 32, 64, 128}; the actual CSV has epsilon in {0.001, 0.01, 0.05, 0.1,
0.25} and L in {10, 25, 50, 100, 200, 500}. per the pilot lesson, we trust the
CSV and plot all available (epsilon, L) combinations.

output: paper/figures/v7/fig1/fig1_panel_a.{pdf,png}

reader takeaway: empirical estimator tracks analytical prediction over four
decades of mutual information (~0.8 to ~1000 bits) with r ~= 1 and slope ~= 1.
"""
import csv
import sys
import subprocess
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
from PIL import Image

# style module lives in same directory
sys.path.insert(0, str(Path(__file__).resolve().parent))
from templating_style import apply_paper_style, OKABE_ITO  # noqa: E402, F401

apply_paper_style()

# paths
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
OUT_DIR = PROJECT_ROOT / "paper" / "figures" / "v7" / "fig1"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CSV_PATH = RESULTS_DIR / "test_a1_results.csv"
PDF_PATH = OUT_DIR / "fig1_panel_a.pdf"
PNG_PATH = OUT_DIR / "fig1_panel_a.png"

# marker per L (6 distinct values in the actual CSV)
L_MARKERS = {10: "o", 25: "s", 50: "^", 100: "D", 200: "v", 500: "P"}


def loadA1():
    """load A1 calibration table; return list of dicts."""
    rows = []
    with open(CSV_PATH, newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            rows.append({
                "epsilon": float(r["epsilon"]),
                "L": int(r["L"]),
                "I_emp": float(r["I_empirical"]),
                "I_thy": float(r["I_theoretical"]),
            })
    return rows


def fitStats(rows):
    """compute pearson r and OLS slope on the full (theoretical, empirical) set."""
    x = np.array([r["I_thy"] for r in rows])
    y = np.array([r["I_emp"] for r in rows])
    r = float(np.corrcoef(x, y)[0, 1])
    # OLS slope through origin would also work, but the dispatch says
    # "slope from a linear fit" — use ordinary least squares with intercept
    slope, intercept = np.polyfit(x, y, 1)
    return r, float(slope), float(intercept)


def plotPanelA(rows, r, slope, intercept):
    fig, ax = plt.subplots(figsize=(60 / 25.4, 60 / 25.4))

    # sequential viridis: one shade per epsilon
    epsilons = sorted({rr["epsilon"] for rr in rows})
    cmap = mpl.colormaps["viridis"]
    # spread shades across [0.15, 0.85] so we keep some contrast at both ends
    eps_color = {
        eps: cmap(0.15 + 0.7 * i / max(1, len(epsilons) - 1))
        for i, eps in enumerate(epsilons)
    }

    # plot in log-log so the four decades are visible (0.7 to ~1000 bits)
    for eps in epsilons:
        for L_val, marker in L_MARKERS.items():
            sel = [rr for rr in rows if rr["epsilon"] == eps and rr["L"] == L_val]
            if not sel:
                continue
            xs = [s["I_thy"] for s in sel]
            ys = [s["I_emp"] for s in sel]
            ax.scatter(
                xs, ys,
                color=eps_color[eps],
                marker=marker,
                s=18,
                edgecolor="black",
                linewidth=0.4,
                zorder=3,
            )

    # y = x reference line
    all_thy = [rr["I_thy"] for rr in rows]
    all_emp = [rr["I_emp"] for rr in rows]
    lo = min(min(all_thy), min(all_emp)) * 0.8
    hi = max(max(all_thy), max(all_emp)) * 1.2
    ax.plot([lo, hi], [lo, hi], color="grey", linestyle="--",
            linewidth=0.8, zorder=1, label="y = x")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel(r"analytical $I_{\mathrm{struct}}^{\mathrm{pop}}$ (bits)")
    ax.set_ylabel(r"empirical $I_{\mathrm{struct}}^{\mathrm{pop}}$ (bits)")

    # annotate r and slope in upper-left
    ax.text(
        0.04, 0.96,
        f"r = {r:.4f}\nslope = {slope:.3f}",
        transform=ax.transAxes,
        fontsize=7,
        ha="left", va="top",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                  edgecolor="black", linewidth=0.5),
    )

    # epsilon legend (color) — compact, lower-right
    eps_handles = [
        plt.Line2D([], [], marker="o", linestyle="",
                   markerfacecolor=eps_color[eps],
                   markeredgecolor="black",
                   markeredgewidth=0.4,
                   markersize=5,
                   label=rf"$\varepsilon$={eps:g}")
        for eps in epsilons
    ]
    ax.legend(
        handles=eps_handles,
        loc="lower right",
        frameon=False,
        handletextpad=0.3,
        labelspacing=0.2,
        borderaxespad=0.3,
        fontsize=6,
    )

    fig.tight_layout(pad=0.3)
    return fig


def validate(pdf_path, png_path, fig):
    """run validation checks; raise RuntimeError on any failure."""
    failures = []

    if not pdf_path.exists():
        failures.append(f"pdf missing: {pdf_path}")
    if not png_path.exists():
        failures.append(f"png missing: {png_path}")
    if failures:
        raise RuntimeError("; ".join(failures))

    try:
        subprocess.run(
            ["pdfinfo", str(pdf_path)],
            check=True, capture_output=True, text=True,
        )
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
    rows = loadA1()
    if not rows:
        print(f"ERROR: no rows loaded from {CSV_PATH}", file=sys.stderr)
        sys.exit(1)

    r, slope, intercept = fitStats(rows)

    fig = plotPanelA(rows, r, slope, intercept)
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=600)

    validate(PDF_PATH, PNG_PATH, fig)
    plt.close(fig)

    print("OK panel A:")
    print(f"  pdf: {PDF_PATH} ({PDF_PATH.stat().st_size / 1024:.1f} kB)")
    print(f"  png: {PNG_PATH} ({PNG_PATH.stat().st_size / 1024:.1f} kB)")
    print(f"  pearson r = {r:.6f}")
    print(f"  OLS slope = {slope:.6f}, intercept = {intercept:.6f}")
    print(f"  n_points  = {len(rows)}")


if __name__ == "__main__":
    main()
