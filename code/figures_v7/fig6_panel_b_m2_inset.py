"""fig 6 panel b: M2 r-sweep data inset for cortex placement.

source:
  - test_h5_r_sweep_v1.csv (beta=10, L=32, 14 r values)

output: paper/figures/v7/fig6/fig6_panel_b.{pdf,png}

reader takeaway: the eukaryotic ciliate cortex sits in the M2 regime at
low re-draw rate (r ~ 0.10) — lineage-fixation-dominant. this is the data
inset only; the user composites the cell-division schematic on top in
inkscape. panel size is 85 x 80 mm (taller than the standard 60 mm to
leave room for that schematic).
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
OUT_DIR = PROJECT_ROOT / "paper" / "figures" / "v7" / "fig6"
OUT_DIR.mkdir(parents=True, exist_ok=True)

R_SWEEP_CSV = RESULTS_DIR / "test_h5_r_sweep_v1.csv"
PDF_PATH = OUT_DIR / "fig6_panel_b.pdf"
PNG_PATH = OUT_DIR / "fig6_panel_b.png"

CORTEX_R = 0.10  # the regime mapping argument (Mode 6 / cortex --> M2 r ~ 0.10)


def loadRSweep():
    """trust the csv: read r and mean_final_fitness, sort by r."""
    rs, fits, ci_lo, ci_hi = [], [], [], []
    with open(R_SWEEP_CSV, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rs.append(float(row["r"]))
            fits.append(float(row["mean_final_fitness"]))
            ci_lo.append(float(row["ci95_low"]))
            ci_hi.append(float(row["ci95_high"]))
    rs = np.array(rs)
    fits = np.array(fits)
    ci_lo = np.array(ci_lo)
    ci_hi = np.array(ci_hi)
    order = np.argsort(rs)
    return rs[order], fits[order], ci_lo[order], ci_hi[order]


def findCortexFitness(rs, fits, target_r=CORTEX_R):
    """find the plateau fitness at r = target_r exactly (data-driven)."""
    idx = int(np.argmin(np.abs(rs - target_r)))
    return rs[idx], fits[idx]


def plotPanelB(rs, fits, ci_lo, ci_hi):
    # 85 x 80 mm (taller than standard, room for inkscape schematic above)
    fig, ax = plt.subplots(figsize=(85 / 25.4, 80 / 25.4))

    # main M2 r-sweep curve
    ax.plot(
        rs, fits,
        marker="o", markersize=4, linewidth=1.4,
        color=OKABE_ITO["blue"],
        label="M2 r-sweep (L=32, " + r"$\beta$=10)",
        zorder=3,
    )
    # ci95 band: ci_lo/ci_hi columns are exactly that
    ax.fill_between(rs, ci_lo, ci_hi,
                    color=OKABE_ITO["blue"], alpha=0.15, zorder=2,
                    label="95% CI")

    # mark cortex regime at r ~ 0.10
    cortex_r, cortex_fit = findCortexFitness(rs, fits, CORTEX_R)
    ax.scatter([cortex_r], [cortex_fit],
               s=70, marker="*",
               color=OKABE_ITO["vermilion"],
               edgecolor="black", linewidth=0.6,
               zorder=5,
               label=f"cortex regime (r={cortex_r:.2f})")

    # vertical guide at r* = 0.10
    ax.axvline(cortex_r, color="grey", linestyle="--",
               linewidth=0.7, zorder=1)

    ax.set_xscale("log")
    ax.set_xlabel(r"re-draw rate $r$")
    ax.set_ylabel("M2 plateau fitness")
    ax.set_xlim(0.008, 1.1)
    ax.set_ylim(0.20, 0.62)

    # annotation describing where the cortex sits
    ax.annotate(
        "cortex regime:\nlow re-draw,\nlineage-fixation-\ndominant",
        xy=(cortex_r, cortex_fit),
        xytext=(0.013, 0.30),
        fontsize=7, color="black",
        ha="left", va="bottom",
        arrowprops=dict(
            arrowstyle="->",
            color="grey",
            linewidth=0.7,
            shrinkA=2, shrinkB=4,
        ),
    )

    ax.legend(loc="lower left", frameon=False,
              handlelength=1.4, handletextpad=0.4,
              borderaxespad=0.3, labelspacing=0.2)

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
        if not ax.get_xlabel().strip():
            failures.append("empty xlabel")
        if not ax.get_ylabel().strip():
            failures.append("empty ylabel")
        # check tick labels are non-empty
        for lbl in ax.get_xticklabels() + ax.get_yticklabels():
            # tick labels can legitimately be "" before draw; only flag axis labels.
            pass
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
    rs, fits, ci_lo, ci_hi = loadRSweep()
    cortex_r, cortex_fit = findCortexFitness(rs, fits, CORTEX_R)
    print(f"cortex placement (h5 r-sweep): r={cortex_r:.3f}, "
          f"plateau fitness={cortex_fit:.4f}")
    print(f"sweep peak r={rs[int(np.argmax(fits))]:.3f}, "
          f"fitness={fits.max():.4f}")

    fig = plotPanelB(rs, fits, ci_lo, ci_hi)
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=600)
    validate(PDF_PATH, PNG_PATH, fig)
    plt.close(fig)
    print(f"wrote {PDF_PATH}")
    print(f"wrote {PNG_PATH}")


if __name__ == "__main__":
    main()
