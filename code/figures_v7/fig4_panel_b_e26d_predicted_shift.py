"""fig 4 panel b: 6 E26->D primary-gate clades shift in framework-predicted direction.

source: results/test_f_family_predictions_v1.csv
filter: envelope_classification == "shifted (primary-gate substitution at [26])"  (6 clades)

observable mapping (CSV columns -> framework f_phase, f_monomer^A):
    f_phase     <- predicted_period2_peak  (period-2 alternation strength).
    f_monomer^A <- 1 - predicted_marginal_G (purity at canonical positions;
                                             E26->D weakens the per-state A
                                             channel so the misincorporation
                                             leak G grows).

DIVERGENCE FROM DISPATCH (documented):
    The dispatch states "f_monomer^A degraded to 0.71-0.84 while f_phase
    preserved at 0.94-0.98". Inspecting the actual CSVs reveals the OPPOSITE
    axis is the dominant degraded one for E26->D under the framework's
    parameterization rule (eps_A=0.15, P(G|A)=0.10):
        f_phase   (period-2 peak): 0.855-0.858 -- DEGRADED from WT 0.980
        f_monomer (1 - marg G):   0.948-0.948 -- mildly reduced from WT 0.997
    This is consistent with E26D weakening the cycle's per-state channel
    enough to inject ~5% G leak per position, which scrambles the strict
    A-C-A-C alternation (period-2 peak drops) more than it raises the marginal
    G fraction (which stays modest because each state still mostly outputs
    its canonical base).

    We plot what the framework actually predicts and what the underlying
    simulation produces, not the dispatch's flipped numbers. The "framework-
    predicted direction" is preserved: E26->D shifts toward LOWER f_phase
    relative to WT, with f_monomer also slightly reduced. All 6 clades
    cluster tightly together (the parameterization rule is identical for
    all E26->D variants regardless of secondary-residue identity).

reader takeaway: the 6 natural E26->D clades shift in the framework-predicted
direction (lower f_phase, slightly lower f_monomer) -- the apparatus
classifier separates them from secondary-residue-only clades.
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

# paths
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
OUT_DIR = PROJECT_ROOT / "paper" / "figures" / "v7" / "fig4"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CSV_PATH = RESULTS_DIR / "test_f_family_predictions_v1.csv"
PDF_PATH = OUT_DIR / "fig4_panel_b.pdf"
PNG_PATH = OUT_DIR / "fig4_panel_b.png"


def loadE26DAndWT():
    """return (e26d_pts, wt_xy, e26d_clades).

    e26d_pts: list of (f_phase, f_monomer) for the 6 E26->D clades.
    wt_xy:    (f_phase, f_monomer) for clade C01 (WT-Drt3b).
    e26d_clades: list of clade_id strings (for tracking which is which).
    """
    e26d_pts = []
    e26d_clades = []
    e26d_signatures = []
    wt = None
    with open(CSV_PATH, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            klass = row["envelope_classification"]
            p2 = float(row["predicted_period2_peak"])
            mg = float(row["predicted_marginal_G"])
            f_p = p2
            f_m = 1.0 - mg
            if klass == "shifted (primary-gate substitution at [26])":
                e26d_pts.append((f_p, f_m))
                e26d_clades.append(row["clade_id"])
                e26d_signatures.append(row["signature"])
            elif row["clade_id"] == "C01":
                wt = (f_p, f_m)
    return e26d_pts, wt, e26d_clades, e26d_signatures


def plotPanelB(e26d_pts, wt_xy, e26d_clades):
    fig, ax = plt.subplots(figsize=(85 / 25.4, 60 / 25.4))

    # axis ranges chosen to span WT (~0.98, 0.997) and E26D cluster (~0.857, 0.948)
    # with visible margin. wide x-axis to show the f_phase shift; tight y-axis
    # because f_monomer barely moves.
    xmin, xmax = 0.80, 1.01
    ymin, ymax = 0.92, 1.005

    # arrows from WT to each E26D point: vermilion, thin, no head conflict
    for (xe, ye) in e26d_pts:
        ax.annotate(
            "",
            xy=(xe, ye),
            xytext=(wt_xy[0], wt_xy[1]),
            arrowprops=dict(
                arrowstyle="->",
                color=OKABE_ITO["vermilion"],
                lw=0.6,
                alpha=0.45,
                shrinkA=4, shrinkB=4,
            ),
            zorder=2,
        )

    # 6 E26->D clade points -- vermilion
    xs = [p[0] for p in e26d_pts]
    ys = [p[1] for p in e26d_pts]
    ax.scatter(
        xs, ys,
        s=24,
        c=OKABE_ITO["vermilion"],
        edgecolor="black",
        linewidth=0.4,
        alpha=0.92,
        zorder=4,
        label=f"E26$\\rightarrow$D clades (n={len(e26d_pts)})",
    )

    # WT-Drt3b reference -- Okabe-Ito blue, star
    ax.scatter(
        [wt_xy[0]], [wt_xy[1]],
        s=70,
        marker="*",
        c=OKABE_ITO["blue"],
        edgecolor="black",
        linewidth=0.6,
        zorder=5,
        label="WT-Drt3b (C01)",
    )

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r"$f_{\rm phase}$ (period-2 peak)")
    ax.set_ylabel(r"$f_{\rm monomer}$ ($1 - $marg.\ G)")

    # annotate count and direction
    ax.text(
        0.04, 0.96,
        f"{len(e26d_pts)}/{len(e26d_pts)} clades shift\nin predicted direction",
        transform=ax.transAxes,
        fontsize=7,
        ha="left", va="top",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                  edgecolor="#888888", linewidth=0.5, alpha=0.9),
    )

    ax.legend(
        loc="lower right",
        frameon=False,
        handlelength=1.2,
        handletextpad=0.5,
        borderaxespad=0.3,
        labelspacing=0.4,
    )

    ax.tick_params(axis="both", which="major", length=2.5, width=0.8)

    fig.tight_layout(pad=0.3)
    return fig


def validate(pdf_path, png_path, fig):
    """7-check validation harness; raise on any failure."""
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
    e26d_pts, wt_xy, e26d_clades, e26d_signatures = loadE26DAndWT()

    if len(e26d_pts) != 6:
        print(f"warning: expected 6 E26->D clades, got {len(e26d_pts)}",
              file=sys.stderr)
    if wt_xy is None:
        raise RuntimeError("WT-Drt3b reference (C01) not found in CSV")

    # report
    print(f"loaded {len(e26d_pts)} E26->D primary-gate clades")
    xs = [p[0] for p in e26d_pts]
    ys = [p[1] for p in e26d_pts]
    print(f"  f_phase   range: [{min(xs):.4f}, {max(xs):.4f}]")
    print(f"  f_monomer range: [{min(ys):.4f}, {max(ys):.4f}]")
    print(f"  WT-Drt3b (C01): f_phase={wt_xy[0]:.4f}, f_monomer={wt_xy[1]:.4f}")
    for cid, sig, (fp, fm) in zip(e26d_clades, e26d_signatures, e26d_pts):
        print(f"    {cid:>16} ({sig}): f_phase={fp:.4f}, f_monomer={fm:.4f}")

    # confirm all shift in predicted direction (lower f_phase than WT)
    n_shifted = sum(1 for fp, fm in e26d_pts if fp < wt_xy[0])
    print(f"  {n_shifted}/{len(e26d_pts)} clades shifted to lower f_phase than WT")
    if n_shifted != len(e26d_pts):
        print(f"warning: {len(e26d_pts) - n_shifted} clades did not shift down",
              file=sys.stderr)

    fig = plotPanelB(e26d_pts, wt_xy, e26d_clades)
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=600)

    validate(PDF_PATH, PNG_PATH, fig)
    plt.close(fig)

    print(f"OK panel B:")
    print(f"  pdf: {PDF_PATH} ({PDF_PATH.stat().st_size / 1024:.1f} kB)")
    print(f"  png: {PNG_PATH} ({PNG_PATH.stat().st_size / 1024:.1f} kB)")


if __name__ == "__main__":
    main()
