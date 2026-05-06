"""fig 4 panel a: 33 secondary-residue clades fall in the WT-Drt3b envelope.

source: results/test_f_family_predictions_v1.csv
filter: envelope_classification == "in-envelope (secondary-only variation)"  (33 clades)

observable mapping (CSV columns -> framework f_phase, f_monomer):
    f_phase   <- predicted_period2_peak  (period-2 alternation strength;
                                          1.0 = perfect alternation between
                                          A and C in the dimer phase, 0.0 = no
                                          phasing). this is the cycle observable.
    f_monomer <- 1 - predicted_marginal_G (purity at canonical positions:
                                           predicted_marginal_G is the leak
                                           fraction of misincorporated G across
                                           the output, so 1 - G is the fraction
                                           consistent with the templated cycle).

WT-Drt3b reference (clade C01, signature ENHGRYTTKY) sits at
    f_phase ~ 0.980, f_monomer ~ 0.997  (period2_peak = 0.9805, marg_G = 0.0033).

envelope: f_phase > 0.95 AND f_monomer > 0.95 (per dispatch). drawn as a
green-tinted rectangle in the upper-right corner of the apparatus space.

reader takeaway: residue diversity at the 5 secondary positions (R168, Y170,
T335, T338, R408) does not push the apparatus signature outside the WT
envelope -- the cyclic-conformational mode is robust to that diversity.
"""
import csv
import sys
import subprocess
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
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
PDF_PATH = OUT_DIR / "fig4_panel_a.pdf"
PNG_PATH = OUT_DIR / "fig4_panel_a.png"

# envelope thresholds (per dispatch)
F_PHASE_THRESHOLD = 0.95
F_MONOMER_THRESHOLD = 0.95


def loadSecondaryClades():
    """return (f_phase_arr, f_monomer_arr, wt_xy) for the 33 in-envelope clades.

    wt_xy is the (f_phase, f_monomer) tuple of clade C01 (WT-Drt3b
    representative ENHGRYTTKY), used as the reference point.
    """
    f_phase = []
    f_monomer = []
    wt = None
    with open(CSV_PATH, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            klass = row["envelope_classification"]
            p2 = float(row["predicted_period2_peak"])
            mg = float(row["predicted_marginal_G"])
            f_p = p2
            f_m = 1.0 - mg
            if klass == "in-envelope (secondary-only variation)":
                f_phase.append(f_p)
                f_monomer.append(f_m)
                if row["clade_id"] == "C01":
                    wt = (f_p, f_m)
    return np.array(f_phase), np.array(f_monomer), wt


def plotPanelA(f_phase, f_monomer, wt_xy):
    fig, ax = plt.subplots(figsize=(85 / 25.4, 60 / 25.4))

    # axis ranges chosen to include envelope corner and all points with margin.
    # secondary clades cluster at (0.978-0.981, 0.996-0.997); WT envelope
    # threshold is 0.95; we show 0.93-1.0 on both axes so the envelope corner
    # is visible.
    xmin, xmax = 0.93, 1.001
    ymin, ymax = 0.93, 1.001

    # green-tinted envelope rectangle in the upper-right (Okabe-Ito green
    # at low alpha)
    env_rect = Rectangle(
        (F_PHASE_THRESHOLD, F_MONOMER_THRESHOLD),
        xmax - F_PHASE_THRESHOLD,
        ymax - F_MONOMER_THRESHOLD,
        facecolor=OKABE_ITO["green"],
        edgecolor=OKABE_ITO["green"],
        alpha=0.18,
        linewidth=0.8,
        zorder=1,
        label="WT envelope\n($f_{\\rm phase}>0.95$, $f_{\\rm monomer}>0.95$)",
    )
    ax.add_patch(env_rect)

    # 33 secondary-residue clades as dark gray points (slight jitter not needed;
    # the cluster is real and visible at this zoom)
    ax.scatter(
        f_phase, f_monomer,
        s=18,
        c="#404040",
        edgecolor="white",
        linewidth=0.4,
        alpha=0.85,
        zorder=3,
        label=f"secondary-residue clades (n={len(f_phase)})",
    )

    # WT-Drt3b reference in accent color (Okabe-Ito blue), larger marker
    if wt_xy is not None:
        ax.scatter(
            [wt_xy[0]], [wt_xy[1]],
            s=55,
            marker="*",
            c=OKABE_ITO["blue"],
            edgecolor="black",
            linewidth=0.6,
            zorder=4,
            label="WT-Drt3b (C01)",
        )

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r"$f_{\rm phase}$ (period-2 peak)")
    ax.set_ylabel(r"$f_{\rm monomer}$ ($1 - $marg.\ G)")

    # in-envelope count annotation in the upper-left of the plot
    n_in = int(np.sum((f_phase > F_PHASE_THRESHOLD) &
                      (f_monomer > F_MONOMER_THRESHOLD)))
    ax.text(
        0.04, 0.96,
        f"{n_in}/{len(f_phase)} clades in-envelope",
        transform=ax.transAxes,
        fontsize=7,
        ha="left", va="top",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                  edgecolor="#888888", linewidth=0.5, alpha=0.9),
    )

    ax.legend(
        loc="lower left",
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
    f_phase, f_monomer, wt_xy = loadSecondaryClades()

    if len(f_phase) != 33:
        print(f"warning: expected 33 secondary-residue clades, got {len(f_phase)}",
              file=sys.stderr)

    # report ranges
    print(f"loaded {len(f_phase)} secondary-residue clades")
    print(f"  f_phase   range: [{f_phase.min():.4f}, {f_phase.max():.4f}]")
    print(f"  f_monomer range: [{f_monomer.min():.4f}, {f_monomer.max():.4f}]")
    if wt_xy is not None:
        print(f"  WT-Drt3b (C01): f_phase={wt_xy[0]:.4f}, f_monomer={wt_xy[1]:.4f}")

    # check all in envelope
    n_in = int(np.sum((f_phase > F_PHASE_THRESHOLD) &
                      (f_monomer > F_MONOMER_THRESHOLD)))
    print(f"  {n_in}/{len(f_phase)} clades inside envelope "
          f"(f_phase > {F_PHASE_THRESHOLD}, f_monomer > {F_MONOMER_THRESHOLD})")
    if n_in != len(f_phase):
        print(f"warning: {len(f_phase) - n_in} clades fell outside the envelope",
              file=sys.stderr)

    fig = plotPanelA(f_phase, f_monomer, wt_xy)
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=600)

    validate(PDF_PATH, PNG_PATH, fig)
    plt.close(fig)

    print(f"OK panel A:")
    print(f"  pdf: {PDF_PATH} ({PDF_PATH.stat().st_size / 1024:.1f} kB)")
    print(f"  png: {PNG_PATH} ({PNG_PATH.stat().st_size / 1024:.1f} kB)")


if __name__ == "__main__":
    main()
