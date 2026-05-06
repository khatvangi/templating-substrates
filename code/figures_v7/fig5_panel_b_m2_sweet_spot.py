"""fig 5 panel b: M2 sweet-spot regime invariance.

source:
  - test_h5_r_sweep_v1.csv          (main, beta=10, L=32, 14 r values)
  - test_h5_beta_sensitivity_v1.csv (5 r values across beta in {2,5,10,20})
  - test_h5_L_sensitivity_v1.csv    (5 r values across L in {16,32,64,128})
output: paper/figures/v7/fig5/fig5_panel_b.{pdf,png}

reader takeaway: the M2 plateau peaks at r ~ 0.10 across (beta, L) regimes —
the cyclic active-site (R3) bound is set by mechanism, not parameters.
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
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
from templating_style import apply_paper_style, OKABE_ITO  # noqa: E402

apply_paper_style()

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
OUT_DIR = PROJECT_ROOT / "paper" / "figures" / "v7" / "fig5"
OUT_DIR.mkdir(parents=True, exist_ok=True)

R_SWEEP_CSV = RESULTS_DIR / "test_h5_r_sweep_v1.csv"
BETA_CSV = RESULTS_DIR / "test_h5_beta_sensitivity_v1.csv"
L_CSV = RESULTS_DIR / "test_h5_L_sensitivity_v1.csv"

PDF_PATH = OUT_DIR / "fig5_panel_b.pdf"
PNG_PATH = OUT_DIR / "fig5_panel_b.png"


def loadRSweep():
    """main curve: beta=10, L=32, full r grid."""
    rs, fits = [], []
    with open(R_SWEEP_CSV, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rs.append(float(row["r"]))
            fits.append(float(row["mean_final_fitness"]))
    return np.array(rs), np.array(fits)


def loadBetaSweep():
    """{beta -> (rs, fits)} at L=32."""
    out = defaultdict(lambda: ([], []))
    with open(BETA_CSV, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            beta = float(row["beta"])
            out[beta][0].append(float(row["r"]))
            out[beta][1].append(float(row["mean_final_fitness"]))
    return {b: (np.array(rs), np.array(fs)) for b, (rs, fs) in out.items()}


def loadLSweep():
    """{L -> (rs, fits)} at beta=10."""
    out = defaultdict(lambda: ([], []))
    with open(L_CSV, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            L = int(row["L"])
            out[L][0].append(float(row["r"]))
            out[L][1].append(float(row["mean_final_fitness"]))
    return {L: (np.array(rs), np.array(fs)) for L, (rs, fs) in out.items()}


def plotPanelB(r_main, fit_main, beta_data, L_data):
    fig, ax = plt.subplots(figsize=(85 / 25.4, 60 / 25.4))

    # six curves arranged on a viridis sequence; the main curve is the "fiducial"
    # and is drawn last so it sits on top.
    # ordering: low-extremes -> fiducial -> high-extremes for a sensible legend
    curves = [
        # (label, rs, fits)
        (r"$\beta$=2, L=32",  beta_data[2.0][0],  beta_data[2.0][1]),
        (r"$\beta$=20, L=32", beta_data[20.0][0], beta_data[20.0][1]),
        (r"$\beta$=10, L=16", L_data[16][0],      L_data[16][1]),
        (r"$\beta$=10, L=64", L_data[64][0],      L_data[64][1]),
        (r"$\beta$=10, L=128",L_data[128][0],     L_data[128][1]),
        (r"$\beta$=10, L=32 (main)", r_main,      fit_main),
    ]

    cmap = plt.cm.viridis
    n = len(curves)
    # avoid the brightest yellow at the very end (poor contrast on white)
    colors = [cmap(0.10 + 0.75 * i / (n - 1)) for i in range(n)]

    for i, (label, rs, fs) in enumerate(curves):
        # sort by r so log plot draws cleanly
        order = np.argsort(rs)
        rs, fs = rs[order], fs[order]
        is_main = (i == n - 1)
        ax.plot(
            rs, fs,
            marker="o", markersize=3 if not is_main else 4,
            linewidth=1.0 if not is_main else 1.6,
            color=colors[i],
            label=label,
            zorder=2 + i,
        )

    # vertical line at r* = 0.10
    ax.axvline(0.10, color="grey", linestyle="--", linewidth=0.8, zorder=0)

    ax.set_xscale("log")
    ax.set_xlabel(r"re-draw rate $r$")
    ax.set_ylabel("M2 plateau fitness")
    ax.set_xlim(0.008, 1.1)
    ax.set_ylim(0.20, 0.75)

    # annotate sweet spot
    ax.text(0.105, 0.21, r"$r^* = 0.10$", fontsize=7, color="grey",
            ha="left", va="bottom")

    ax.legend(loc="lower left", frameon=False, ncol=1,
              handlelength=1.4, handletextpad=0.4, borderaxespad=0.3,
              labelspacing=0.2)

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
    r_main, fit_main = loadRSweep()
    beta_data = loadBetaSweep()
    L_data = loadLSweep()

    # sanity: peak of main curve is at r=0.10
    peak_idx = int(np.argmax(fit_main))
    if abs(r_main[peak_idx] - 0.10) > 0.001:
        print(f"warning: main-curve peak is at r={r_main[peak_idx]} not 0.10",
              file=sys.stderr)

    # check peak invariance across regimes (within +/- 0.05 r range)
    for label, (rs, fs) in [
        ("beta=2", beta_data[2.0]), ("beta=5", beta_data[5.0]),
        ("beta=10", beta_data[10.0]), ("beta=20", beta_data[20.0]),
        ("L=16", L_data[16]), ("L=32", L_data[32]),
        ("L=64", L_data[64]), ("L=128", L_data[128]),
    ]:
        order = np.argsort(rs)
        rs_s, fs_s = rs[order], fs[order]
        peak_r = rs_s[int(np.argmax(fs_s))]
        if abs(peak_r - 0.10) > 0.06:
            print(f"warning: {label} peaks at r={peak_r}, not near 0.10",
                  file=sys.stderr)

    fig = plotPanelB(r_main, fit_main, beta_data, L_data)
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=600)
    validate(PDF_PATH, PNG_PATH, fig)
    plt.close(fig)
    print(f"OK panel B:")
    print(f"  pdf: {PDF_PATH} ({PDF_PATH.stat().st_size / 1024:.1f} kB)")
    print(f"  png: {PNG_PATH} ({PNG_PATH.stat().st_size / 1024:.1f} kB)")


if __name__ == "__main__":
    main()
