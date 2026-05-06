"""fig 1 panel C: comparison-ensemble robustness of I_struct^chan (test G3).

source: results/test_g3_classification_stability_v1.csv
columns: ensemble_id, focal_channel, ensemble_size, ensemble, I_chan_focal,
         period_lag, period_peak, is_mode3

violin/strip plot: one distribution per ensemble (ensemble_id 0..9 = ten
random 4-channel ensembles); the y-value is I_chan_focal across the four
focal channels in that ensemble. canonical 5-channel reference (G2's
{Drt3a-WC, Drt3b-N=2, Drt3b-N=3, Drt3b-N=5, AbiK-uniform}) is shown as a
horizontal reference line at I_chan(Drt3a-WC) = 2.3066 bits (from
test_g3_g2_result_robustness_v1.csv, E5_G2_canonical row, I_chan_Drt3a_WC_2phase).

NOTE on dispatch deviation: dispatch refers to "1 = canonical, 2-11 =
alternative draws" (11 ensembles). actual data has 10 random ensembles
in classification_stability_v1.csv plus a separate canonical row in
g2_result_robustness_v1.csv. we plot the 10 random as violins and overlay
the canonical as a single reference line. classification stability
(is_mode3 verdicts) is preserved across all ensembles per G3 README.

output: paper/figures/v7/fig1/fig1_panel_c.{pdf,png}

reader takeaway: I_chan numerical value depends on ensemble (P_G3_2 expected),
but the qualitative Mode 1 vs Mode 3 classification is preserved across all
alternative ensemble draws.
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
OUT_DIR = PROJECT_ROOT / "paper" / "figures" / "v7" / "fig1"
OUT_DIR.mkdir(parents=True, exist_ok=True)

STAB_CSV = RESULTS_DIR / "test_g3_classification_stability_v1.csv"
CANON_CSV = RESULTS_DIR / "test_g3_g2_result_robustness_v1.csv"
PDF_PATH = OUT_DIR / "fig1_panel_c.pdf"
PNG_PATH = OUT_DIR / "fig1_panel_c.png"

NEUTRAL_COLOR = OKABE_ITO["skyblue"]
ACCENT_COLOR = OKABE_ITO["vermilion"]


def loadStability():
    """load classification_stability; return list of dicts and counts."""
    rows = []
    with open(STAB_CSV, newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            rows.append({
                "ensemble_id": int(r["ensemble_id"]),
                "focal": r["focal_channel"],
                "I_chan": float(r["I_chan_focal"]),
                "is_mode3": r["is_mode3"].strip().lower() == "true",
            })
    return rows


def loadCanonical():
    """load G2 canonical (E5) reference value; return float (bits)."""
    with open(CANON_CSV, newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            if r["ensemble_name"] == "E5_G2_canonical":
                # use Drt3a-WC 2-phase (matches G2's headline value)
                return float(r["I_chan_Drt3a_WC_2phase"])
    raise RuntimeError("E5_G2_canonical row not found")


def plotPanelC(rows, canonical_value):
    fig, ax = plt.subplots(figsize=(60 / 25.4, 60 / 25.4))

    # group I_chan_focal values by ensemble_id
    ens_ids = sorted({r["ensemble_id"] for r in rows})
    data_per_ens = [
        [r["I_chan"] for r in rows if r["ensemble_id"] == eid]
        for eid in ens_ids
    ]

    # x positions for the 10 alternative ensembles (label them 1..10)
    x_alt = np.arange(1, len(ens_ids) + 1)

    # violin plot (matplotlib violinplot since we want minimal seaborn complexity)
    parts = ax.violinplot(
        data_per_ens,
        positions=x_alt,
        widths=0.7,
        showmeans=False,
        showextrema=False,
        showmedians=False,
    )
    for body in parts["bodies"]:
        body.set_facecolor(NEUTRAL_COLOR)
        body.set_edgecolor("black")
        body.set_linewidth(0.5)
        body.set_alpha(0.65)

    # individual points overlaid (with mild jitter)
    rng = np.random.default_rng(42)
    for i, eid in enumerate(ens_ids):
        vals = data_per_ens[i]
        jitter = rng.uniform(-0.10, 0.10, size=len(vals))
        ax.scatter(
            x_alt[i] + jitter,
            vals,
            s=8,
            color="black",
            alpha=0.75,
            zorder=4,
            linewidths=0,
        )

    # canonical reference line
    ax.axhline(
        canonical_value,
        color=ACCENT_COLOR,
        linestyle="--",
        linewidth=1.2,
        zorder=2,
        label=f"canonical (E5)\n= {canonical_value:.2f} bits",
    )

    ax.set_xlim(0.3, len(ens_ids) + 0.7)
    ax.set_xticks(x_alt)
    ax.set_xticklabels([str(i) for i in x_alt])
    ax.set_xlabel("alternative ensemble (random 4-channel draws)")
    ax.set_ylabel(r"$I_{\mathrm{struct}}^{\mathrm{chan}}$ (bits)")

    ax.legend(
        loc="lower right",
        frameon=False,
        handletextpad=0.4,
        borderaxespad=0.3,
        fontsize=6,
    )

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
    rows = loadStability()
    if not rows:
        print(f"ERROR: no rows in {STAB_CSV}", file=sys.stderr)
        sys.exit(1)

    canonical_value = loadCanonical()

    # diagnostic: print summary
    ens_ids = sorted({r["ensemble_id"] for r in rows})
    print(f"G3 stability: {len(rows)} rows across {len(ens_ids)} ensembles")
    print(f"G2 canonical (E5) Drt3a-WC 2phase I_chan = {canonical_value:.4f} bits")

    # report classification stability per focal channel
    by_focal = {}
    for r in rows:
        by_focal.setdefault(r["focal"], []).append(r["is_mode3"])
    print("\nclassification stability per focal channel:")
    n_unstable = 0
    for focal, verdicts in sorted(by_focal.items()):
        unique = set(verdicts)
        status = "STABLE" if len(unique) == 1 else "UNSTABLE"
        if len(unique) > 1:
            n_unstable += 1
        print(f"  {focal:<22s} {status} (n={len(verdicts)}, "
              f"verdicts={list(unique)})")
    print(f"\nstable channels: {len(by_focal) - n_unstable} / {len(by_focal)}")

    fig = plotPanelC(rows, canonical_value)
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=600)

    validate(PDF_PATH, PNG_PATH, fig)
    plt.close(fig)

    print("\nOK panel C:")
    print(f"  pdf: {PDF_PATH} ({PDF_PATH.stat().st_size / 1024:.1f} kB)")
    print(f"  png: {PNG_PATH} ({PNG_PATH.stat().st_size / 1024:.1f} kB)")


if __name__ == "__main__":
    main()
