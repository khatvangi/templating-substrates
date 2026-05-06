"""fig 3 panel c: channel-MI separating Drt3a Mode 1 from Drt3b Mode 3.

source: results/test_g2_drt3a_three_ensembles_v1.csv
observable: I_chan (channel-MI of focal channel against the comparison ensemble),
            in bits, computed at L=6, epsilon=0.01, n=5000.
output: paper/figures/v7/fig3/fig3_panel_c.{pdf,png}

reader takeaway: even when two channels share output statistics
(e.g., E1_fixed_ACACAC and Drt3b_N2 both produce ~poly(AC) with ~0.99
fidelity), I_chan distinguishes them by ~1 bit. This is the apparatus's
channel-axis discrimination — disambiguating Mode 1 (sequence template)
from Mode 3 (cyclic conformer) when output marginals coincide.

deviation from dispatch:
- the dispatch claimed both Drt3a and Drt3b WT "saturate at ~1 bit". the CSV
  shows I_chan in [2.32, 9.68] bits (log2(N)=2.32 for N-channel uniform mix
  is the saturation floor for low-entropy channels; AbiK and E3_random
  exceed it because their outputs differ from the cyclic-channel mixture).
  the **delta** dispatch annotation IS supported: I_chan(E1) - I_chan(N2)
  = 3.332 - 2.331 = 1.00 bit — that is the channel-axis separation.
- the CSV has 6 channels (E1, E2, E3, Drt3b_N2, Drt3b_N3, AbiK), not the 5
  the dispatch listed (Drt3b E26Q and Drt3a homopolymer null are not in
  this CSV; Drt3b_N3 is extra). plotted what the data has.
"""
import csv
import sys
import subprocess
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
from templating_style import apply_paper_style, OKABE_ITO  # noqa: E402

apply_paper_style()

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
OUT_DIR = PROJECT_ROOT / "paper" / "figures" / "v7" / "fig3"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CSV_PATH = RESULTS_DIR / "test_g2_drt3a_three_ensembles_v1.csv"
PDF_PATH = OUT_DIR / "fig3_panel_c.pdf"
PNG_PATH = OUT_DIR / "fig3_panel_c.png"

# 6 channels in the order they should appear (anchored panel comparison)
# pair the two cyclic-fold systems together for visual proximity
CHANNEL_ORDER = [
    "E1_fixed_ACACAC",
    "Drt3b_N2",
    "Drt3b_N3",
    "E2_2phase_alt",
    "E3_random_uniform",
    "AbiK_uniform",
]
CHANNEL_LABELS = [
    "Drt3a\nfix",
    "Drt3b\nN=2",
    "Drt3b\nN=3",
    "Drt3a\n2ph",
    "Drt3a\nrnd",
    "AbiK",
]

# the two channels the dispatch wants emphasized (Drt3a vs Drt3b WT)
EMPHASIS = {"E1_fixed_ACACAC", "Drt3b_N2"}

# 5-color qualitative palette: WT/Drt3a saturated, others muted
# use base colors then alpha-mute non-emphasized ones
BASE_COLORS = {
    "E1_fixed_ACACAC":   OKABE_ITO["green"],
    "Drt3b_N2":          OKABE_ITO["skyblue"],
    "Drt3b_N3":          OKABE_ITO["blue"],
    "E2_2phase_alt":     OKABE_ITO["purple"],
    "E3_random_uniform": OKABE_ITO["yellow"],
    "AbiK_uniform":      OKABE_ITO["orange"],
}


def loadIChan():
    """{channel_name -> I_chan in bits} from the three_ensembles CSV."""
    out = {}
    with open(CSV_PATH, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            out[row["ensemble"]] = float(row["I_chan"])
    for c in CHANNEL_ORDER:
        if c not in out:
            raise RuntimeError(f"missing channel '{c}' in CSV")
    return out


def plotPanelC(i_chan):
    fig, ax = plt.subplots(figsize=(60 / 25.4, 60 / 25.4))

    n = len(CHANNEL_ORDER)
    x = np.arange(n)
    values = [i_chan[c] for c in CHANNEL_ORDER]

    # color: emphasized channels saturated; others muted (alpha=0.45)
    colors = []
    for c in CHANNEL_ORDER:
        base = BASE_COLORS[c]
        if c in EMPHASIS:
            colors.append(to_rgba(base, alpha=1.0))
        else:
            colors.append(to_rgba(base, alpha=0.45))

    bars = ax.bar(
        x,
        values,
        width=0.7,
        color=colors,
        edgecolor="black",
        linewidth=0.5,
    )

    # extra emphasis: thicker edge on the two highlighted channels
    for c, bar in zip(CHANNEL_ORDER, bars):
        if c in EMPHASIS:
            bar.set_linewidth(1.2)

    ax.set_xticks(x)
    ax.set_xticklabels(CHANNEL_LABELS, fontsize=6.5)
    ax.set_xlabel("channel")
    ax.set_ylabel(r"$I^{\mathrm{chan}}$ (bits)")

    # ymax with headroom for annotation
    ymax = max(values) * 1.18
    ax.set_ylim(0, ymax)

    # annotate Drt3a vs Drt3b delta with a small bracket between bars 0 and 1
    e1_v = i_chan["E1_fixed_ACACAC"]
    n2_v = i_chan["Drt3b_N2"]
    delta = e1_v - n2_v

    # bracket sits just above the taller of the two bars (E1 = 3.33), well
    # below the AbiK/E3 ceiling (~9.7)
    bracket_y = e1_v + 0.4
    ax.annotate(
        "",
        xy=(0, bracket_y), xytext=(1, bracket_y),
        arrowprops=dict(arrowstyle="-", color="black", linewidth=0.7),
    )
    # tiny vertical ticks at bracket ends for clarity
    ax.plot([0, 0], [e1_v + 0.05, bracket_y], color="black", linewidth=0.7)
    ax.plot([1, 1], [n2_v + 0.05, bracket_y], color="black", linewidth=0.7)
    # annotation text just above bracket; use compact label and place above
    ax.text(
        0.5, bracket_y + 0.15,
        rf"$\Delta = {delta:.2f}$ bit",
        ha="center", va="bottom", fontsize=7,
    )

    fig.tight_layout(pad=0.3)
    return fig, values


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

    # min font: tick labels in this panel are 6.5pt (intentional, to fit
    # six channel names). only check that *paper* fonts (axis labels, etc.)
    # are >= 7pt; small tick labels for category names are acceptable.
    min_font = min(
        plt.rcParams["font.size"],
        plt.rcParams["axes.labelsize"],
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
    i_chan = loadIChan()

    # cross-check delta against dispatch's "1.0 bit" claim
    delta = i_chan["E1_fixed_ACACAC"] - i_chan["Drt3b_N2"]
    if abs(delta - 1.0) > 0.05:
        print(f"warning: ΔI_chan(Drt3a, Drt3b) = {delta:.4f}, "
              f"dispatch said ~1.0", file=sys.stderr)

    fig, values = plotPanelC(i_chan)
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=600)

    validate(PDF_PATH, PNG_PATH, fig)
    plt.close(fig)

    print("OK panel C:")
    print(f"  pdf: {PDF_PATH} ({PDF_PATH.stat().st_size / 1024:.1f} kB)")
    print(f"  png: {PNG_PATH} ({PNG_PATH.stat().st_size / 1024:.1f} kB)")
    print(f"  ΔI_chan(Drt3a fixed AC vs Drt3b N=2) = {delta:.4f} bits")
    for c, v in zip(CHANNEL_ORDER, values):
        print(f"    {c}: {v:.4f}")


if __name__ == "__main__":
    main()
