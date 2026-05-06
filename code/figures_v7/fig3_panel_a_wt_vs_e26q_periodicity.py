"""fig 3 panel a: Drt3b WT vs E26Q periodicity peak across L.

source: results/test_e_v2_results.csv
filter: system_label in {Drt3b_WT, Drt3b_E26Q}, L in {4, 10, 20, 100, 500}
observable: periodicity_peak_value (autocorrelation peak at the dominant lag)
output: paper/figures/v7/fig3/fig3_panel_a.{pdf,png}

reader takeaway: WT periodicity peak ~0.98 (Mode 3 cyclic, perfect alternation
in product), E26Q peak ~0.83 (gate-broken — dG misincorporation occasionally
breaks the alternation). The gap is stable across L, so the periodicity
observable is a robust apparatus signature for distinguishing intact vs
broken cyclic-conformer templating.

deviation from dispatch: the dispatch suggested an autocorrelation overlay
(lag 0..6), but only the peak value at the dominant lag is in the CSV
(no full ACF curve). using paired bars across L preserves data fidelity and
matches the "paired bar plot" option in the dispatch.
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
OUT_DIR = PROJECT_ROOT / "paper" / "figures" / "v7" / "fig3"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CSV_PATH = RESULTS_DIR / "test_e_v2_results.csv"
PDF_PATH = OUT_DIR / "fig3_panel_a.pdf"
PNG_PATH = OUT_DIR / "fig3_panel_a.png"

# L values to plot (skip L=2 — periodicity not yet established at this length)
L_VALUES = [4, 10, 20, 100, 500]
SYSTEMS = ["Drt3b_WT", "Drt3b_E26Q"]
SYSTEM_LABELS = ["Drt3b WT", "Drt3b E26Q"]

# okabe-ito two-color per dispatch: sky-blue WT, vermilion E26Q
COLORS = [OKABE_ITO["skyblue"], OKABE_ITO["vermilion"]]


def loadPeriodicity():
    """load periodicity_peak_value for both systems across L; returns
    {system_label -> {L -> peak_value}}."""
    data = {s: {} for s in SYSTEMS}
    with open(CSV_PATH, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sys_label = row["system_label"]
            if sys_label not in SYSTEMS:
                continue
            L = int(row["L"])
            if L not in L_VALUES:
                continue
            data[sys_label][L] = float(row["periodicity_peak_value"])
    # sanity check
    for s in SYSTEMS:
        for L in L_VALUES:
            if L not in data[s]:
                raise RuntimeError(f"missing periodicity for system={s} L={L}")
    return data


def plotPanelA(data):
    fig, ax = plt.subplots(figsize=(60 / 25.4, 60 / 25.4))

    n_sys = len(SYSTEMS)
    n_L = len(L_VALUES)
    x = np.arange(n_L)
    bar_w = 0.38
    offsets = np.array([-bar_w / 2, bar_w / 2])

    for i, sys_label in enumerate(SYSTEMS):
        peaks = [data[sys_label][L] for L in L_VALUES]
        ax.bar(
            x + offsets[i],
            peaks,
            width=bar_w,
            color=COLORS[i],
            edgecolor="black",
            linewidth=0.5,
            label=SYSTEM_LABELS[i],
        )

    # annotation: WT vs E26Q at L=100 (representative long-output case)
    wt_100 = data["Drt3b_WT"][100]
    e26q_100 = data["Drt3b_E26Q"][100]
    ax.axhline(wt_100, color=COLORS[0], linestyle=":", linewidth=0.6,
               alpha=0.7, zorder=0)
    ax.axhline(e26q_100, color=COLORS[1], linestyle=":", linewidth=0.6,
               alpha=0.7, zorder=0)

    ax.set_xticks(x)
    ax.set_xticklabels([str(L) for L in L_VALUES])
    ax.set_xlabel("output length L")
    ax.set_ylabel("periodicity peak (lag-2 ACF)")
    ax.set_ylim(0.0, 1.05)

    # text annotations: WT alternation: 0.98 / E26Q: 0.83 (matches dispatch)
    ax.text(0.02, 0.96, f"WT: {wt_100:.2f}",
            transform=ax.transAxes, fontsize=7, color=COLORS[0],
            ha="left", va="top", weight="bold")
    ax.text(0.02, 0.86, f"E26Q: {e26q_100:.2f}",
            transform=ax.transAxes, fontsize=7, color=COLORS[1],
            ha="left", va="top", weight="bold")

    ax.legend(loc="lower right", frameon=False, handlelength=1.0,
              handletextpad=0.4, borderaxespad=0.3)

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
    data = loadPeriodicity()

    # cross-check against dispatch: WT peak ~0.98, E26Q ~0.83 at L>=4
    wt_100 = data["Drt3b_WT"][100]
    e26q_100 = data["Drt3b_E26Q"][100]
    if abs(wt_100 - 0.98) > 0.02:
        print(f"warning: WT L=100 peak={wt_100:.4f} expected ~0.98", file=sys.stderr)
    if abs(e26q_100 - 0.83) > 0.02:
        print(f"warning: E26Q L=100 peak={e26q_100:.4f} expected ~0.83",
              file=sys.stderr)

    fig = plotPanelA(data)
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=600)

    validate(PDF_PATH, PNG_PATH, fig)
    plt.close(fig)

    print("OK panel A:")
    print(f"  pdf: {PDF_PATH} ({PDF_PATH.stat().st_size / 1024:.1f} kB)")
    print(f"  png: {PNG_PATH} ({PNG_PATH.stat().st_size / 1024:.1f} kB)")
    print(f"  WT L=100 peak: {wt_100:.4f}")
    print(f"  E26Q L=100 peak: {e26q_100:.4f}")


if __name__ == "__main__":
    main()
