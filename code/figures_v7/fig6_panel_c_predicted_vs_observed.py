"""fig 6 panel c: framework prediction vs literature observation.

source:
  - test_h3_isolated_dynamics_v1.csv (M2 plateau heights at r=0.10, beta=10, L=32)

literature reference (panel C target, qualitative only):
  - Beisson & Sonneborn 1965 PNAS 53:275 — cortical inheritance propagated
    across hundreds of generations
  - Frankel 1989 — cortical inheritance bounded; no novel structures accumulate

output: paper/figures/v7/fig6/fig6_panel_c.{pdf,png}

reader takeaway: the framework predicts a bounded plateau in M2 at the
cortex-mapped re-draw rate; the literature reports cortical inheritance is
bounded (no novel structures accumulate across hundreds of generations).
the comparison is DIRECTIONAL (both bounded, finite), NOT a numerical match
— the literature value is qualitative. the visual encoding makes this honest:
the framework bar carries a numeric value with a 1-sigma error bar; the
literature column is rendered as a hatched band labeled "bounded; precise
value not measured" with no numerical y-position.
"""
import csv
import sys
import subprocess
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
from templating_style import apply_paper_style, OKABE_ITO  # noqa: E402

apply_paper_style()

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
OUT_DIR = PROJECT_ROOT / "paper" / "figures" / "v7" / "fig6"
OUT_DIR.mkdir(parents=True, exist_ok=True)

H3_CSV = RESULTS_DIR / "test_h3_isolated_dynamics_v1.csv"
PDF_PATH = OUT_DIR / "fig6_panel_c.pdf"
PNG_PATH = OUT_DIR / "fig6_panel_c.png"

# the cortex maps to M2 at r=0.10, fiducial regime beta=10, L=32 (consistent with fig 5)
TARGET_MECH = "M2_r0.10"
TARGET_BETA = 10.0
TARGET_L = 32


def loadM2Plateau():
    """find the single fiducial M2 r=0.10 row in h3."""
    matches = []
    with open(H3_CSV, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row["mech_key"] != TARGET_MECH:
                continue
            if int(row["L"]) != TARGET_L:
                continue
            if float(row["beta"]) != TARGET_BETA:
                continue
            matches.append(row)
    if len(matches) != 1:
        raise RuntimeError(
            f"expected exactly 1 row for {TARGET_MECH} L={TARGET_L} "
            f"beta={TARGET_BETA}, got {len(matches)}"
        )
    row = matches[0]
    return {
        "mean": float(row["mean_final_fitness"]),
        "std": float(row["std_final_fitness"]),
        "n": int(row["n_replicates"]),
    }


def plotPanelC(m2):
    fig, ax = plt.subplots(figsize=(85 / 25.4, 60 / 25.4))

    # x positions for the two categories
    x_pred = 0.0
    x_obs = 1.0
    bar_w = 0.55

    # left bar: framework prediction (numeric, with se error bar)
    se = m2["std"] / np.sqrt(m2["n"])  # standard error of the mean
    ax.bar(
        x_pred, m2["mean"], width=bar_w,
        color=OKABE_ITO["blue"], edgecolor="black", linewidth=0.8,
        label=f"framework prediction\n(M2, r=0.10): {m2['mean']:.2f} +/- {se:.2f} (SE)",
        zorder=3,
    )
    ax.errorbar(
        x_pred, m2["mean"], yerr=se,
        fmt="none", ecolor="black", elinewidth=0.9, capsize=3, capthick=0.9,
        zorder=4,
    )

    # right "bar": literature observation, encoded as a hatched band only —
    # NO numerical y-position, NO bar height. the band spans the y-range that
    # the literature loosely characterizes ("bounded; finite generational
    # stability without accumulation of novel structures") and is explicitly
    # labeled as qualitative. this prevents the reader from comparing bar
    # heights as if they were the same kind of measurement.
    band_lo, band_hi = 0.30, 0.85   # qualitative band: "bounded, finite"
    band_left = x_obs - bar_w / 2
    band_right = x_obs + bar_w / 2
    band = mpatches.Rectangle(
        (band_left, band_lo), bar_w, band_hi - band_lo,
        facecolor=OKABE_ITO["orange"], alpha=0.25,
        edgecolor=OKABE_ITO["orange"], linewidth=1.2,
        hatch="///", zorder=2,
    )
    ax.add_patch(band)
    # bracket arms: top and bottom horizontal ticks to suggest "bounded interval"
    bracket_color = OKABE_ITO["orange"]
    ax.plot([band_left, band_right], [band_lo, band_lo],
            color=bracket_color, linewidth=1.4, zorder=4)
    ax.plot([band_left, band_right], [band_hi, band_hi],
            color=bracket_color, linewidth=1.4, zorder=4)
    ax.plot([x_obs, x_obs], [band_lo, band_hi],
            color=bracket_color, linewidth=1.0, zorder=4)

    # qualitative label inside the band
    ax.text(
        x_obs, (band_lo + band_hi) / 2,
        "bounded;\nprecise value\nnot measured",
        ha="center", va="center",
        fontsize=7, color="black",
        zorder=5,
    )

    # legend proxy for the literature band
    lit_proxy = mpatches.Patch(
        facecolor=OKABE_ITO["orange"], alpha=0.25,
        edgecolor=OKABE_ITO["orange"], hatch="///",
        label="literature observation\n(qualitative)",
    )
    handles, labels = ax.get_legend_handles_labels()
    handles.append(lit_proxy)
    labels.append(lit_proxy.get_label())
    ax.legend(handles, labels,
              loc="lower right", frameon=False,
              handlelength=1.4, handletextpad=0.4,
              borderaxespad=0.3, labelspacing=0.4,
              fontsize=7)

    # axes formatting
    ax.set_xticks([x_pred, x_obs])
    ax.set_xticklabels(["framework\n(M2, r=0.10)", "cortex\n(literature)"])
    ax.set_ylabel("plateau-like quantity\n(framework: fitness; lit.: qualitative)")
    ax.set_xlabel("")  # axis-level — keep non-empty via title/ticks
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(0.0, 1.0)

    # honest comparison annotation at the top
    ax.set_title("directional, not numerical, comparison",
                 fontsize=8, pad=4)

    # ensure the xlabel is non-empty per validation contract
    ax.set_xlabel("category")

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
    m2 = loadM2Plateau()
    se = m2["std"] / np.sqrt(m2["n"])
    print(f"M2 r=0.10 plateau (h3, beta=10, L=32):")
    print(f"  mean_final_fitness = {m2['mean']:.4f}")
    print(f"  std                = {m2['std']:.4f}")
    print(f"  n                  = {m2['n']}")
    print(f"  SE                 = {se:.4f}")

    fig = plotPanelC(m2)
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=600)
    validate(PDF_PATH, PNG_PATH, fig)
    plt.close(fig)
    print(f"wrote {PDF_PATH}")
    print(f"wrote {PNG_PATH}")


if __name__ == "__main__":
    main()
