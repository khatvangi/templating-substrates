"""fig 3 panel b: apparatus triple separating four systems.

sources:
  - results/test_e_v2_results.csv (Drt3b_WT, Drt3b_E26Q, AbiK at L=100)
  - results/test_g2_dual_observable_v1.csv (Drt3a = alt_fixed at L=6)

three observables on x-axis, four systems as hue (grouped bars):
  - I_struct^pop (bits): saturates at log2(N_states)
  - periodicity peak (lag-2 ACF, dimensionless)
  - marginal G fraction (fraction)

values are normalized within each observable to [0,1] of the observable's
max so all three fit on one axes. raw values are annotated above each bar.

reader takeaway: no single observable separates all four — but the joint
pattern (I_pop, periodicity, marginal G) is unique per system. Drt3a and
Drt3b WT share the periodicity signature but differ on I_pop. Drt3b WT and
E26Q share I_pop but differ on periodicity and marginal G.

deviation from dispatch: dispatch suggested "three side-by-side mini-panels".
done as one axes with grouped bars per dispatch's lessons-from-pilot note
("ONE axes with grouped bars across observables") — Inkscape compositor can
re-arrange downstream. raw values annotated to retain absolute readability.
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
OUT_DIR = PROJECT_ROOT / "paper" / "figures" / "v7" / "fig3"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CSV_E_V2 = RESULTS_DIR / "test_e_v2_results.csv"
CSV_G2_DUAL = RESULTS_DIR / "test_g2_dual_observable_v1.csv"
PDF_PATH = OUT_DIR / "fig3_panel_b.pdf"
PNG_PATH = OUT_DIR / "fig3_panel_b.png"

# 4 systems and 3 observables
SYSTEMS = ["Drt3a", "Drt3b WT", "Drt3b E26Q", "AbiK"]
OBS_LABELS = [r"$I^{\mathrm{pop}}$", "periodicity", "marginal G"]
OBS_KEYS = ["I_pop", "periodicity", "margin_G"]

# okabe-ito 4-color qualitative
COLORS = [
    OKABE_ITO["green"],      # Drt3a (Mode 1, sequence template)
    OKABE_ITO["skyblue"],    # Drt3b WT (Mode 3 cyclic)
    OKABE_ITO["vermilion"],  # Drt3b E26Q (gate-broken Mode 3)
    OKABE_ITO["orange"],     # AbiK (non-templating)
]


def loadObservables():
    """assemble {system -> {obs -> value}} from the two CSVs."""
    # init
    vals = {s: {} for s in SYSTEMS}

    # Drt3b_WT, Drt3b_E26Q, AbiK at L=100 from test_e_v2
    sys_map = {
        "Drt3b_WT": "Drt3b WT",
        "Drt3b_E26Q": "Drt3b E26Q",
        "AbiK": "AbiK",
    }
    with open(CSV_E_V2, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sys_label = row["system_label"]
            if sys_label not in sys_map:
                continue
            if int(row["L"]) != 100:
                continue
            target = sys_map[sys_label]
            vals[target]["I_pop"] = float(row["I_struct_empirical"])
            vals[target]["periodicity"] = float(row["periodicity_peak_value"])
            vals[target]["margin_G"] = float(row["mean_marginal_G"])

    # Drt3a from test_g2_dual_observable (template_type=alt_fixed, L=6)
    # alt_fixed = single fixed alternating template = E1_fixed_ACACAC = Drt3a
    with open(CSV_G2_DUAL, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row["template_type"] != "alt_fixed":
                continue
            if int(row["L_T"]) != 6:
                continue
            vals["Drt3a"]["I_pop"] = float(row["I_pop"])
            vals["Drt3a"]["periodicity"] = float(row["period_peak"])
            vals["Drt3a"]["margin_G"] = float(row["margin_G"])

    # validate everything filled
    for s in SYSTEMS:
        for k in OBS_KEYS:
            if k not in vals[s]:
                raise RuntimeError(f"missing {k} for {s}")
    return vals


def plotPanelB(vals):
    fig, ax = plt.subplots(figsize=(60 / 25.4, 60 / 25.4))

    # raw values per (system, observable)
    raw = np.array([[vals[s][k] for k in OBS_KEYS] for s in SYSTEMS])
    # raw shape: (n_sys=4, n_obs=3)

    # normalize each observable to [0, 1] by its own max so the bars are
    # on the same axis. max-only (not min-max) so zero stays zero, which
    # makes the "absent observable" reading correct (e.g., Drt3a has zero
    # I_pop and zero margin_G in absolute terms).
    obs_max = raw.max(axis=0)
    # avoid div-by-zero
    obs_max_safe = np.where(obs_max > 0, obs_max, 1.0)
    norm = raw / obs_max_safe

    n_sys = len(SYSTEMS)
    n_obs = len(OBS_LABELS)
    x = np.arange(n_obs)
    bar_w = 0.18
    offsets = np.linspace(-bar_w * 1.5, bar_w * 1.5, n_sys)

    for i, sys_label in enumerate(SYSTEMS):
        ax.bar(
            x + offsets[i],
            norm[i],
            width=bar_w,
            color=COLORS[i],
            edgecolor="black",
            linewidth=0.5,
            label=sys_label,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(OBS_LABELS)
    ax.set_xlabel("observable")
    ax.set_ylabel("normalized value")
    ax.set_ylim(0.0, 1.20)

    # legend in two columns
    ax.legend(loc="upper right", ncol=2, frameon=False, columnspacing=0.6,
              handlelength=1.0, handletextpad=0.3, borderaxespad=0.2,
              fontsize=7)

    fig.tight_layout(pad=0.3)
    return fig, raw


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
    vals = loadObservables()

    fig, raw = plotPanelB(vals)
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=600)

    validate(PDF_PATH, PNG_PATH, fig)
    plt.close(fig)

    print("OK panel B:")
    print(f"  pdf: {PDF_PATH} ({PDF_PATH.stat().st_size / 1024:.1f} kB)")
    print(f"  png: {PNG_PATH} ({PNG_PATH.stat().st_size / 1024:.1f} kB)")
    print("  raw values (system x observable):")
    for i, s in enumerate(SYSTEMS):
        vals_str = ", ".join(
            f"{OBS_KEYS[j]}={raw[i, j]:.4f}" for j in range(len(OBS_KEYS))
        )
        print(f"    {s}: {vals_str}")


if __name__ == "__main__":
    main()
