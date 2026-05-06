"""fig 2 panel a: Mode 3 (cyclic active site) saturation at log_2(N).

source: results/test_b_results.csv
filter: epsilon=0.001 (low-noise lines), N in {2,3,4,5,6,8,10}, L in {N,5N,10N,50N}
output: paper/figures/v7/fig2/fig2_panel_a.{pdf,png}

reader takeaway: I_struct(X;Y) for Mode 3 saturates at log_2(N) regardless
of L, because the cyclic active site has only N distinguishable phase states.
Curves at different L collapse onto the same N-dependent ceiling.
"""
import csv
import sys
import subprocess
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from PIL import Image

# style module lives in same directory
sys.path.insert(0, str(Path(__file__).resolve().parent))
from templating_style import apply_paper_style, OKABE_ITO  # noqa: E402

apply_paper_style()

# paths
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
OUT_DIR = PROJECT_ROOT / "paper" / "figures" / "v7" / "fig2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CSV_PATH = RESULTS_DIR / "test_b_results.csv"
PDF_PATH = OUT_DIR / "fig2_panel_a.pdf"
PNG_PATH = OUT_DIR / "fig2_panel_a.png"

# dispatch: L in {N, 5N, 10N, 50N}; data has {N, 2N, 5N, 10N, 50N} so we
# pick the four requested multipliers
L_MULTIPLIERS = [1, 5, 10, 50]
L_LABELS = [r"$L=N$", r"$L=5N$", r"$L=10N$", r"$L=50N$"]
N_VALUES = [2, 3, 4, 5, 6, 8, 10]
EPS_FIXED = 0.001  # low-noise so the saturation regime is cleanly visible


def loadB():
    """load (N, L, I) tuples for epsilon=EPS_FIXED into a nested dict
    {multiplier -> {N -> I_empirical}}."""
    series = {m: {} for m in L_MULTIPLIERS}
    with open(CSV_PATH, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            eps = float(row["epsilon"])
            if abs(eps - EPS_FIXED) > 1e-9:
                continue
            N = int(row["N"])
            L = int(row["L"])
            if N not in N_VALUES:
                continue
            # accept only exact multiples in our requested set
            if L % N != 0:
                continue
            mult = L // N
            if mult not in L_MULTIPLIERS:
                continue
            series[mult][N] = float(row["I_empirical"])
    # sanity: every (mult, N) must be present
    for m in L_MULTIPLIERS:
        for N in N_VALUES:
            if N not in series[m]:
                raise RuntimeError(
                    f"missing data for L={m}*N at N={N} eps={EPS_FIXED}")
    return series


def plotPanelA(series):
    fig, ax = plt.subplots(figsize=(85 / 25.4, 60 / 25.4))

    # viridis sequential — darker = larger L
    viridis = cm.get_cmap("viridis")
    colors = [viridis(t) for t in np.linspace(0.15, 0.85, len(L_MULTIPLIERS))]
    # different markers per L so the four collapsed lines remain visible
    markers = ["o", "s", "^", "D"]

    Ns = np.array(N_VALUES, dtype=float)

    # one line per L multiplier — they collapse onto log_2(N) (the panel's
    # take-home: saturation independent of L). Offset markers by tiny x-jitter
    # so the four overlapping markers can be told apart.
    n_lines = len(L_MULTIPLIERS)
    log_jitter = 0.012  # multiplicative jitter on log-scale x
    for i, mult in enumerate(L_MULTIPLIERS):
        ys = [series[mult][N] for N in N_VALUES]
        # symmetric jitter spread around 1.0 in log space
        jit = np.exp((i - (n_lines - 1) / 2.0) * log_jitter)
        ax.plot(Ns * jit, ys, color=colors[i], marker=markers[i],
                markersize=3.5, linewidth=0.8, alpha=0.9,
                label=L_LABELS[i], zorder=3 + i)

    # theoretical saturation: log_2(N), drawn ON TOP so it stays visible
    n_dense = np.linspace(min(Ns), max(Ns), 200)
    ax.plot(n_dense, np.log2(n_dense), color="black",
            linestyle="--", linewidth=1.2, label=r"$\log_2 N$ (theory)",
            zorder=10)

    ax.set_xscale("log")
    ax.set_xticks(N_VALUES)
    ax.set_xticklabels([str(N) for N in N_VALUES])
    # silence the minor-tick labeling that log scale adds by default
    ax.minorticks_off()
    ax.set_xlabel(r"cycle states $N$")
    ax.set_ylabel(r"$I_{\mathrm{struct}}(X;Y)$ (bits)")

    ax.legend(loc="lower right", frameon=False,
              handlelength=1.5, handletextpad=0.5, borderaxespad=0.3,
              labelspacing=0.3)

    # annotate the saturation collapse
    ax.text(0.04, 0.95, "all $L$ collapse onto $\\log_2 N$",
            transform=ax.transAxes, fontsize=6.5, color="grey",
            ha="left", va="top")

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
        subprocess.run(["pdfinfo", str(pdf_path)],
                       check=True, capture_output=True, text=True)
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
    series = loadB()
    fig = plotPanelA(series)
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=600)
    validate(PDF_PATH, PNG_PATH, fig)
    plt.close(fig)
    print(f"OK panel A:")
    print(f"  pdf: {PDF_PATH} ({PDF_PATH.stat().st_size / 1024:.1f} kB)")
    print(f"  png: {PNG_PATH} ({PNG_PATH.stat().st_size / 1024:.1f} kB)")


if __name__ == "__main__":
    main()
