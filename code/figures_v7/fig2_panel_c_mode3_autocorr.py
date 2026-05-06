"""fig 2 panel c: Mode 3 output autocorrelation peaks at lag = N.

The Mode 3 channel produces y_k = (phi_0 + k) mod N + noise. With phi_0 fixed
within a trajectory, the monomer indicator y_k = y_{k+N} (perfect period-N
periodicity). The autocorrelation of any monomer-indicator function of Y
therefore peaks at every multiple of N — most informatively at the first
nonzero peak, lag = N.

This panel simulates Mode 3 on the fly (the test_b CSV stores only the
peak-lag scalar, not the full lag vs autocorr curve), with epsilon=0.001 to
keep periodicity sharp. Source simulator parameters mirror the Mode 3
definition in code/test_b_mode3_capacity.py (independently re-implemented
here per the repo's "tests are self-contained" rule, applied to the panel
script as well so it does not break if test_b is moved).

output: paper/figures/v7/fig2/fig2_panel_c.{pdf,png}

reader takeaway: the autocorrelation of the Mode 3 output sequence peaks at
lag = N for each cycle-state count N -- the empirical signature of the
N-bounded saturation in panel A and B.
"""
import sys
import subprocess
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
from templating_style import apply_paper_style, OKABE_ITO  # noqa: E402

apply_paper_style()

# paths
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
OUT_DIR = PROJECT_ROOT / "paper" / "figures" / "v7" / "fig2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PDF_PATH = OUT_DIR / "fig2_panel_c.pdf"
PNG_PATH = OUT_DIR / "fig2_panel_c.png"

# panel parameters per dispatch
N_VALUES = [2, 3, 4, 6]
EPSILON = 0.001          # low-noise so periodicity is clean and visible
L_TRAJECTORY = 800       # plenty of cycles to average over
N_TRAJECTORIES = 200     # independent (phi_0 draws) trajectories to average autocorr
MAX_LAG = 12             # x-axis range 0..MAX_LAG (covers two full periods even at N=6)
SEED = 42


def simulateMode3(N, L, epsilon, n_traj, rng):
    """re-implementation of the Mode 3 channel from test_b_mode3_capacity.py.

    Y[i, k] = (phi_0[i] + k) mod N with iid per-position miscoding at rate eps.
    """
    phi0 = rng.integers(0, N, size=n_traj, dtype=np.int64)
    pos = np.arange(L, dtype=np.int64)[None, :]
    Y_correct = (phi0[:, None] + pos) % N

    if epsilon == 0.0 or N == 1:
        return Y_correct

    is_error = rng.random(size=(n_traj, L)) < epsilon
    # uniform among N-1 wrong choices (skip the correct one)
    offset = rng.integers(0, N - 1, size=(n_traj, L), dtype=np.int64)
    Y_wrong = np.where(offset >= Y_correct, offset + 1, offset)
    return np.where(is_error, Y_wrong, Y_correct)


def autocorrelation(Y, max_lag):
    """sample-wise autocorrelation of an integer-valued trajectory, averaged
    over independent trajectories.

    For an integer sequence we use the natural Pearson correlation between
    Y[k] and Y[k+lag] (centered, normalized). For Mode 3 with N>2 this is
    only one of several valid choices; cross-checking with one-hot indicator
    autocorrelation gives the same period-N peak structure, so we use the
    simpler scalar form.
    """
    n_traj, L = Y.shape
    Yf = Y.astype(np.float64)
    # subtract per-trajectory mean (each trajectory has its own phi_0)
    mu = Yf.mean(axis=1, keepdims=True)
    Yc = Yf - mu
    var = (Yc ** 2).mean(axis=1)  # per-trajectory variance

    lags = np.arange(0, max_lag + 1)
    rs = np.zeros_like(lags, dtype=np.float64)

    for j, lag in enumerate(lags):
        if lag == 0:
            rs[j] = 1.0
            continue
        prod = Yc[:, :L - lag] * Yc[:, lag:]
        cov_per = prod.mean(axis=1)
        # normalize per trajectory then average across trajectories
        valid = var > 1e-12
        if not valid.any():
            rs[j] = 0.0
            continue
        r_per = np.zeros(n_traj)
        r_per[valid] = cov_per[valid] / var[valid]
        rs[j] = r_per[valid].mean()
    return lags, rs


def plotPanelC(curves):
    fig, ax = plt.subplots(figsize=(85 / 25.4, 60 / 25.4))

    viridis = cm.get_cmap("viridis")
    colors = [viridis(t) for t in np.linspace(0.15, 0.85, len(N_VALUES))]

    # zero-line reference
    ax.axhline(0.0, color="grey", linestyle=":", linewidth=0.6, zorder=0)

    for i, N in enumerate(N_VALUES):
        lags, rs = curves[N]
        ax.plot(lags, rs, color=colors[i], marker="o", markersize=2.8,
                linewidth=1.3, label=f"$N={N}$", zorder=3 + i)
        # mark the lag=N peak
        peak_idx = N
        if peak_idx <= MAX_LAG:
            ax.plot([peak_idx], [rs[peak_idx]], marker="o", markersize=5.5,
                    markerfacecolor="none", markeredgecolor=colors[i],
                    markeredgewidth=1.0, zorder=10)

    ax.set_xlabel("lag (positions)")
    ax.set_ylabel("autocorrelation")
    ax.set_xlim(-0.5, MAX_LAG + 0.5)
    ax.set_xticks(range(0, MAX_LAG + 1, 2))

    # tighten y range based on data, leave headroom above for legend
    all_y = np.concatenate([rs for _, rs in curves.values()])
    ymin = min(-1.05, all_y.min() - 0.05)
    ymax = 1.50
    ax.set_ylim(ymin, ymax)
    ax.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])

    ax.legend(loc="upper center", frameon=False,
              handlelength=1.2, handletextpad=0.4, borderaxespad=0.3,
              labelspacing=0.2, ncol=4, columnspacing=0.8,
              bbox_to_anchor=(0.5, 1.02))

    # annotate "peak at lag = N" in a corner that doesn't collide with curves
    ax.text(0.03, 0.05, "rings: peak at lag $= N$",
            transform=ax.transAxes, fontsize=6.5, color="grey",
            ha="left", va="bottom")

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
    rng = np.random.default_rng(SEED)
    curves = {}
    for N in N_VALUES:
        Y = simulateMode3(N, L_TRAJECTORY, EPSILON, N_TRAJECTORIES, rng)
        lags, rs = autocorrelation(Y, MAX_LAG)
        curves[N] = (lags, rs)
        # sanity: lag=N should be near +1 (perfect periodicity at low eps)
        peak_at_N = rs[N] if N <= MAX_LAG else float("nan")
        print(f"  N={N}: r(lag=N)={peak_at_N:.4f}")

    fig = plotPanelC(curves)
    fig.savefig(PDF_PATH)
    fig.savefig(PNG_PATH, dpi=600)
    validate(PDF_PATH, PNG_PATH, fig)
    plt.close(fig)
    print(f"OK panel C:")
    print(f"  pdf: {PDF_PATH} ({PDF_PATH.stat().st_size / 1024:.1f} kB)")
    print(f"  png: {PNG_PATH} ({PNG_PATH.stat().st_size / 1024:.1f} kB)")


if __name__ == "__main__":
    main()
