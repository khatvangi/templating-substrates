"""
test A.1 — Mode 1 (Watson-Crick) classification via mutual information scaling

Templating Substrates Framework, Test A.1
-----------------------------------------
Purpose: Verify that for a Mode 1 templating system (DNA-replication-like
Markov chain), the transferred structural information I_struct(X; Y) scales
linearly with template length L, with slope ≈ 2 - H(Y|X) bits per position.

Mode 1 channel:
    X is uniform over {A,C,G,T} = {0,1,2,3} at every position
    Watson-Crick pairing: A↔T, C↔G    →    wc = [3, 2, 1, 0]
    P(Y_i = wc[X_i] | X_i)        = 1 - epsilon       (correct incorporation)
    P(Y_i = each wrong | X_i)     = epsilon / 3       (uniform misincorporation)

Per-position information (closed form, since X uniform makes Y uniform):
    H(Y|X)            = -(1-eps) log2(1-eps) - eps log2(eps/3)
    I_per_position    = log2(4) - H(Y|X) = 2 - H(Y|X)
    I_total           = L * I_per_position

PASS criterion: for all (epsilon, L) with L >= 25,
    |I_empirical - I_theoretical| / I_theoretical < 0.10
"""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless backend; we save to PNG
import matplotlib.pyplot as plt
import numpy as np


# ----------------------------------------------------------------------------
# directory layout — script lives in code/, results/figures siblings
# ----------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = PROJECT_ROOT / "figures"

# Watson-Crick complement table:
#   A=0 ↔ T=3
#   C=1 ↔ G=2
WC_PAIR = np.array([3, 2, 1, 0], dtype=np.int8)

EPS_CLIP = 1e-12  # for log2 stability


# ----------------------------------------------------------------------------
# simulator
# ----------------------------------------------------------------------------
def simulate_mode1(L, epsilon, n_samples, rng=None):
    """
    generate n_samples independent (X, Y) pairs from the Mode 1 channel.

    returns
    -------
    X, Y : ndarray of shape (n_samples, L), dtype int8, values in {0,1,2,3}
    """
    if rng is None:
        rng = np.random.default_rng()

    # template uniformly random
    X = rng.integers(0, 4, size=(n_samples, L), dtype=np.int8)

    # the "intended" product is the WC complement at every position
    Y_correct = WC_PAIR[X]

    # for each position decide correct vs error
    is_error = rng.random(size=(n_samples, L)) < epsilon

    # for error positions: pick uniformly among the 3 wrong nucleotides
    # trick: pick offset in {0,1,2}, then if offset >= correct value, shift up
    # by 1 so we skip exactly the correct value, leaving the 3 wrong choices
    offset = rng.integers(0, 3, size=(n_samples, L), dtype=np.int8)
    Y_wrong = np.where(offset >= Y_correct, offset + 1, offset).astype(np.int8)

    Y = np.where(is_error, Y_wrong, Y_correct).astype(np.int8)
    return X, Y


# ----------------------------------------------------------------------------
# empirical MI estimator (plug-in, position by position)
# ----------------------------------------------------------------------------
def estimate_mutual_information(X, Y):
    """
    plug-in estimator of I(X;Y) decomposed across positions.

    for each column i, build the empirical 4x4 joint over (x_i, y_i) and
    compute I_i = sum P(x,y) log2 [ P(x,y) / (P(x) P(y)) ].
    return the position sum and per-position mean.
    """
    n_samples, L = X.shape
    # encode (x,y) as a single int 0..15 for fast bincount
    pair_codes = X.astype(np.int64) * 4 + Y.astype(np.int64)  # (n_samples, L)

    I_total = 0.0
    I_per_pos = np.zeros(L)
    for i in range(L):
        counts = np.bincount(pair_codes[:, i], minlength=16).reshape(4, 4)
        joint = counts / n_samples
        p_x = joint.sum(axis=1, keepdims=True)
        p_y = joint.sum(axis=0, keepdims=True)
        denom = p_x * p_y
        mask = (joint > 0) & (denom > 0)
        # plug-in MI contribution from this position
        I_i = np.sum(joint[mask] * np.log2(joint[mask] / denom[mask]))
        I_per_pos[i] = I_i
        I_total += I_i
    return I_total, I_total / L, I_per_pos


# ----------------------------------------------------------------------------
# closed-form theoretical reference
# ----------------------------------------------------------------------------
def theoretical_per_position_info(epsilon):
    """
    closed form: I_per_position = 2 - H(Y|X)
    with H(Y|X) = -(1-eps) log2(1-eps) - eps log2(eps/3).
    """
    if epsilon <= 0.0:
        H_y_given_x = 0.0
    else:
        one_minus = max(1.0 - epsilon, EPS_CLIP)
        eps_third = max(epsilon / 3.0, EPS_CLIP)
        H_y_given_x = -(1 - epsilon) * np.log2(one_minus) - epsilon * np.log2(eps_third)
    return 2.0 - H_y_given_x


# ----------------------------------------------------------------------------
# experiment driver
# ----------------------------------------------------------------------------
def run_experiment():
    np.random.seed(42)
    rng = np.random.default_rng(42)

    epsilons = [0.001, 0.01, 0.05, 0.10, 0.25]
    Ls = [10, 25, 50, 100, 200, 500]
    n_samples = 5000

    rows = []  # for CSV
    print(f"{'epsilon':>8s}  {'L':>4s}  "
          f"{'I_emp':>10s}  {'I_theo':>10s}  "
          f"{'ipp_emp':>10s}  {'ipp_theo':>10s}  "
          f"{'%err':>8s}")
    print("-" * 70)

    for eps in epsilons:
        ipp_theo = theoretical_per_position_info(eps)
        for L in Ls:
            X, Y = simulate_mode1(L, eps, n_samples, rng=rng)
            I_emp_total, I_emp_pp, _ = estimate_mutual_information(X, Y)
            I_theo_total = L * ipp_theo

            # relative error guarded against tiny denom (won't happen for our eps grid)
            denom = max(abs(I_theo_total), EPS_CLIP)
            pct_err = 100.0 * (I_emp_total - I_theo_total) / denom

            rows.append({
                "epsilon": eps,
                "L": L,
                "I_empirical": I_emp_total,
                "I_theoretical": I_theo_total,
                "I_per_position_empirical": I_emp_pp,
                "I_per_position_theoretical": ipp_theo,
            })
            print(f"{eps:8.4f}  {L:4d}  "
                  f"{I_emp_total:10.4f}  {I_theo_total:10.4f}  "
                  f"{I_emp_pp:10.4f}  {ipp_theo:10.4f}  "
                  f"{pct_err:+7.2f}%")

    return rows, epsilons, Ls


# ----------------------------------------------------------------------------
# I/O — CSV
# ----------------------------------------------------------------------------
def save_csv(rows, path):
    fieldnames = [
        "epsilon", "L",
        "I_empirical", "I_theoretical",
        "I_per_position_empirical", "I_per_position_theoretical",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


# ----------------------------------------------------------------------------
# plotting (matplotlib only, per task constraint)
# ----------------------------------------------------------------------------
def plot_scaling(rows, epsilons, Ls, path):
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    cmap = plt.get_cmap("viridis")

    for k, eps in enumerate(epsilons):
        xs, y_emp, y_theo = [], [], []
        for r in rows:
            if r["epsilon"] == eps:
                xs.append(r["L"])
                y_emp.append(r["I_empirical"])
                y_theo.append(r["I_theoretical"])
        color = cmap(k / max(1, len(epsilons) - 1))
        # theoretical: solid line
        L_dense = np.linspace(min(Ls), max(Ls), 100)
        ipp = theoretical_per_position_info(eps)
        ax.plot(L_dense, ipp * L_dense, color=color, linestyle="-", alpha=0.7,
                label=f"theory  eps={eps}")
        # empirical: scatter
        ax.scatter(xs, y_emp, color=color, marker="o", s=45,
                   edgecolor="black", linewidth=0.6, zorder=3,
                   label=f"empirical eps={eps}")

    ax.set_xlabel("template length L")
    ax.set_ylabel(r"$I_{\mathrm{struct}}(X;Y)$  [bits]")
    ax.set_title("Test A.1 — Mode 1 Watson-Crick: I(X;Y) scales linearly with L")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7, ncol=2, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_per_position(rows, epsilons, path):
    fig, ax = plt.subplots(figsize=(7.0, 5.0))

    # empirical points: average per-position info across L for each eps,
    # with min/max as error bar to show consistency
    for eps in epsilons:
        ipps = [r["I_per_position_empirical"]
                for r in rows if r["epsilon"] == eps]
        ipps = np.array(ipps)
        ax.errorbar(eps, ipps.mean(),
                    yerr=[[ipps.mean() - ipps.min()],
                          [ipps.max() - ipps.mean()]],
                    fmt="o", color="C0", markersize=8, capsize=4,
                    label="empirical (mean ± range over L)" if eps == epsilons[0] else None)

    # theoretical curve
    eps_dense = np.linspace(1e-4, 0.74, 400)
    ipp_curve = np.array([theoretical_per_position_info(e) for e in eps_dense])
    ax.plot(eps_dense, ipp_curve, "r-", lw=1.7, label="theory: 2 - H(Y|X)")

    ax.set_xlabel(r"misincorporation rate $\varepsilon$")
    ax.set_ylabel(r"$I$ per position  [bits]")
    ax.set_title("Test A.1 — per-position information vs misincorporation rate")
    ax.set_xscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ----------------------------------------------------------------------------
# pass / fail logic
# ----------------------------------------------------------------------------
def evaluate_pass(rows, tol=0.10, min_L=25):
    """
    PASS if for every (eps, L) with L >= min_L:
        |I_empirical - I_theoretical| / I_theoretical < tol
    """
    failures = []
    for r in rows:
        if r["L"] < min_L:
            continue
        denom = max(abs(r["I_theoretical"]), EPS_CLIP)
        rel = abs(r["I_empirical"] - r["I_theoretical"]) / denom
        if rel >= tol:
            failures.append((r["epsilon"], r["L"], rel))
    return failures


# ----------------------------------------------------------------------------
# README
# ----------------------------------------------------------------------------
README_TEXT = """\
# Test A.1 — Mode 1 (Watson-Crick) Classification via Mutual Information Scaling

## What this test does

Validates the framework's diagnostic apparatus on the textbook Mode 1 channel
(DNA-replication-like Watson-Crick base pairing). For a template X of length L
drawn uniformly over {A,C,G,T} and a product Y obtained by independently
incorporating the WC complement with probability `1 - epsilon` and a uniformly
chosen wrong nucleotide with probability `epsilon` (split equally over the 3
wrong choices), the transferred structural information should be:

    I_struct(X; Y)  =  L * (2 - H(Y|X))   bits
    H(Y|X)          =  -(1-eps) log2(1-eps) - eps log2(eps/3)

That is: I scales **linearly in L**, with slope set entirely by the
per-position channel.

## How the test works

1. `simulate_mode1(L, epsilon, n_samples)` draws n_samples independent
   (X, Y) pairs from the Mode 1 channel, vectorized over n_samples and L.
2. `estimate_mutual_information(X, Y)` plug-in-estimates I(x_i; y_i) per
   position from the empirical 4x4 joint and sums across positions.
3. `theoretical_per_position_info(epsilon)` returns the closed-form
   I_per_position = 2 - H(Y|X).
4. We sweep epsilon ∈ {0.001, 0.01, 0.05, 0.1, 0.25} and
   L ∈ {10, 25, 50, 100, 200, 500} with n_samples = 5000.
5. We compare empirical vs theoretical at every (epsilon, L), generate
   plots, and decide PASS/FAIL.

## How to interpret the results

- `figures/test_a1_scaling.png` — empirical points should fall on the
  theoretical lines on a log-log plot, with slope ~1 in log-log (linear
  scaling in L). Different epsilon values give parallel lines offset by
  the per-position information.
- `figures/test_a1_per_position.png` — empirical mean (with min/max range
  across L) should sit on the theoretical 2 - H(Y|X) curve.
- `results/test_a1_results.csv` — tabular data, one row per (epsilon, L).

## PASS / FAIL criterion

PASS if for every (epsilon, L) with L >= 25,
    |I_empirical - I_theoretical| / I_theoretical < 10%.

The L >= 25 cutoff filters out the small-L regime where finite-sample
bias of the plug-in MI estimator can dominate the true signal.

## What a PASS means

The diagnostic apparatus correctly captures length-scaling templating
information for the canonical Mode 1 channel. We can now apply it to
harder cases: Test A.2 (bulk-matched control) and Test B (Mode 3
capacity prediction).

## What a FAIL means

Either the simulator is wrong, the estimator is biased beyond tolerance,
or the framework's specification of Mode 1 is inconsistent with the
test. The script's stdout summary localizes which (epsilon, L) cells
failed, which usually narrows the diagnosis.

## Reproducibility

`np.random.seed(42)` and an explicitly seeded `default_rng(42)` are set
at the start of `run_experiment()`. Re-running the script reproduces
the table exactly.
"""


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main():
    rows, epsilons, Ls = run_experiment()

    csv_path = RESULTS_DIR / "test_a1_results.csv"
    save_csv(rows, csv_path)

    scaling_png = FIGURES_DIR / "test_a1_scaling.png"
    pp_png = FIGURES_DIR / "test_a1_per_position.png"
    plot_scaling(rows, epsilons, Ls, scaling_png)
    plot_per_position(rows, epsilons, pp_png)

    readme_path = RESULTS_DIR / "test_a1_README.md"
    readme_path.write_text(README_TEXT)

    failures = evaluate_pass(rows)
    print("-" * 70)
    print(f"results saved to: {csv_path}")
    print(f"figures saved to: {scaling_png}, {pp_png}")
    print(f"readme  saved to: {readme_path}")
    print("-" * 70)
    if failures:
        print(f"FAIL  ({len(failures)} cells outside 10% tolerance, L >= 25):")
        for eps, L, rel in failures:
            print(f"    epsilon={eps}, L={L}  →  rel_err={rel:.3%}")
    else:
        print("PASS  — all (epsilon, L) with L >= 25 are within 10% of theory")


if __name__ == "__main__":
    main()
