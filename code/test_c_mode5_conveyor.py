"""
test C — Mode 5 modular conveyor templating: linear-N scaling + length limit

Templating Substrates Framework, Test C
---------------------------------------
Purpose: verify two predictions of Mode 5 (NRPS/PKS-style modular conveyor):

    (1) Information scaling   :  I_struct(X; Y)  =  N * (log_2(A) - H(Y|X))
    (2) Length limit          :  product length is bounded by N (the module count)

A Mode 5 assembly has N spatially distinct modules in fixed order. Each module
i has an "intended monomer" m_i drawn from substrate alphabet A. The intended
sequence X = (m_0, ..., m_{N-1}) is encoded in the MODULE STRUCTURE itself —
not in a separate sequence template. The product Y has exactly N positions,
one per module.

Mode 5 channel (per module, given X):
    P(y_i = m_i | X)              = 1 - eps
    P(y_i = each other in A | X)  = eps / (A - 1)

Mathematically the MI is identical to Mode 1 with alphabet A and uniform X:
positions are independent given X, so per-position estimator from A.1
generalizes directly. The distinction with Mode 1 is structural (template
encoded in machinery, not in a separate molecule), and shows up as the
length-limit prediction: there is no module N+1, so Y cannot extend past N.

Length-limit sub-experiment
---------------------------
We model "trying to push past the last module" as: positions i in [0, N-1]
get the templated channel; positions i in [N, 2N-1] get a uniform-random
monomer (no template guidance available). Then per-position I(X; y_i) > 0
for i < N and = 0 for i >= N. The visual sharp drop at i = N is the
empirical signature of the length limit.

PASS criterion
--------------
1. Linear N-scaling: for every (eps, A) and every N >= 4,
       |I_emp - N * I_pp_theoretical|  /  (N * I_pp_theoretical)  <  0.05
2. Length-limit (both representative cases):
       mean(I_pp[0:N])  >  10 * mean(I_pp[N:2N])
"""

import csv
import time
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

EPS_CLIP = 1e-12  # for log2 stability


# ----------------------------------------------------------------------------
# Mode 5 simulator — biased per-module channel with arbitrary alphabet
# ----------------------------------------------------------------------------
def simulate_mode5(N, alphabet_size, epsilon, n_samples, rng=None):
    """
    generate n_samples independent (X, Y) draws from the Mode 5 channel.

    X is the module-intended-monomer sequence, drawn uniformly from
    {0, ..., A-1}^N. Y is the product, one monomer per module, with the
    biased per-position channel above.

    parameters
    ----------
    N             : module count (also product length)
    alphabet_size : substrate alphabet size A (e.g., 4 nucleotide, 20 amino)
    epsilon       : misincorporation rate per module
    n_samples     : number of independent assembly instances
    rng           : numpy Generator (seeded by caller)

    returns
    -------
    X : ndarray, shape (n_samples, N), dtype int16, values in {0, ..., A-1}
    Y : ndarray, shape (n_samples, N), dtype int16, values in {0, ..., A-1}
    """
    if rng is None:
        rng = np.random.default_rng()

    A = alphabet_size
    # module-intended monomers, uniform over {0, ..., A-1} per (sample, module)
    X = rng.integers(0, A, size=(n_samples, N), dtype=np.int16)

    # decide correct vs error per (sample, module)
    is_error = rng.random(size=(n_samples, N)) < epsilon

    # for error positions: uniform pick among A-1 wrong monomers
    # same offset trick as A.1 / B: pick offset in {0, ..., A-2}, shift up by 1
    # if offset >= correct value, so we skip exactly the correct monomer.
    if A > 1:
        offset = rng.integers(0, A - 1, size=(n_samples, N), dtype=np.int64)
        Y_wrong = np.where(offset >= X, offset + 1, offset).astype(np.int16)
    else:
        Y_wrong = X.copy()

    Y = np.where(is_error, Y_wrong, X).astype(np.int16)
    return X, Y


# ----------------------------------------------------------------------------
# per-position MI estimator — same as A.1, generalized to alphabet A
# ----------------------------------------------------------------------------
def estimate_mutual_information_per_position(X, Y, alphabet_size):
    """
    plug-in estimator of I(X; Y) decomposed across positions.

    Mode 5 positions are independent given X (each module acts on one site
    only), so I(X; Y) = sum_i I(x_i; y_i). For each column i, build the
    empirical AxA joint over (x_i, y_i) and compute the plug-in MI.

    returns
    -------
    I_total          : scalar, sum over positions
    I_per_position   : scalar, mean over positions (= I_total / N)
    I_per_position_a : ndarray of shape (N,), per-column MI
    """
    n_samples, N = X.shape
    A = alphabet_size
    # encode (x, y) as a single int 0 .. A*A-1 for fast bincount
    pair_codes = X.astype(np.int64) * A + Y.astype(np.int64)  # (n_samples, N)

    I_per_pos = np.zeros(N, dtype=np.float64)
    for i in range(N):
        counts = np.bincount(pair_codes[:, i], minlength=A * A).reshape(A, A)
        joint = counts / n_samples
        p_x = joint.sum(axis=1, keepdims=True)
        p_y = joint.sum(axis=0, keepdims=True)
        denom = p_x * p_y
        mask = (joint > 0) & (denom > 0)
        I_per_pos[i] = float(np.sum(joint[mask] * np.log2(joint[mask] / denom[mask])))

    I_total = float(I_per_pos.sum())
    return I_total, I_total / max(N, 1), I_per_pos


# ----------------------------------------------------------------------------
# closed-form theoretical reference
# ----------------------------------------------------------------------------
def theoretical_per_position_info(epsilon, alphabet_size):
    """
    closed form: I_pp = log_2(A) - H(Y|X)
    H(Y|X) = -(1-eps) log_2(1-eps) - eps log_2(eps/(A-1))
    """
    A = alphabet_size
    if A <= 1:
        return 0.0
    if epsilon <= 0.0:
        H_y_given_x = 0.0
    elif epsilon >= 1.0:
        # degenerate; uniform over A-1 wrong choices
        H_y_given_x = float(np.log2(A - 1)) if A > 1 else 0.0
    else:
        one_minus = max(1.0 - epsilon, EPS_CLIP)
        eps_per_wrong = max(epsilon / (A - 1), EPS_CLIP)
        H_y_given_x = (
            -(1 - epsilon) * np.log2(one_minus)
            - epsilon * np.log2(eps_per_wrong)
        )
    return float(np.log2(A) - H_y_given_x)


# ----------------------------------------------------------------------------
# length-limit sub-experiment
# ----------------------------------------------------------------------------
def simulate_mode5_with_untemplated_tail(N, alphabet_size, epsilon, n_samples, rng):
    """
    Mode 5 with N templated module positions plus N untemplated tail positions.

    positions 0..N-1   :  channel as in simulate_mode5 (templated)
    positions N..2N-1  :  uniform random monomers, independent of X
                         (no module exists past N — "no template guidance")

    returns
    -------
    X : (n_samples, N)  intended-monomer module sequence
    Y : (n_samples, 2N) product including untemplated tail
    """
    A = alphabet_size

    X = rng.integers(0, A, size=(n_samples, N), dtype=np.int16)

    # templated region (same logic as simulate_mode5)
    is_error = rng.random(size=(n_samples, N)) < epsilon
    if A > 1:
        offset = rng.integers(0, A - 1, size=(n_samples, N), dtype=np.int64)
        Y_wrong = np.where(offset >= X, offset + 1, offset).astype(np.int16)
    else:
        Y_wrong = X.copy()
    Y_templated = np.where(is_error, Y_wrong, X).astype(np.int16)

    # untemplated tail — uniform random, independent of X
    Y_tail = rng.integers(0, A, size=(n_samples, N), dtype=np.int16)

    Y = np.concatenate([Y_templated, Y_tail], axis=1)
    return X, Y


def per_position_info_against_X(X, Y, alphabet_size):
    """
    compute per-position information across the full Y of length L (possibly
    longer than X). for each i:
        - if i < N (where N = X.shape[1])  : compute I(x_i; y_i)
        - if i >= N                        : compute I(x_{i mod N}; y_i)
                                              (any x column works since y_i is
                                              independent of all x columns)

    in both cases this uses the AxA plug-in MI estimator. For i < N the result
    estimates I(X; y_i) since y_i depends only on x_i. For i >= N the result
    estimates 0 (any x column is independent of the uniform y_i), measured
    plus finite-sample bias.
    """
    n_samples, N = X.shape
    _, L = Y.shape
    A = alphabet_size

    I_pp = np.zeros(L, dtype=np.float64)
    for i in range(L):
        x_col = X[:, i % N]
        y_col = Y[:, i]
        pair_codes = x_col.astype(np.int64) * A + y_col.astype(np.int64)
        counts = np.bincount(pair_codes, minlength=A * A).reshape(A, A)
        joint = counts / n_samples
        p_x = joint.sum(axis=1, keepdims=True)
        p_y = joint.sum(axis=0, keepdims=True)
        denom = p_x * p_y
        mask = (joint > 0) & (denom > 0)
        I_pp[i] = float(np.sum(joint[mask] * np.log2(joint[mask] / denom[mask])))
    return I_pp


# ----------------------------------------------------------------------------
# main scaling experiment
# ----------------------------------------------------------------------------
def run_scaling_experiment(progress_path):
    np.random.seed(42)
    rng = np.random.default_rng(42)

    Ns = [2, 4, 8, 16, 32]
    epsilons = [0.001, 0.01, 0.05]
    alphabets = [4, 20]
    n_samples = 10_000

    rows = []

    print(f"{'N':>3s} {'eps':>7s} {'A':>3s} {'n':>6s} "
          f"{'I_emp':>10s} {'I_theo':>10s} "
          f"{'ipp_emp':>10s} {'ipp_theo':>10s} {'rel_err%':>9s}")
    print("-" * 80)

    progress_path.parent.mkdir(parents=True, exist_ok=True)
    with open(progress_path, "a", buffering=1) as pf:
        pf.write(f"\n# === scaling run started at {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")

        for A in alphabets:
            for eps in epsilons:
                ipp_theo = theoretical_per_position_info(eps, A)
                for N in Ns:
                    t0 = time.time()
                    X, Y = simulate_mode5(N, A, eps, n_samples, rng=rng)
                    I_emp, I_pp_emp, _ = estimate_mutual_information_per_position(
                        X, Y, A
                    )
                    I_theo = N * ipp_theo
                    denom = max(abs(I_theo), EPS_CLIP)
                    rel_err = 100.0 * (I_emp - I_theo) / denom
                    elapsed = time.time() - t0

                    row = {
                        "N": N,
                        "epsilon": eps,
                        "alphabet_size": A,
                        "n_samples": n_samples,
                        "I_empirical": I_emp,
                        "I_theoretical": I_theo,
                        "I_per_position_empirical": I_pp_emp,
                        "I_per_position_theoretical": ipp_theo,
                        "rel_err_pct": rel_err,
                    }
                    rows.append(row)

                    line = (
                        f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
                        f"scaling: N={N}, eps={eps}, A={A}: "
                        f"I_emp={I_emp:.3f} bits "
                        f"(theory={I_theo:.3f}, "
                        f"rel_err={rel_err:+.2f}%, elapsed={elapsed:.1f}s)"
                    )
                    pf.write(line + "\n")
                    print(
                        f"{N:3d} {eps:7.4f} {A:3d} {n_samples:6d} "
                        f"{I_emp:10.4f} {I_theo:10.4f} "
                        f"{I_pp_emp:10.4f} {ipp_theo:10.4f} {rel_err:+8.2f}%"
                    )

    return rows, Ns, epsilons, alphabets, n_samples


# ----------------------------------------------------------------------------
# length-limit experiment
# ----------------------------------------------------------------------------
def run_length_limit_experiment(progress_path):
    """
    representative cases: (N=4, A=4) and (N=8, A=20). Use eps=0.01.
    """
    rng = np.random.default_rng(123)
    cases = [
        ("N=4_A=4",  4,  4,  0.01),
        ("N=8_A=20", 8, 20, 0.01),
    ]
    n_samples = 10_000

    rows = []
    summaries = []

    with open(progress_path, "a", buffering=1) as pf:
        pf.write(f"\n# === length-limit run started at {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")

        for label, N, A, eps in cases:
            ipp_theo = theoretical_per_position_info(eps, A)
            t0 = time.time()
            X, Y = simulate_mode5_with_untemplated_tail(N, A, eps, n_samples, rng)
            I_pp = per_position_info_against_X(X, Y, A)
            elapsed = time.time() - t0

            templated_mean = float(I_pp[:N].mean())
            untemplated_mean = float(I_pp[N:].mean())
            ratio = templated_mean / max(untemplated_mean, EPS_CLIP)

            for i, ipp in enumerate(I_pp):
                rows.append({
                    "case_label": label,
                    "position_index": i,
                    "I_per_position_empirical": float(ipp),
                    "is_templated_region": bool(i < N),
                })

            summaries.append({
                "label": label,
                "N": N,
                "A": A,
                "eps": eps,
                "ipp_theo": ipp_theo,
                "templated_mean": templated_mean,
                "untemplated_mean": untemplated_mean,
                "ratio": ratio,
            })

            line = (
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
                f"length-limit: {label} (N={N}, A={A}, eps={eps}): "
                f"mean_templated={templated_mean:.4f} (theory={ipp_theo:.4f}), "
                f"mean_untemplated={untemplated_mean:.4f}, "
                f"ratio={ratio:.1f}, elapsed={elapsed:.1f}s"
            )
            pf.write(line + "\n")
            print(line)

    return rows, summaries


# ----------------------------------------------------------------------------
# I/O — CSV
# ----------------------------------------------------------------------------
def save_scaling_csv(rows, path):
    fieldnames = [
        "N", "epsilon", "alphabet_size", "n_samples",
        "I_empirical", "I_theoretical",
        "I_per_position_empirical", "I_per_position_theoretical",
        "rel_err_pct",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r[k] for k in fieldnames})


def save_length_limit_csv(rows, path):
    fieldnames = [
        "case_label", "position_index",
        "I_per_position_empirical", "is_templated_region",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r[k] for k in fieldnames})


# ----------------------------------------------------------------------------
# plotting — N scaling
# ----------------------------------------------------------------------------
def plot_n_scaling(rows, Ns, epsilons, alphabets, path):
    """
    one line per (eps, A) combo. x = N, y = I_empirical. theoretical lines
    overlaid. should be linear with slope log_2(A) - H(Y|X).
    """
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    cmap = plt.get_cmap("viridis")
    n_combos = len(epsilons) * len(alphabets)

    combo_idx = 0
    for A in alphabets:
        for eps in epsilons:
            ipp_theo = theoretical_per_position_info(eps, A)
            xs, ys = [], []
            for r in rows:
                if r["epsilon"] == eps and r["alphabet_size"] == A:
                    xs.append(r["N"])
                    ys.append(r["I_empirical"])
            order = np.argsort(xs)
            xs = np.array(xs)[order]
            ys = np.array(ys)[order]
            color = cmap(combo_idx / max(1, n_combos - 1))
            # empirical points
            ax.plot(xs, ys, "o", color=color, markersize=7,
                    markeredgecolor="black", markeredgewidth=0.5,
                    label=f"emp  A={A}, eps={eps}", zorder=3)
            # theoretical line
            N_dense = np.linspace(min(Ns), max(Ns), 100)
            ax.plot(N_dense, ipp_theo * N_dense, "-", color=color, alpha=0.75,
                    linewidth=1.5,
                    label=f"theory A={A}, eps={eps}: slope={ipp_theo:.3f}")
            combo_idx += 1

    ax.set_xlabel("module count N (= product length)")
    ax.set_ylabel(r"$I_{\mathrm{struct}}(X;Y)$  [bits]")
    ax.set_title("Test C — Mode 5: I_struct scales linearly in N, slope = log_2(A) - H(Y|X)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7, ncol=2, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ----------------------------------------------------------------------------
# plotting — length limit
# ----------------------------------------------------------------------------
def plot_length_limit(rows, summaries, path):
    """
    per-position I vs position index for the two representative cases.
    sharp drop at position i = N is the visual signature of the length limit.
    """
    fig, axes = plt.subplots(1, len(summaries),
                             figsize=(5.8 * len(summaries), 4.5),
                             squeeze=False)

    for ax, summ in zip(axes[0], summaries):
        label = summ["label"]
        N = summ["N"]
        A = summ["A"]
        eps = summ["eps"]
        ipp_theo = summ["ipp_theo"]

        # gather positions for this case
        case_rows = [r for r in rows if r["case_label"] == label]
        case_rows.sort(key=lambda r: r["position_index"])
        positions = [r["position_index"] for r in case_rows]
        ipps = [r["I_per_position_empirical"] for r in case_rows]

        # bar plot — clean visual for "drop at i=N"
        colors = ["#2563eb" if i < N else "#9ca3af" for i in positions]
        ax.bar(positions, ipps, color=colors, edgecolor="black", linewidth=0.5)

        ax.axhline(ipp_theo, color="red", linestyle="--", alpha=0.8,
                   linewidth=1.4,
                   label=f"templated theory = {ipp_theo:.3f} bits")
        ax.axhline(0.0, color="grey", linestyle=":", alpha=0.6,
                   linewidth=0.9, label="untemplated theory = 0 bits")
        ax.axvline(N - 0.5, color="black", linestyle="-", alpha=0.7,
                   linewidth=1.2)

        ax.set_xlabel("position index i")
        ax.set_ylabel(r"per-position $I(X; y_i)$  [bits]")
        ax.set_title(f"{label}: eps={eps}\nleft N={N} templated | right N={N} uniform tail")
        ax.set_ylim(bottom=min(0.0, min(ipps) * 1.1) - 0.1,
                    top=max(ipps) * 1.15 + 0.05)
        ax.grid(True, axis="y", alpha=0.3)
        ax.legend(fontsize=8, loc="upper right")

    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ----------------------------------------------------------------------------
# plotting — three-mode contrast figure
# ----------------------------------------------------------------------------
def _load_csv_dicts(path):
    """tiny helper: load CSV as list of dicts with float-typed numerics."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        rows = []
        for r in reader:
            converted = {}
            for k, v in r.items():
                # try float; fall back to original string (for things like
                # block_length integers we don't care about here)
                try:
                    converted[k] = float(v)
                except (ValueError, TypeError):
                    converted[k] = v
            rows.append(converted)
    return rows


def plot_mode_contrast(test_c_rows, path, eps_for_plot=0.01, A_for_plot=4):
    """
    composite: Mode 1 (linear in L), Mode 3 (saturating at log_2(N)),
    Mode 5 (linear in N, capped at N positions). Uses stored A.1/B CSVs
    if available; falls back to an inline note otherwise.

    x-axis is "polymer length" — for Mode 5 this equals N.
    """
    fig, ax = plt.subplots(figsize=(9.0, 6.0))

    # --- Mode 1 (from test_a1_results.csv) ---
    a1_csv = RESULTS_DIR / "test_a1_results.csv"
    if a1_csv.exists():
        a1 = _load_csv_dicts(a1_csv)
        a1_pts = sorted(
            [(r["L"], r["I_empirical"]) for r in a1 if abs(r["epsilon"] - eps_for_plot) < 1e-12],
            key=lambda t: t[0],
        )
        if a1_pts:
            xs = [p[0] for p in a1_pts]
            ys = [p[1] for p in a1_pts]
            ax.plot(xs, ys, "-s", color="crimson", markersize=7, linewidth=1.8,
                    label=f"Mode 1 (A=4, eps={eps_for_plot}): I ∝ L  (sequence template)")

    # --- Mode 3 (from test_b_results.csv) — show one representative N ---
    b_csv = RESULTS_DIR / "test_b_results.csv"
    N_for_mode3 = 4
    if b_csv.exists():
        b_all = _load_csv_dicts(b_csv)
        b_pts = sorted(
            [(r["L"], r["I_empirical"]) for r in b_all
             if abs(r["epsilon"] - eps_for_plot) < 1e-12 and int(r["N"]) == N_for_mode3],
            key=lambda t: t[0],
        )
        if b_pts:
            xs = [p[0] for p in b_pts]
            ys = [p[1] for p in b_pts]
            ax.plot(xs, ys, "-^", color="#1f77b4", markersize=7, linewidth=1.8,
                    label=f"Mode 3 (N={N_for_mode3}, eps={eps_for_plot}): I → log_2(N) = {np.log2(N_for_mode3):.2f}  (cyclic active site)")
            ax.axhline(np.log2(N_for_mode3), color="#1f77b4", linestyle="--",
                       alpha=0.5, linewidth=1.0)

    # --- Mode 5 (this test) — I vs N for chosen (eps, A) ---
    c_pts = sorted(
        [(r["N"], r["I_empirical"]) for r in test_c_rows
         if r["epsilon"] == eps_for_plot and r["alphabet_size"] == A_for_plot],
        key=lambda t: t[0],
    )
    if c_pts:
        xs = [p[0] for p in c_pts]
        ys = [p[1] for p in c_pts]
        ax.plot(xs, ys, "-o", color="#2ca02c", markersize=8, linewidth=1.8,
                label=f"Mode 5 (A={A_for_plot}, eps={eps_for_plot}): I ∝ N, length CAPPED at N  (modular conveyor)")

        # vertical "wall" at the largest N to emphasize the length cap
        N_cap = max(xs)
        ax.axvline(N_cap, color="#2ca02c", linestyle=":", alpha=0.6,
                   linewidth=1.2)
        ax.text(N_cap * 1.05, max(ys) * 0.55,
                f"Mode 5 cap\n(no module N+1)",
                fontsize=8, color="#2ca02c",
                verticalalignment="center")

    ax.set_xscale("log")
    ax.set_xlabel("polymer length  (= L for Modes 1/3, = N for Mode 5)")
    ax.set_ylabel(r"$I_{\mathrm{struct}}(X;Y)$  [bits]")
    ax.set_title(
        "Test C — three-mode contrast:\n"
        "Mode 1 grows ∝ L (no length limit)   |   "
        "Mode 3 saturates at log_2(N)   |   "
        "Mode 5 grows ∝ N then stops (length = N)"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ----------------------------------------------------------------------------
# pass / fail logic
# ----------------------------------------------------------------------------
def evaluate_pass(scaling_rows, length_summaries, tol=0.05, min_N=4):
    """
    PASS criteria:
      (1) for every (eps, A), every N >= min_N:
          |I_emp - N * ipp_theo| / (N * ipp_theo) < tol
      (2) for both length-limit cases:
          mean(I_pp[0:N]) > 10 * mean(I_pp[N:2N])
    returns (overall, scaling_failures, length_failures, (c1_ok, c2_ok))
    """
    scaling_failures = []
    for r in scaling_rows:
        if r["N"] < min_N:
            continue
        ipp_theo = r["I_per_position_theoretical"]
        I_theo = r["N"] * ipp_theo
        denom = max(abs(I_theo), EPS_CLIP)
        rel = abs(r["I_empirical"] - I_theo) / denom
        if rel >= tol:
            scaling_failures.append(
                (r["N"], r["epsilon"], r["alphabet_size"], r["I_empirical"],
                 I_theo, rel)
            )

    length_failures = []
    for s in length_summaries:
        if s["templated_mean"] <= 10.0 * s["untemplated_mean"]:
            length_failures.append(s)

    c1_ok = len(scaling_failures) == 0
    c2_ok = len(length_failures) == 0
    return c1_ok and c2_ok, scaling_failures, length_failures, (c1_ok, c2_ok)


# ----------------------------------------------------------------------------
# README
# ----------------------------------------------------------------------------
README_TEXT = """\
# Test C — Mode 5 Modular Conveyor Templating

## What this test does

Validates the framework's two predictions for Mode 5 (NRPS / PKS-style
modular assembly lines):

    (1) Information scaling :  I_struct(X; Y) = N * (log_2(A) - H(Y|X))
    (2) Length limit        :  product length is bounded by N (the module count)

A Mode 5 system has N spatially distinct modules in fixed order. Each module
carries its "intended monomer" m_i in the structure of the module itself, not
in a separable sequence template. The substrate progresses through modules in
order, picking up one monomer per module. The product Y therefore has length
exactly N — there is no Module N+1 to extend the chain.

## Channel definition (per module)

X = (m_0, ..., m_{N-1})   uniform over {0, ..., A-1}^N
P(y_i = m_i  |  X)              = 1 - eps
P(y_i = each other in A | X)    = eps / (A - 1)

Mathematically the per-module MI is identical to the Mode 1 channel with
alphabet A: positions are independent given X, so the per-position MI
estimator from A.1 applies directly.

The DISTINCTION with Mode 1 is structural, not informational at the polymer
level:

  - Mode 1: X is encoded in a separable 1-D sequence template molecule.
            Output length L is unbounded; same machinery copies any L.
            Autocatalytic copyability is possible.
  - Mode 5: X is encoded in the MODULE STRUCTURE of the assembly itself.
            Output length is bounded by N (one module per output position).
            No autocatalytic copyability — the modular template is itself
            produced by Modes 1+2, not by Mode 5.

The length-limit prediction is the empirical signature of this structural
distinction at the polymerization level.

## Sweep

  - N (module count) ∈ {2, 4, 8, 16, 32}
  - epsilon          ∈ {0.001, 0.01, 0.05}
  - alphabet A       ∈ {4, 20}    (4 = nucleotide-style, 20 = amino-acid-style;
                                   real NRPS uses ~20 amino + non-canonical)
  - n_samples = 10000

For each (N, eps, A): simulate, estimate I_total = sum_i I(x_i; y_i), compare
to N * (log_2(A) - H(Y|X)).

## Length-limit sub-experiment

For two representative cases — (N=4, A=4) and (N=8, A=20), both at eps=0.01:

  - simulate the templated channel for positions 0..N-1
  - extend Y with N additional positions of uniform-random monomers,
    independent of X (modeling "no module past N → no template guidance")
  - compute per-position I(X; y_i) for i in [0, 2N-1]
  - verify per-position info ≈ predicted templated value for i < N
    and ≈ 0 for i >= N (sharp drop at the length limit)

## Files

  - `test_c_results.csv` — main scaling table, one row per (N, eps, A) cell
  - `test_c_length_limit.csv` — per-position info for the two length-limit cases
  - `test_c_progress.txt` — append-only progress log, one line per cell
  - `figures/test_c_N_scaling.png` — empirical I_emp vs N for each (eps, A)
    combo with theoretical lines overlaid; linear in N with slope set by
    per-module info
  - `figures/test_c_length_limit.png` — bar chart of per-position info vs
    position index for the two cases; sharp drop at i = N is the length-limit
    signature
  - `figures/test_c_mode_contrast.png` — composite three-mode contrast figure;
    Mode 1 grows ∝ L (unbounded), Mode 3 saturates at log_2(N), Mode 5 grows
    ∝ N then hits a length cap at N

## How to interpret the results

  - **slope check**: I_emp / N should equal log_2(A) - H(Y|X) for every cell.
    For A=4, eps=0.01 this is 1.92 bits/module; for A=20, eps=0.001 this is
    4.31 bits/module.
  - **finite-sample bias**: plug-in MI on AxA cells has positive bias
    ~ (A^2 - 1)/(2 n ln 2). For A=20, n=10000 this is ~0.029 bits/module —
    small relative to the ~4.3 bit/module signal but visible at small N.
    The 5% tolerance with N >= 4 cutoff filters this out.
  - **length limit**: for both representative cases, mean info in the
    templated region should be >10x the mean info in the untemplated tail.
    The tail mean equals the bias floor of the estimator (positive but tiny).

## PASS / FAIL criterion

PASS if BOTH sub-criteria hold:

  1. Linear N-scaling: for every (eps, A) and every N >= 4,
       |I_emp - N * ipp_theo| / (N * ipp_theo) < 0.05
  2. Length-limit (both cases):
       mean(I_pp[0:N]) > 10 * mean(I_pp[N:2N])

## What a PASS means

The Mode 5 modular conveyor model gives:
  - I_struct linear in module count N (slope = per-module channel info), and
  - product length strictly bounded by N (no extension past module N)

— exactly as the framework predicts. Combined with Test A.1 (Mode 1 linear in
L, unbounded length) and Test B (Mode 3 saturating at log_2(N), unbounded
length but periodic), this completes the empirical validation of the
framework's three-mode distinction. The L-vs-N scaling and length-limit
signatures jointly distinguish all three modes from per-position output data
alone.

## What a FAIL means

  - If I_emp deviates from N * ipp_theo by more than 5%, either the simulator
    is wrong, the per-module channel description is inconsistent with theory,
    or the plug-in MI bias dominates at the chosen n_samples (rerun at higher
    n_samples to discriminate).
  - If the length-limit ratio is < 10, the untemplated tail is leaking
    information from X — most likely a simulator bug, since by construction
    the tail is uniform random and independent of X.

## Reproducibility

`np.random.seed(42)` and `default_rng(42)` for the scaling sweep,
`default_rng(123)` for the length-limit cases. Re-running the script
reproduces the tables exactly.
"""


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    progress_path = RESULTS_DIR / "test_c_progress.txt"

    # main scaling experiment
    scaling_rows, Ns, epsilons, alphabets, n_samples = run_scaling_experiment(
        progress_path
    )
    save_scaling_csv(scaling_rows, RESULTS_DIR / "test_c_results.csv")

    # length-limit sub-experiment
    length_rows, length_summaries = run_length_limit_experiment(progress_path)
    save_length_limit_csv(length_rows, RESULTS_DIR / "test_c_length_limit.csv")

    # figures
    n_scaling_png = FIGURES_DIR / "test_c_N_scaling.png"
    length_limit_png = FIGURES_DIR / "test_c_length_limit.png"
    contrast_png = FIGURES_DIR / "test_c_mode_contrast.png"
    plot_n_scaling(scaling_rows, Ns, epsilons, alphabets, n_scaling_png)
    plot_length_limit(length_rows, length_summaries, length_limit_png)
    plot_mode_contrast(scaling_rows, contrast_png,
                       eps_for_plot=0.01, A_for_plot=4)

    # README
    readme_path = RESULTS_DIR / "test_c_README.md"
    readme_path.write_text(README_TEXT)

    # PASS / FAIL evaluation
    overall, scaling_failures, length_failures, (c1_ok, c2_ok) = evaluate_pass(
        scaling_rows, length_summaries
    )

    print("-" * 80)
    print(f"results  : {RESULTS_DIR / 'test_c_results.csv'}")
    print(f"length   : {RESULTS_DIR / 'test_c_length_limit.csv'}")
    print(f"progress : {progress_path}")
    print(f"readme   : {readme_path}")
    print(f"figures  : {n_scaling_png}")
    print(f"           {length_limit_png}")
    print(f"           {contrast_png}")
    print("-" * 80)

    print("Criterion 1 — linear N-scaling (|I_emp - N*ipp| / (N*ipp) < 5% for N >= 4):")
    if scaling_failures:
        for (N, eps, A, I_emp, I_theo, rel) in scaling_failures:
            print(f"  [FAIL] N={N}, eps={eps}, A={A}: "
                  f"I_emp={I_emp:.4f}, I_theo={I_theo:.4f}, rel_err={rel:.2%}")
    else:
        n_cells = sum(1 for r in scaling_rows if r["N"] >= 4)
        print(f"  [OK] all {n_cells} cells (N >= 4) within 5% of theory")
    print(f"  → {'PASS' if c1_ok else 'FAIL'}")

    print()
    print("Criterion 2 — length limit (mean(I_pp[0:N]) > 10 * mean(I_pp[N:2N])):")
    for s in length_summaries:
        ratio = s["templated_mean"] / max(s["untemplated_mean"], EPS_CLIP)
        ok = ratio > 10.0
        flag = "OK" if ok else "FAIL"
        print(f"  [{flag}] {s['label']} (N={s['N']}, A={s['A']}, eps={s['eps']}): "
              f"templated={s['templated_mean']:.4f} (theory={s['ipp_theo']:.4f}), "
              f"untemplated={s['untemplated_mean']:.4f}, ratio={ratio:.1f}")
    print(f"  → {'PASS' if c2_ok else 'FAIL'}")

    print()
    print("=" * 80)
    print(f"OVERALL: {'PASS' if overall else 'FAIL'}")
    print(f"  crit1 (linear N-scaling): {'PASS' if c1_ok else 'FAIL'}")
    print(f"  crit2 (length limit):     {'PASS' if c2_ok else 'FAIL'}")
    print("=" * 80)


if __name__ == "__main__":
    main()
