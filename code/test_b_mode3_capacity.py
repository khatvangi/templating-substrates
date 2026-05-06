"""
test B — Mode 3 cyclic active-site: information-capacity saturation prediction
==============================================================================

Templating Substrates Framework, Test B
---------------------------------------
Verifies the framework's central quantitative prediction for Mode 3 templating:

    I_struct(X; Y)  →  log_2(N)   bits, INDEPENDENT of polymer length L

This is in qualitative contrast to Mode 1 (test A.1), where I_struct grows
linearly in L. The L-scaling difference is the framework's central
diagnostic distinction between sequence-templated (Mode 1) and
structure-templated (Mode 3) information transfer.

Biological motivation: Drt3b (Rousseau et al., bioRxiv 2024) is a polymerase
that produces (AC)_n repetitive output without a sequence template, by
cycling through N=2 active-site conformational states. The framework places
this in Mode 3 and predicts its information capacity is exactly 1 bit
(= log_2(2)) regardless of polymer length.

Mode 3 channel
--------------
- Hidden "template" X = phi_0 ∈ {0, ..., N-1}, uniform, the cycle's start phase.
  H(X) = log_2(N) bits.
- At polymer position k, the active site is in state s_{(phi_0 + k) mod N}
  and is selective for substrate a_{(phi_0 + k) mod N}.
- P(Y_k = correct | X) = 1 - eps;  P(Y_k = each wrong | X) = eps / (N - 1).
- Y is a length-L sequence over alphabet {0, ..., N-1}.

Why we estimate I(X;Y) at the JOINT level (not per-position)
-------------------------------------------------------------
In Mode 1, Y positions are independent given X; per-position MI summed across
positions equals the joint MI exactly. In Mode 3, positions are NOT
independent — they are perfectly cyclically correlated. The per-position MI
estimator therefore underestimates joint I(X;Y). The cleanest, spec-aligned
approach is to treat the entire Y sequence as a single random object and use
the empirical joint distribution:

    I(X;Y) = H(Y) - H(Y|X)

For epsilon = 0, there are exactly N distinct Y trajectories (one per phase),
so H(Y) = log_2(N) and H(Y|X) = 0, giving I = log_2(N) exactly — independent
of L. For epsilon > 0, errors expand the support of Y and add randomness
given X; both H(Y) and H(Y|X) grow but their difference stays bounded
above by log_2(N).

Self-contained: no imports from other test scripts.
Uses only numpy, matplotlib, csv, pathlib.
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
# Mode 3 simulator
# ----------------------------------------------------------------------------
def simulate_mode3(N, L, epsilon, n_samples, rng):
    """
    generate n_samples independent (X, Y) draws from the Mode 3 channel.

    parameters
    ----------
    N         : number of cyclic active-site states (= alphabet size)
    L         : polymer length (positions in Y)
    epsilon   : misincorporation rate
    n_samples : number of independent (X, Y) draws
    rng       : numpy Generator (seeded by caller — must not be None)

    returns
    -------
    X : ndarray, shape (n_samples,), dtype int64, values in {0, ..., N-1}
        the starting cycle phase phi_0 for each trajectory.
    Y : ndarray, shape (n_samples, L), dtype int64, values in {0, ..., N-1}
        the polymer sequence emitted under that phase.
    """
    if N < 1:
        raise ValueError("N must be >= 1")
    if L < 1:
        raise ValueError("L must be >= 1")

    # phi_0 uniform over {0, ..., N-1}
    X = rng.integers(0, N, size=n_samples, dtype=np.int64)

    # intended monomer at each position: (phi_0 + k) mod N
    pos = np.arange(L, dtype=np.int64)[None, :]                # (1, L)
    Y_correct = (X[:, None] + pos) % N                         # (n_samples, L)

    if epsilon == 0.0 or N == 1:
        return X, Y_correct.astype(np.int64)

    # decide correct vs error per (sample, position)
    is_error = rng.random(size=(n_samples, L)) < epsilon

    # for error positions: pick uniformly among the N-1 wrong choices.
    # trick: pick offset in {0, ..., N-2}, then if offset >= correct value
    # shift up by 1 so we skip exactly the correct value, leaving the N-1
    # wrong choices uniformly.
    offset = rng.integers(0, N - 1, size=(n_samples, L), dtype=np.int64)
    Y_wrong = np.where(offset >= Y_correct, offset + 1, offset)

    Y = np.where(is_error, Y_wrong, Y_correct).astype(np.int64)
    return X, Y


# ----------------------------------------------------------------------------
# joint mutual information estimator
# ----------------------------------------------------------------------------
def _entropy_from_counts(counts):
    """plug-in entropy in bits from an array of nonneg integer counts."""
    total = counts.sum()
    if total == 0:
        return 0.0
    p = counts.astype(np.float64) / total
    nz = p[p > 0]
    return float(-(nz * np.log2(nz)).sum())


def _encode_rows_as_keys(Y):
    """
    encode each row of Y as a hashable bytes key. this avoids overflow risk
    that would arise from base-N integer encoding when N**L is large.

    returns a 1-D ndarray of bytes objects, one per row.
    """
    # contiguous, fixed-stride byte view of each row
    Yc = np.ascontiguousarray(Y, dtype=np.int64)
    n_rows = Yc.shape[0]
    row_bytes = Yc.view(np.uint8).reshape(n_rows, -1)
    # list of bytes is faster than np.unique on a 2-D array for our sizes
    return [bytes(row) for row in row_bytes]


def _entropy_from_keys(keys):
    """plug-in entropy in bits from an iterable of hashable keys."""
    if len(keys) == 0:
        return 0.0
    counts = {}
    for k in keys:
        counts[k] = counts.get(k, 0) + 1
    total = len(keys)
    p = np.array(list(counts.values()), dtype=np.float64) / total
    nz = p[p > 0]
    return float(-(nz * np.log2(nz)).sum())


def estimate_joint_information(X, Y, N):
    """
    estimate I(X;Y) = H(Y) - H(Y|X) on the FULL Y sequence using the
    empirical distribution over Y-strings.

    H(Y)        : plug-in entropy of the marginal distribution of Y rows.
    H(Y | X=x)  : plug-in entropy on the rows of Y restricted to X=x.
    H(Y | X)    : average of H(Y|X=x) weighted by P(X=x) ≈ count_x / n.

    This estimator is unbiased for epsilon=0 (support is just N strings) and
    has the standard plug-in positive bias for epsilon>0; we use it because
    the spec asks for direct empirical estimation at the sequence level, and
    the framework's PASS criteria for epsilon>0 are loose enough to absorb
    the small finite-sample bias at our n_samples calibration.

    parameters
    ----------
    X : (n_samples,) int, values in {0, ..., N-1}
    Y : (n_samples, L) int, values in {0, ..., N-1}
    N : alphabet / cyclic-state count

    returns
    -------
    I_emp, H_y, H_y_given_x : floats (bits)
    """
    n_samples = Y.shape[0]
    keys = _encode_rows_as_keys(Y)
    H_y = _entropy_from_keys(keys)

    H_y_given_x = 0.0
    X_arr = np.asarray(X, dtype=np.int64)
    counts_x = np.bincount(X_arr, minlength=N)
    for x in range(N):
        nx = counts_x[x]
        if nx == 0:
            continue
        # gather keys for this phase
        idx = np.flatnonzero(X_arr == x)
        sub_keys = [keys[i] for i in idx]
        H_y_given_x += (nx / n_samples) * _entropy_from_keys(sub_keys)

    I_emp = H_y - H_y_given_x
    return I_emp, H_y, H_y_given_x


# ----------------------------------------------------------------------------
# theoretical reference
# ----------------------------------------------------------------------------
def theoretical_capacity(N, epsilon):
    """
    framework prediction: I(X;Y) <= log_2(N) for all epsilon, with equality
    at epsilon=0. For epsilon>0 we don't have a clean closed form (the cycle
    couples positions), so we return log_2(N) as the upper-bound reference.
    """
    return float(np.log2(N))


# ----------------------------------------------------------------------------
# periodicity check (sanity diagnostic — not part of PASS criteria here,
# but used to validate the simulator)
# ----------------------------------------------------------------------------
def verify_periodicity_eps0(X, Y, N):
    """
    for epsilon=0 trajectories: assert that for every sample,
        Y[sample, k]  ==  (X[sample] + k) mod N
    returns the fraction of (sample, position) pairs that match (should be 1.0).
    """
    n_samples, L = Y.shape
    pos = np.arange(L, dtype=np.int64)[None, :]
    expected = (np.asarray(X, dtype=np.int64)[:, None] + pos) % N
    return float((Y == expected).mean())


# ----------------------------------------------------------------------------
# experiment driver — main sweep (N x epsilon at L = 5N)
# ----------------------------------------------------------------------------
def run_main_sweep(rng):
    """
    main sweep per spec:
      N ∈ {2, 3, 4, 5, 6, 8, 10}
      epsilon ∈ {0.0, 0.01, 0.05}
      L = N * 5
      n_samples = 20_000
    """
    Ns = [2, 3, 4, 5, 6, 8, 10]
    epsilons = [0.0, 0.01, 0.05]
    n_samples = 20_000

    rows = []
    print(f"\n{'='*88}")
    print("Test B — main sweep: I(X;Y) at L = 5N for varying N, epsilon")
    print(f"{'='*88}")
    print(f"{'N':>3s} {'eps':>6s} {'L':>4s} {'n':>6s} "
          f"{'H(Y)':>8s} {'H(Y|X)':>8s} {'I_emp':>8s} "
          f"{'log2N':>8s} {'I/log2N':>8s} {'time':>6s}")
    print("-" * 88)

    for N in Ns:
        log2N = theoretical_capacity(N, 0.0)
        L = 5 * N
        for eps in epsilons:
            t0 = time.time()
            X, Y = simulate_mode3(N, L, eps, n_samples, rng)
            I_emp, H_y, H_yx = estimate_joint_information(X, Y, N)
            elapsed = time.time() - t0

            # for epsilon=0 verify perfect periodicity as a sanity check
            if eps == 0.0:
                frac_match = verify_periodicity_eps0(X, Y, N)
                if abs(frac_match - 1.0) > EPS_CLIP:
                    raise RuntimeError(
                        f"periodicity check failed for N={N}, L={L}: "
                        f"frac_match={frac_match}"
                    )

            ratio = I_emp / log2N if log2N > 0 else float("nan")
            row = {
                "N": N,
                "L": L,
                "epsilon": eps,
                "n_samples": n_samples,
                "H_Y": H_y,
                "H_Y_given_X": H_yx,
                "I_empirical": I_emp,
                "I_theoretical_eps0": log2N,
                "ratio_to_logN": ratio,
            }
            rows.append(row)
            print(f"{N:3d} {eps:6.3f} {L:4d} {n_samples:6d} "
                  f"{H_y:8.4f} {H_yx:8.4f} {I_emp:8.4f} "
                  f"{log2N:8.4f} {ratio:8.4f} {elapsed:5.1f}s")
    return rows, Ns, epsilons, n_samples


# ----------------------------------------------------------------------------
# L-scaling sub-experiment (fixed N, varying L)
# ----------------------------------------------------------------------------
def run_l_scaling(rng, N=4, L_multipliers=(2, 3, 5, 10), epsilons=(0.0, 0.01, 0.05),
                  n_samples=20_000):
    """
    KEY plot: at fixed N, varying L should NOT change I(X;Y).

    returns a list of dict rows: {N, L, epsilon, n_samples, H_Y, H_Y_given_X,
    I_empirical, I_theoretical_eps0, ratio_to_logN}.
    """
    log2N = theoretical_capacity(N, 0.0)
    rows = []
    print(f"\n{'='*88}")
    print(f"Test B — L-scaling sub-test at fixed N = {N}")
    print(f"{'='*88}")
    print(f"{'N':>3s} {'eps':>6s} {'L':>4s} {'n':>6s} "
          f"{'I_emp':>8s} {'log2N':>8s} {'I/log2N':>8s} {'time':>6s}")
    print("-" * 88)

    for eps in epsilons:
        for mult in L_multipliers:
            L = N * mult
            t0 = time.time()
            X, Y = simulate_mode3(N, L, eps, n_samples, rng)
            I_emp, H_y, H_yx = estimate_joint_information(X, Y, N)
            elapsed = time.time() - t0
            ratio = I_emp / log2N
            row = {
                "N": N,
                "L": L,
                "epsilon": eps,
                "n_samples": n_samples,
                "H_Y": H_y,
                "H_Y_given_X": H_yx,
                "I_empirical": I_emp,
                "I_theoretical_eps0": log2N,
                "ratio_to_logN": ratio,
            }
            rows.append(row)
            print(f"{N:3d} {eps:6.3f} {L:4d} {n_samples:6d} "
                  f"{I_emp:8.4f} {log2N:8.4f} {ratio:8.4f} {elapsed:5.1f}s")
    return rows


# ----------------------------------------------------------------------------
# I/O — CSV
# ----------------------------------------------------------------------------
CSV_FIELDS = [
    "N",
    "L",
    "epsilon",
    "n_samples",
    "H_Y",
    "H_Y_given_X",
    "I_empirical",
    "I_theoretical_eps0",
    "ratio_to_logN",
]


def save_csv(rows, path):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r[k] for k in CSV_FIELDS})


# ----------------------------------------------------------------------------
# plotting
# ----------------------------------------------------------------------------
def plot_capacity(main_rows, Ns, epsilons, path):
    """
    figures/test_b_capacity.png
    scatter of I_empirical vs N for each epsilon. dashed line = log_2(N) for
    each N, showing saturation.
    """
    fig, ax = plt.subplots(figsize=(8.0, 5.5))

    # log_2(N) reference line
    Ns_arr = np.array(Ns, dtype=float)
    ax.plot(Ns_arr, np.log2(Ns_arr), "--", color="black", linewidth=1.4,
            alpha=0.8, label=r"theory: $\log_2 N$")

    eps_colors = {0.0: "C0", 0.01: "C1", 0.05: "C2"}
    eps_markers = {0.0: "o", 0.01: "s", 0.05: "^"}

    for eps in epsilons:
        xs, ys = [], []
        for r in main_rows:
            if r["epsilon"] == eps:
                xs.append(r["N"])
                ys.append(r["I_empirical"])
        order = np.argsort(xs)
        xs_a = np.array(xs)[order]
        ys_a = np.array(ys)[order]
        ax.plot(xs_a, ys_a, eps_markers.get(eps, "o"),
                color=eps_colors.get(eps, None),
                markersize=8, linewidth=0,
                label=f"empirical (eps={eps})")

    ax.set_xlabel("N (number of cyclic active-site states)")
    ax.set_ylabel(r"$I(X;Y)$  [bits]")
    ax.set_title(
        r"Test B — Mode 3 capacity saturates at $\log_2 N$  (L = 5N, n = 20k)"
    )
    ax.set_xticks(Ns)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_l_scaling(l_rows, N, path):
    """
    figures/test_b_l_scaling.png
    one line per epsilon, x = L, y = I_empirical. should be flat.
    horizontal dashed line at log_2(N).
    """
    fig, ax = plt.subplots(figsize=(8.0, 5.5))

    log2N = float(np.log2(N))
    ax.axhline(log2N, color="black", linestyle="--", linewidth=1.4,
               alpha=0.8, label=fr"theory: $\log_2 {N} = {log2N:.3f}$")

    eps_colors = {0.0: "C0", 0.01: "C1", 0.05: "C2"}
    eps_markers = {0.0: "o", 0.01: "s", 0.05: "^"}

    epsilons_seen = sorted({r["epsilon"] for r in l_rows})
    for eps in epsilons_seen:
        xs, ys = [], []
        for r in l_rows:
            if r["epsilon"] == eps:
                xs.append(r["L"])
                ys.append(r["I_empirical"])
        order = np.argsort(xs)
        xs_a = np.array(xs)[order]
        ys_a = np.array(ys)[order]
        ax.plot(xs_a, ys_a, "-" + eps_markers.get(eps, "o"),
                color=eps_colors.get(eps, None),
                markersize=8, linewidth=1.6,
                label=f"empirical (eps={eps})")

    ax.set_xlabel("polymer length L")
    ax.set_ylabel(r"$I(X;Y)$  [bits]")
    ax.set_title(
        rf"Test B — L-scaling at fixed N = {N}: $I$ is L-INDEPENDENT (Mode 3 signature)"
    )
    ax.grid(True, alpha=0.3)
    # widen y-axis a bit to make the flatness visually obvious
    ymin = min([r["I_empirical"] for r in l_rows] + [log2N]) - 0.15
    ymax = max([r["I_empirical"] for r in l_rows] + [log2N]) + 0.15
    ax.set_ylim(ymin, ymax)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ----------------------------------------------------------------------------
# PASS / FAIL evaluation
# ----------------------------------------------------------------------------
def evaluate_pass(main_rows, l_rows):
    """
    PASS criteria per spec:
      (a) for epsilon=0: |I_empirical - log_2(N)| < 0.05 bits, all N tested.
      (b) for epsilon>0: I_empirical < log_2(N) + 0.05 bits AND
                         I_empirical > log_2(N) * 0.5  (some info, but bounded).
      (c) L-scaling sub-test: I at L=N*10 differs from I at L=N*2 by < 10%
          (per epsilon, in the l_rows).
    """
    detail_a = []
    detail_b = []
    detail_c = []

    for r in main_rows:
        log2N = float(np.log2(r["N"]))
        I = r["I_empirical"]
        if r["epsilon"] == 0.0:
            ok_a = abs(I - log2N) < 0.05
            detail_a.append((r["N"], r["epsilon"], r["L"], I, log2N,
                             abs(I - log2N), ok_a))
        else:
            upper_ok = I < log2N + 0.05
            lower_ok = I > 0.5 * log2N
            ok_b = upper_ok and lower_ok
            detail_b.append((r["N"], r["epsilon"], r["L"], I, log2N,
                             upper_ok, lower_ok, ok_b))

    # criterion (c): L-scaling — for each epsilon, ratio at L=N*10 vs L=N*2
    by_eps = {}
    for r in l_rows:
        by_eps.setdefault(r["epsilon"], {})[r["L"]] = r["I_empirical"]
    N_used = l_rows[0]["N"] if l_rows else None
    for eps, by_L in by_eps.items():
        L_low = N_used * 2
        L_high = N_used * 10
        if L_low not in by_L or L_high not in by_L:
            continue
        I_low = by_L[L_low]
        I_high = by_L[L_high]
        denom = max(abs(I_low), abs(I_high), EPS_CLIP)
        rel_diff = abs(I_high - I_low) / denom
        ok_c = rel_diff < 0.10
        detail_c.append((N_used, eps, L_low, L_high, I_low, I_high,
                         rel_diff, ok_c))

    crit_a_ok = all(t[-1] for t in detail_a)
    crit_b_ok = all(t[-1] for t in detail_b)
    crit_c_ok = all(t[-1] for t in detail_c)
    overall = crit_a_ok and crit_b_ok and crit_c_ok
    return overall, (crit_a_ok, crit_b_ok, crit_c_ok), {
        "a": detail_a, "b": detail_b, "c": detail_c,
    }


def print_pass_report(overall, crits, details):
    crit_a_ok, crit_b_ok, crit_c_ok = crits
    print()
    print("=" * 88)
    print("PASS / FAIL evaluation")
    print("=" * 88)

    print()
    print("(a) eps=0 saturation: |I_emp - log_2(N)| < 0.05 bits")
    for N, eps, L, I, log2N, abs_err, ok in details["a"]:
        flag = "OK  " if ok else "FAIL"
        print(f"   [{flag}] N={N:2d}, eps={eps:.3f}, L={L:3d}: "
              f"I={I:.4f}, log2N={log2N:.4f}, |diff|={abs_err:.5f}")
    print(f"  → {'PASS' if crit_a_ok else 'FAIL'}")

    print()
    print("(b) eps>0 bounds: I_emp < log_2(N) + 0.05 AND I_emp > log_2(N)/2")
    for N, eps, L, I, log2N, upper_ok, lower_ok, ok in details["b"]:
        flag = "OK  " if ok else "FAIL"
        u = "✓" if upper_ok else "✗"
        l = "✓" if lower_ok else "✗"
        print(f"   [{flag}] N={N:2d}, eps={eps:.3f}, L={L:3d}: "
              f"I={I:.4f}, log2N={log2N:.4f}  upper{u} lower{l}")
    print(f"  → {'PASS' if crit_b_ok else 'FAIL'}")

    print()
    print("(c) L-scaling: |I(L=10N) - I(L=2N)| / I < 10%  (Mode 3 signature)")
    for N, eps, L_lo, L_hi, I_lo, I_hi, rel, ok in details["c"]:
        flag = "OK  " if ok else "FAIL"
        print(f"   [{flag}] N={N}, eps={eps:.3f}: I(L={L_lo})={I_lo:.4f}, "
              f"I(L={L_hi})={I_hi:.4f}, rel_diff={rel:.2%}")
    print(f"  → {'PASS' if crit_c_ok else 'FAIL'}")

    print()
    print("=" * 88)
    print(f"OVERALL: {'PASS' if overall else 'FAIL'}")
    print(f"  (a) eps=0 saturation:        {'PASS' if crit_a_ok else 'FAIL'}")
    print(f"  (b) eps>0 bounds:            {'PASS' if crit_b_ok else 'FAIL'}")
    print(f"  (c) L-independence:          {'PASS' if crit_c_ok else 'FAIL'}")
    print("=" * 88)


# ----------------------------------------------------------------------------
# README
# ----------------------------------------------------------------------------
README_TEXT = """\
# Test B — Mode 3 Cyclic Active-Site Information Capacity

## What this test does

Validates the framework's central distinguishing prediction for Mode 3
templating: a system whose information sits in the *structure* of an N-state
cyclic active site (rather than in a sequence template) has a transferred
information that saturates at log_2(N) bits, INDEPENDENT of the polymer
output length L.

Compared with Test A.1, where Mode 1 gives `I_struct ∝ L`, this is the
framework's central diagnostic distinction:

| mode | scaling of I(X;Y) with L            |
|------|-------------------------------------|
| 1    | linear in L                          |
| 3    | bounded by log_2(N), L-independent   |

## Biological motivation

Drt3b (Rousseau et al., bioRxiv 2024) is a polymerase that produces (AC)_n
without a sequence template, by cycling between N=2 active-site
conformational states. The framework predicts its information content is
exactly 1 bit (= log_2(2)), regardless of how long the (AC)_n polymer is.
Test B with N=2 reproduces this prediction; the test with N ∈ {3,...,10}
generalizes it.

## Channel definition

- Hidden phase X = phi_0 ∈ {0, ..., N-1}, uniform; H(X) = log_2(N) bits.
- Position k intended monomer = (phi_0 + k) mod N.
- P(Y_k = correct | X) = 1 - eps;  P(Y_k = each wrong | X) = eps / (N - 1).

Determinism note (epsilon=0): given phi_0, the entire Y is fixed. So Y has
exactly N possible values across the population. H(Y) = log_2(N), H(Y|X) = 0,
I(X;Y) = log_2(N) — and this is independent of L.

## Estimator

We estimate the joint mutual information at the SEQUENCE level:

    I(X; Y) = H(Y) - H(Y|X)

with both terms computed as plug-in entropies on empirical distributions:
- H(Y): hash each Y row to a bytes key; entropy of the key distribution.
- H(Y|X): for each phase x, restrict samples to X==x, compute H(Y) on that
  subset, weight by P(X=x), sum.

This differs from Test A.1's per-position MI estimator because positions in
Mode 3 are NOT independent given X — they are perfectly cyclically
correlated. Per-position MI summed across positions would underestimate the
true joint information here.

For epsilon=0 the estimator is exact: only N distinct Y rows exist, so the
plug-in entropy is precise. For epsilon>0 the estimator has a small positive
plug-in bias which is negligible at the n_samples used here for the loose
upper-bound criterion.

## Sweep

Main capacity sweep (figures/test_b_capacity.png):
- N ∈ {2, 3, 4, 5, 6, 8, 10}
- epsilon ∈ {0.0, 0.01, 0.05}
- L = 5N (5 full cycles, enough to span periodicity but small enough that
  empirical estimation is reliable)
- n_samples = 20,000

L-scaling sub-test (figures/test_b_l_scaling.png) — the KEY diagnostic:
- N = 4 fixed
- L ∈ {N*2, N*3, N*5, N*10} = {8, 12, 20, 40}
- epsilon ∈ {0.0, 0.01, 0.05}
- n_samples = 20,000
- shows that I(X;Y) does NOT scale with L

## Files

- `test_b_results.csv` — every cell from BOTH the main sweep and the
  L-scaling sub-test. Columns: N, L, epsilon, n_samples, H_Y, H_Y_given_X,
  I_empirical, I_theoretical_eps0, ratio_to_logN.
- `figures/test_b_capacity.png` — I_empirical vs N for each epsilon, with
  dashed log_2(N) reference line. The points fall on or below the line.
- `figures/test_b_l_scaling.png` — at fixed N=4, I_empirical vs L for each
  epsilon. The lines should be flat (this is the Mode 3 signature).

## PASS / FAIL criteria

PASS if all three sub-criteria hold:

  (a) Low-noise saturation. For epsilon = 0:
      |I_empirical - log_2(N)| < 0.05 bits, for all N tested.

  (b) Noisy bounds. For epsilon > 0:
      I_empirical  <  log_2(N) + 0.05  bits   (cannot exceed the bound)
      I_empirical  >  0.5 * log_2(N)          (some structural info remains)

  (c) L-independence. In the L-scaling sub-test (N=4):
      |I(L=10N) - I(L=2N)| / I  <  10%   (per epsilon)

## What a PASS means

- The Mode 3 cyclic-active-site model carries exactly log_2(N) bits of
  information about the start phase, regardless of how long the polymer
  output is.
- Combined with Test A.1 (I ∝ L for Mode 1), this confirms the framework's
  two-mode distinction has empirical bite: the L-scaling slope distinguishes
  sequence-templated information (Mode 1) from structure-templated
  information (Mode 3).
- For Drt3b specifically (N=2): I_struct = 1 bit, regardless of polymer
  length. This is the framework's prediction for a real biological system.

## What a FAIL means

- If I_empirical grows with L for Mode 3, the saturation prediction is
  refuted — the framework would need an additional information channel
  beyond the cyclic state count.
- If I_empirical exceeds log_2(N) significantly, the cyclic-state count is
  NOT the actual information capacity of the model.
- If I_empirical at epsilon=0 deviates from log_2(N), the simulator or the
  estimator has an error.

## Reproducibility

Random seed: `np.random.seed(42)` and `np.random.default_rng(42)`. Re-running
the script reproduces the table exactly.
"""


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main():
    np.random.seed(42)
    rng = np.random.default_rng(42)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # main sweep at L = 5N
    main_rows, Ns, epsilons, n_samples = run_main_sweep(rng)

    # L-scaling sub-test at fixed N=4
    l_rows = run_l_scaling(
        rng,
        N=4,
        L_multipliers=(2, 3, 5, 10),
        epsilons=(0.0, 0.01, 0.05),
        n_samples=20_000,
    )

    # write CSV (combined: main sweep then L-scaling rows)
    csv_path = RESULTS_DIR / "test_b_results.csv"
    save_csv(main_rows + l_rows, csv_path)

    # figures
    capacity_png = FIGURES_DIR / "test_b_capacity.png"
    l_scaling_png = FIGURES_DIR / "test_b_l_scaling.png"
    plot_capacity(main_rows, Ns, epsilons, capacity_png)
    plot_l_scaling(l_rows, N=4, path=l_scaling_png)

    # README
    readme_path = RESULTS_DIR / "test_b_README.md"
    readme_path.write_text(README_TEXT)

    print()
    print(f"results: {csv_path}")
    print(f"figures: {capacity_png}")
    print(f"         {l_scaling_png}")
    print(f"readme : {readme_path}")

    # PASS / FAIL
    overall, crits, details = evaluate_pass(main_rows, l_rows)
    print_pass_report(overall, crits, details)


if __name__ == "__main__":
    main()
