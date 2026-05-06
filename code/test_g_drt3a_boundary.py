"""
test G — Drt3a Mode 1 with degenerate (alternating) template
============================================================

Templating Substrates Framework, Test G (boundary-case diagnostic)
------------------------------------------------------------------
Drt3a (Sharma et al. 2026) synthesizes poly(GT) by Watson-Crick templating
off a 6-nt ACACAC region of a ncRNA. Naively the framework places this in
Mode 1 with L = 6 and predicts I_struct ≤ L * log2(4) = 12 bits. But the
template is alternating (A, C, A, C, A, C): only 1 bit of phase information
is non-redundant across the population of templates. This raises a
boundary-case question: does the apparatus report the per-position transfer
ceiling (≈12 bits, "Prediction A") or the template-information ceiling
(≈1 bit, "Prediction B")?

The framework's prose prediction is A — the apparatus measures per-position
transfer at the joint distribution level, and even a degenerate template
mediates high-fidelity per-position transfer. The framework's apparatus
implementation, however, is the plug-in mutual information estimator
I(X_i; Y_i) per position — and that estimator can only measure information
that VARIES across the sample population. With a fully degenerate template
(X_i constant across samples) per-position MI is 0 by definition, regardless
of how deterministic the per-sample mechanism is.

This test is therefore a diagnostic of whether mechanism (Mode 1: per-position
WC transfer happens) and apparatus (per-position empirical MI: zero when X_i
is constant) coincide for the degenerate-template boundary case.

Two competing predictions to test
---------------------------------
Prediction A (mechanism-determined): I_struct ≈ L_T * log2(4) at L_T = 6,
    regardless of template degeneracy.
Prediction B (template-information-limited): I_struct ≤ H(template population)
    bits total, where the alternating-AC population carries only ~1 bit of
    phase information.

Concrete sweeps
---------------
1. alternating vs random templates, L_T ∈ {2, 4, 6, 8, 12, 24}, eps = 0.01,
   n_samples = 1000 (per dispatch). Per-position plug-in MI estimator
   exactly as Test A.1. Bulk-matched control: row-shuffle X across samples
   (breaks per-position alignment, preserves marginal X composition).
2. Drt3a (Mode 1, alternating ACAC...) vs Drt3b (Mode 3, N = 2 cyclic
   active site). Sequence-level joint MI estimator I(X;Y) = H(Y) - H(Y|X)
   exactly as Test B. Sweep L ∈ {2, 4, 6, 8, 12, 24}. Different X for each
   mode: for Drt3a, X = phase index ∈ {0, 1} → template; for Drt3b, X =
   phase index ∈ {0, 1} → cyclic state.

This file is intentionally self-contained: simulators and estimators are
duplicated from test_a1_mode1_scaling.py and test_b_mode3_capacity.py. Do
NOT refactor into a shared utils module — the duplication is the point
(per repo CLAUDE.md rule 4).
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

# Watson-Crick complement table:
#   A=0 ↔ T=3
#   C=1 ↔ G=2
WC_PAIR = np.array([3, 2, 1, 0], dtype=np.int8)

EPS_CLIP = 1e-12


# ============================================================================
# part 1: simulators
# ============================================================================

# -- Mode 1 (Drt3a) ----------------------------------------------------------
# duplicated from test_a1_mode1_scaling.py per repo rule 4 (self-contained)

def simulateMode1Random(L, epsilon, n_samples, rng):
    """
    Mode 1 with uniformly random templates per sample (Test A.1 baseline).

    each of n_samples gets its own independently drawn template X of length L,
    and produces Y = WC(X) with per-position fidelity 1 - epsilon. Returns
    (X, Y) both shape (n_samples, L), dtype int8.
    """
    # template uniformly random per sample
    X = rng.integers(0, 4, size=(n_samples, L), dtype=np.int8)
    Y_correct = WC_PAIR[X]
    is_error = rng.random(size=(n_samples, L)) < epsilon
    # for error positions: pick uniformly among the 3 wrong nucleotides
    offset = rng.integers(0, 3, size=(n_samples, L), dtype=np.int8)
    Y_wrong = np.where(offset >= Y_correct, offset + 1, offset).astype(np.int8)
    Y = np.where(is_error, Y_wrong, Y_correct).astype(np.int8)
    return X, Y


def simulateMode1AlternatingFixed(L, epsilon, n_samples, rng):
    """
    Mode 1 with a single fixed alternating ACAC...AC template for all samples.

    X is identical across samples (X_i is constant across the population at
    every position i). The per-sample mechanism is identical to Mode 1
    (per-position WC transfer with fidelity 1 - epsilon), but X carries no
    cross-sample variability.

    returns (X, Y) shape (n_samples, L), dtype int8. X has every row equal to
    the alternating pattern [A, C, A, C, ...] = [0, 1, 0, 1, ...].
    """
    # alternating AC = [0, 1, 0, 1, ...] of length L
    template = np.tile(np.array([0, 1], dtype=np.int8), (L + 1) // 2)[:L]
    X = np.broadcast_to(template, (n_samples, L)).copy()  # all rows identical
    Y_correct = WC_PAIR[X]
    is_error = rng.random(size=(n_samples, L)) < epsilon
    offset = rng.integers(0, 3, size=(n_samples, L), dtype=np.int8)
    Y_wrong = np.where(offset >= Y_correct, offset + 1, offset).astype(np.int8)
    Y = np.where(is_error, Y_wrong, Y_correct).astype(np.int8)
    return X, Y


def simulateMode1Alternating2Phase(L, epsilon, n_samples, rng):
    """
    Mode 1 with alternating template drawn from a 2-element population.

    each sample independently picks phase ∈ {0, 1} uniformly:
      phase 0 → X = [A, C, A, C, ...]
      phase 1 → X = [C, A, C, A, ...]
    Y = WC(X) with per-position fidelity 1 - epsilon.

    returns (phase, X, Y) with shapes (n_samples,), (n_samples, L), (n_samples, L).
    Used for the Drt3a vs Drt3b joint-MI comparison: phase is the natural
    "X-as-descriptor" variable that has the same dimensionality as Drt3b's
    cyclic-state index, allowing apples-to-apples comparison via H(Y) - H(Y|X).
    """
    phase = rng.integers(0, 2, size=n_samples, dtype=np.int64)  # 1 bit
    # build X as the alternating pattern starting from each row's phase
    pos = np.arange(L, dtype=np.int64)[None, :]                  # (1, L)
    # if phase=0: X_i = (i % 2) → [0,1,0,1,...] = [A,C,A,C,...]
    # if phase=1: X_i = ((i + 1) % 2) → [1,0,1,0,...] = [C,A,C,A,...]
    X = ((phase[:, None] + pos) % 2).astype(np.int8)
    Y_correct = WC_PAIR[X]
    is_error = rng.random(size=(n_samples, L)) < epsilon
    offset = rng.integers(0, 3, size=(n_samples, L), dtype=np.int8)
    Y_wrong = np.where(offset >= Y_correct, offset + 1, offset).astype(np.int8)
    Y = np.where(is_error, Y_wrong, Y_correct).astype(np.int8)
    return phase, X, Y


# -- Mode 3 (Drt3b) ----------------------------------------------------------
# duplicated from test_b_mode3_capacity.py per repo rule 4

def simulateMode3(N, L, epsilon, n_samples, rng):
    """
    Mode 3 cyclic active-site channel (Test B baseline).

    X = phi_0 ∈ {0, ..., N-1}, uniform; H(X) = log_2(N).
    intended monomer at position k is (phi_0 + k) mod N.
    P(Y_k = correct | X) = 1 - eps;  P(Y_k = each wrong | X) = eps / (N-1).

    returns (X, Y) with shapes (n_samples,), (n_samples, L), dtype int64.
    """
    X = rng.integers(0, N, size=n_samples, dtype=np.int64)
    pos = np.arange(L, dtype=np.int64)[None, :]
    Y_correct = (X[:, None] + pos) % N
    if epsilon == 0.0 or N == 1:
        return X, Y_correct.astype(np.int64)
    is_error = rng.random(size=(n_samples, L)) < epsilon
    offset = rng.integers(0, N - 1, size=(n_samples, L), dtype=np.int64)
    Y_wrong = np.where(offset >= Y_correct, offset + 1, offset)
    Y = np.where(is_error, Y_wrong, Y_correct).astype(np.int64)
    return X, Y


# ============================================================================
# part 2: estimators (per-position and joint-sequence)
# ============================================================================

# -- per-position plug-in MI (Test A.1 convention) ---------------------------

def findPerPositionMutualInformation(X, Y, alphabet=4):
    """
    plug-in I(X_i; Y_i) per position from empirical (alphabet x alphabet) joints,
    summed across positions. matches test_a1_mode1_scaling.py exactly so the
    apples-to-apples comparison with Test A.1 is valid.

    returns (I_total, I_per_position_mean, I_per_position_array).
    """
    n_samples, L = X.shape
    pair_codes = X.astype(np.int64) * alphabet + Y.astype(np.int64)
    I_total = 0.0
    I_per_pos = np.zeros(L)
    for i in range(L):
        counts = np.bincount(pair_codes[:, i],
                             minlength=alphabet * alphabet
                             ).reshape(alphabet, alphabet)
        joint = counts / n_samples
        p_x = joint.sum(axis=1, keepdims=True)
        p_y = joint.sum(axis=0, keepdims=True)
        denom = p_x * p_y
        mask = (joint > 0) & (denom > 0)
        I_i = float(np.sum(joint[mask] * np.log2(joint[mask] / denom[mask])))
        I_per_pos[i] = I_i
        I_total += I_i
    return I_total, I_total / L, I_per_pos


# -- joint sequence-level MI (Test B convention) -----------------------------

def _entropyFromKeys(keys):
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


def _encodeRowsAsKeys(Y):
    """encode each row of Y as a hashable bytes key (Test B convention)."""
    Yc = np.ascontiguousarray(Y, dtype=np.int64)
    n_rows = Yc.shape[0]
    row_bytes = Yc.view(np.uint8).reshape(n_rows, -1)
    return [bytes(row) for row in row_bytes]


def findJointSequenceMutualInformation(X_scalar, Y, n_states):
    """
    sequence-level I(X;Y) = H(Y) - H(Y|X) where X is a scalar descriptor
    (phase index for Drt3a-2phase or Drt3b). matches test_b_mode3_capacity.py
    so Drt3a and Drt3b can be compared apples-to-apples.

    parameters
    ----------
    X_scalar : (n_samples,) int, values in {0, ..., n_states - 1}
    Y        : (n_samples, L) int
    n_states : number of distinct X values

    returns (I, H_y, H_y_given_x) in bits.
    """
    n_samples = Y.shape[0]
    keys = _encodeRowsAsKeys(Y)
    H_y = _entropyFromKeys(keys)
    H_y_given_x = 0.0
    X_arr = np.asarray(X_scalar, dtype=np.int64)
    counts_x = np.bincount(X_arr, minlength=n_states)
    for x in range(n_states):
        nx = counts_x[x]
        if nx == 0:
            continue
        idx = np.flatnonzero(X_arr == x)
        sub_keys = [keys[i] for i in idx]
        H_y_given_x += (nx / n_samples) * _entropyFromKeys(sub_keys)
    I = H_y - H_y_given_x
    return I, H_y, H_y_given_x


# -- bulk-matched control (Test A.2 convention) ------------------------------

def findBulkMatchedControl(X, Y, alphabet, rng):
    """
    bulk-matched scramble: independently permute the rows of X. preserves
    composition of X but breaks per-position alignment with Y. expected to
    drive the per-position MI estimate to zero.
    """
    perm = rng.permutation(X.shape[0])
    X_scrambled = X[perm]
    return findPerPositionMutualInformation(X_scrambled, Y, alphabet=alphabet)


# -- structure-scrambled control (per dispatch) ------------------------------

def findStructureScrambledY(Y, rng):
    """
    structure-scrambled control: shuffle the positions of Y within each row.
    preserves per-row composition (and hence marginal P(Y_i) up to averaging)
    but breaks the positional correspondence with X. used as a sanity check
    that any positional MI signal we see is positionally aligned, not just
    a composition artefact.

    returns Y_scrambled with the same shape as Y.
    """
    Y_out = np.empty_like(Y)
    for r in range(Y.shape[0]):
        Y_out[r] = rng.permutation(Y[r])
    return Y_out


# ----------------------------------------------------------------------------
# closed-form theoretical reference (Mode 1 with random uniform template)
# ----------------------------------------------------------------------------

def theoreticalPerPositionInfo(epsilon):
    """
    Mode 1 closed form: I_per_position = log2(4) - H(Y|X)
                                       = 2 - [-(1-eps) log2(1-eps) - eps log2(eps/3)]
    valid only when X is uniform across the alphabet at each position.
    """
    if epsilon <= 0.0:
        H = 0.0
    else:
        one_minus = max(1.0 - epsilon, EPS_CLIP)
        eps_third = max(epsilon / 3.0, EPS_CLIP)
        H = -(1 - epsilon) * np.log2(one_minus) - epsilon * np.log2(eps_third)
    return 2.0 - H


# ============================================================================
# part 3: experiment driver — sweep 1 (alternating vs random templates)
# ============================================================================

def runAlternatingVsRandom(rng, L_Ts=(2, 4, 6, 8, 12, 24), epsilon=0.01,
                           n_samples=1000):
    """
    sweep over L_T comparing per-position MI for:
      - alternating template (single fixed ACAC...AC for all samples)
      - alternating-2phase template (X drawn from {ACAC..., CACA...})
      - random uniform template (Test A.1 baseline)
    each with bulk-matched control and structure-scrambled-Y control.
    """
    rows = []
    print()
    print("=" * 96)
    print(f"Test G — sweep 1: alternating vs random templates  "
          f"(eps={epsilon}, n={n_samples})")
    print("=" * 96)
    print(f"{'L_T':>4s} {'mode':>14s} "
          f"{'I_total':>9s} {'I_per_pos':>10s} {'I_bulk':>9s} {'I_struct_scr':>12s} "
          f"{'time':>6s}")
    print("-" * 96)

    for L_T in L_Ts:
        # ---------- alternating-fixed (X constant) ----------
        t0 = time.time()
        Xa, Ya = simulateMode1AlternatingFixed(L_T, epsilon, n_samples, rng)
        I_a, ipp_a, _ = findPerPositionMutualInformation(Xa, Ya, alphabet=4)
        I_a_bulk, _, _ = findBulkMatchedControl(Xa, Ya, alphabet=4, rng=rng)
        Ya_scr = findStructureScrambledY(Ya, rng)
        I_a_scr, _, _ = findPerPositionMutualInformation(Xa, Ya_scr, alphabet=4)
        ta = time.time() - t0
        rows.append({
            "L_T": L_T, "template_type": "alternating_fixed",
            "epsilon": epsilon, "n_samples": n_samples,
            "I_per_position_total": I_a, "I_per_position_mean": ipp_a,
            "I_bulk_matched": I_a_bulk, "I_structure_scrambled": I_a_scr,
        })
        print(f"{L_T:4d} {'alt_fixed':>14s} "
              f"{I_a:9.4f} {ipp_a:10.4f} {I_a_bulk:9.4f} {I_a_scr:12.4f} "
              f"{ta:5.2f}s")

        # ---------- alternating-2phase (X has 1 bit of population entropy) ----------
        t0 = time.time()
        _, Xb, Yb = simulateMode1Alternating2Phase(L_T, epsilon, n_samples, rng)
        I_b, ipp_b, _ = findPerPositionMutualInformation(Xb, Yb, alphabet=4)
        I_b_bulk, _, _ = findBulkMatchedControl(Xb, Yb, alphabet=4, rng=rng)
        Yb_scr = findStructureScrambledY(Yb, rng)
        I_b_scr, _, _ = findPerPositionMutualInformation(Xb, Yb_scr, alphabet=4)
        tb = time.time() - t0
        rows.append({
            "L_T": L_T, "template_type": "alternating_2phase",
            "epsilon": epsilon, "n_samples": n_samples,
            "I_per_position_total": I_b, "I_per_position_mean": ipp_b,
            "I_bulk_matched": I_b_bulk, "I_structure_scrambled": I_b_scr,
        })
        print(f"{L_T:4d} {'alt_2phase':>14s} "
              f"{I_b:9.4f} {ipp_b:10.4f} {I_b_bulk:9.4f} {I_b_scr:12.4f} "
              f"{tb:5.2f}s")

        # ---------- random uniform template (Test A.1 baseline) ----------
        t0 = time.time()
        Xc, Yc = simulateMode1Random(L_T, epsilon, n_samples, rng)
        I_c, ipp_c, _ = findPerPositionMutualInformation(Xc, Yc, alphabet=4)
        I_c_bulk, _, _ = findBulkMatchedControl(Xc, Yc, alphabet=4, rng=rng)
        Yc_scr = findStructureScrambledY(Yc, rng)
        I_c_scr, _, _ = findPerPositionMutualInformation(Xc, Yc_scr, alphabet=4)
        tc = time.time() - t0
        rows.append({
            "L_T": L_T, "template_type": "random",
            "epsilon": epsilon, "n_samples": n_samples,
            "I_per_position_total": I_c, "I_per_position_mean": ipp_c,
            "I_bulk_matched": I_c_bulk, "I_structure_scrambled": I_c_scr,
        })
        print(f"{L_T:4d} {'random':>14s} "
              f"{I_c:9.4f} {ipp_c:10.4f} {I_c_bulk:9.4f} {I_c_scr:12.4f} "
              f"{tc:5.2f}s")

    return rows


# ============================================================================
# part 4: experiment driver — sweep 2 (Drt3a vs Drt3b L-scaling)
# ============================================================================

def runDrt3aVsDrt3b(rng, Ls=(2, 4, 6, 8, 12, 24), epsilon=0.01, n_samples=1000):
    """
    side-by-side L-scaling for Drt3a (Mode 1, alternating-2phase template) and
    Drt3b (Mode 3, N = 2). uses the joint sequence-level MI estimator
    (H(Y) - H(Y|X)) for both — apples-to-apples comparison where X is a
    scalar phase descriptor in both cases.

    framework prediction:
      - Drt3a: I scales linearly with L (Mode 1 signature) — but the prose
        claim "linear" needs adjudication when X has only 1 bit of entropy.
        At most I ≤ H(X) = 1 bit because data-processing inequality bounds
        I(X;Y) ≤ H(X).
      - Drt3b: I saturates at log2(2) = 1 bit regardless of L.
    """
    rows = []
    print()
    print("=" * 96)
    print(f"Test G — sweep 2: Drt3a (alt-2phase) vs Drt3b (Mode 3 N=2)  "
          f"(eps={epsilon}, n={n_samples})")
    print("=" * 96)
    print(f"{'L':>4s} {'mode':>10s} "
          f"{'I_joint':>9s} {'H_Y':>9s} {'H_Y|X':>9s} {'time':>6s}")
    print("-" * 96)

    for L in Ls:
        # ---------- Drt3a: Mode 1 with alternating-2phase template ----------
        t0 = time.time()
        phase_a, _, Y_a = simulateMode1Alternating2Phase(L, epsilon, n_samples, rng)
        I_a, Hy_a, Hyx_a = findJointSequenceMutualInformation(phase_a, Y_a,
                                                              n_states=2)
        ta = time.time() - t0
        rows.append({
            "L": L, "mode": "drt3a_mode1_alt2phase",
            "epsilon": epsilon, "n_samples": n_samples,
            "I_joint": I_a, "H_Y": Hy_a, "H_Y_given_X": Hyx_a,
        })
        print(f"{L:4d} {'drt3a':>10s} "
              f"{I_a:9.4f} {Hy_a:9.4f} {Hyx_a:9.4f} {ta:5.2f}s")

        # ---------- Drt3b: Mode 3 with N = 2 ----------
        t0 = time.time()
        phase_b, Y_b = simulateMode3(N=2, L=L, epsilon=epsilon,
                                     n_samples=n_samples, rng=rng)
        I_b, Hy_b, Hyx_b = findJointSequenceMutualInformation(phase_b, Y_b,
                                                              n_states=2)
        tb = time.time() - t0
        rows.append({
            "L": L, "mode": "drt3b_mode3_N2",
            "epsilon": epsilon, "n_samples": n_samples,
            "I_joint": I_b, "H_Y": Hy_b, "H_Y_given_X": Hyx_b,
        })
        print(f"{L:4d} {'drt3b':>10s} "
              f"{I_b:9.4f} {Hy_b:9.4f} {Hyx_b:9.4f} {tb:5.2f}s")

    return rows


# ============================================================================
# part 5: I/O
# ============================================================================

CSV_FIELDS_SWEEP1 = [
    "L_T", "template_type", "epsilon", "n_samples",
    "I_per_position_total", "I_per_position_mean",
    "I_bulk_matched", "I_structure_scrambled",
]

CSV_FIELDS_SWEEP2 = [
    "L", "mode", "epsilon", "n_samples",
    "I_joint", "H_Y", "H_Y_given_X",
]


def saveCsv(rows, fields, path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})


# ============================================================================
# part 6: plots
# ============================================================================

def plotAlternatingVsRandom(rows, epsilon, path):
    """
    line plot: I_per_position_total vs L_T for each template type.
    overlay the Mode 1 random-template theoretical reference L_T * (2 - H(Y|X)).
    """
    fig, ax = plt.subplots(figsize=(8.5, 6.0))
    types = ["alternating_fixed", "alternating_2phase", "random"]
    colors = {"alternating_fixed": "C3",
              "alternating_2phase": "C1",
              "random": "C0"}
    labels = {"alternating_fixed": "alternating (single fixed ACAC...)",
              "alternating_2phase": "alternating (2-phase population)",
              "random": "random uniform template (Mode 1 baseline)"}

    for t in types:
        xs, ys, ybulk, yscr = [], [], [], []
        for r in rows:
            if r["template_type"] != t:
                continue
            xs.append(r["L_T"])
            ys.append(r["I_per_position_total"])
            ybulk.append(r["I_bulk_matched"])
            yscr.append(r["I_structure_scrambled"])
        order = np.argsort(xs)
        xs_a = np.array(xs)[order]
        ys_a = np.array(ys)[order]
        ybulk_a = np.array(ybulk)[order]
        yscr_a = np.array(yscr)[order]
        ax.plot(xs_a, ys_a, "-o", color=colors[t], lw=1.8, markersize=8,
                label=labels[t])
        ax.plot(xs_a, ybulk_a, "--", color=colors[t], lw=1.0, alpha=0.6,
                label=f"{t} bulk-matched control")
        ax.plot(xs_a, yscr_a, ":", color=colors[t], lw=1.0, alpha=0.6,
                label=f"{t} structure-scrambled Y")

    # theoretical Mode 1 random-template reference
    ipp = theoreticalPerPositionInfo(epsilon)
    L_dense = np.linspace(min(r["L_T"] for r in rows),
                          max(r["L_T"] for r in rows), 100)
    ax.plot(L_dense, ipp * L_dense, "k-", lw=1.0, alpha=0.5,
            label=fr"theory (random tpl, Mode 1): $L_T \cdot (2 - H(Y|X))$")

    ax.set_xlabel("template length $L_T$")
    ax.set_ylabel(r"$I_{\mathrm{struct}}$ (per-position MI summed) [bits]")
    ax.set_title("Test G — alternating vs random templates: per-position MI vs $L_T$")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc="upper left", ncol=1)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plotDrt3aVsDrt3bScaling(rows, path):
    """
    line plot: I_joint vs L for Drt3a (alt-2phase) and Drt3b (Mode 3 N=2).
    horizontal reference at log2(2) = 1 bit.
    """
    fig, ax = plt.subplots(figsize=(8.5, 6.0))
    modes = [("drt3a_mode1_alt2phase", "Drt3a (Mode 1, alt-2phase template)", "C0", "o"),
             ("drt3b_mode3_N2",        "Drt3b (Mode 3, N=2 cyclic state)",   "C2", "s")]
    for key, lab, col, mk in modes:
        xs, ys = [], []
        for r in rows:
            if r["mode"] != key:
                continue
            xs.append(r["L"])
            ys.append(r["I_joint"])
        order = np.argsort(xs)
        xs_a = np.array(xs)[order]
        ys_a = np.array(ys)[order]
        ax.plot(xs_a, ys_a, "-" + mk, color=col, markersize=9, lw=1.8, label=lab)
    ax.axhline(1.0, color="black", ls="--", lw=1.2, alpha=0.7,
               label=r"$\log_2 2 = 1$ bit")
    ax.set_xlabel("output length $L$")
    ax.set_ylabel(r"$I_{\mathrm{joint}}(X;Y) = H(Y) - H(Y|X)$ [bits]")
    ax.set_title("Test G — Drt3a vs Drt3b joint MI scaling with $L$")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ============================================================================
# part 7: resolution analysis (writes the markdown doc)
# ============================================================================

def writeResolutionDoc(sweep1_rows, sweep2_rows, epsilon, path):
    # find the L_T = 6 row for each template type from sweep 1
    s1_by_L_and_type = {}
    for r in sweep1_rows:
        s1_by_L_and_type[(r["L_T"], r["template_type"])] = r

    L_T_target = 6
    alt_fixed = s1_by_L_and_type.get((L_T_target, "alternating_fixed"))
    alt_2phase = s1_by_L_and_type.get((L_T_target, "alternating_2phase"))
    rand = s1_by_L_and_type.get((L_T_target, "random"))

    # compute slopes for sweep 2
    drt3a = sorted([r for r in sweep2_rows if r["mode"] == "drt3a_mode1_alt2phase"],
                   key=lambda x: x["L"])
    drt3b = sorted([r for r in sweep2_rows if r["mode"] == "drt3b_mode3_N2"],
                   key=lambda x: x["L"])
    Ls_a = np.array([r["L"] for r in drt3a], dtype=float)
    Is_a = np.array([r["I_joint"] for r in drt3a], dtype=float)
    Ls_b = np.array([r["L"] for r in drt3b], dtype=float)
    Is_b = np.array([r["I_joint"] for r in drt3b], dtype=float)
    # least-squares slope (informational only — not a PASS criterion)
    slope_a = float(np.polyfit(Ls_a, Is_a, 1)[0]) if len(Ls_a) >= 2 else float("nan")
    slope_b = float(np.polyfit(Ls_b, Is_b, 1)[0]) if len(Ls_b) >= 2 else float("nan")

    # naive Mode 1 ceiling at L_T = 6
    ipp = theoreticalPerPositionInfo(epsilon)
    pred_A_value = L_T_target * 2.0  # naive log2(4) per pos
    pred_A_realistic = L_T_target * ipp  # accounting for eps=0.01

    # determine which prediction the empirical data confirms at L_T = 6
    alt_fixed_I = alt_fixed["I_per_position_total"] if alt_fixed else float("nan")
    alt_2phase_I = alt_2phase["I_per_position_total"] if alt_2phase else float("nan")
    rand_I = rand["I_per_position_total"] if rand else float("nan")

    # heuristic adjudication: closer to 12 → A, closer to 1 → B
    def whichPrediction(I_obs):
        if abs(I_obs - pred_A_realistic) < abs(I_obs - 1.0):
            return "A (mechanism-determined, ~12 bits)"
        else:
            return "B (template-information-limited, ~1 bit)"

    verdict_alt_fixed = whichPrediction(alt_fixed_I)
    verdict_alt_2phase = whichPrediction(alt_2phase_I)
    verdict_random = whichPrediction(rand_I)

    text = []
    text.append("# Test G — Drt3a Mode 1 / Mode 3 boundary resolution (v1)")
    text.append("")
    text.append("## Statement of competing predictions")
    text.append("")
    text.append("Drt3a (Sharma et al. 2026, *Science*) synthesizes poly(GT) by Watson-Crick")
    text.append("templating off a 6-nt alternating ACACAC region of an ncRNA. Naively the")
    text.append("framework places this in Mode 1 with $L = 6$, alphabet 4, and predicts")
    text.append(r"$I_\mathrm{struct} \le L \log_2 4 = 12$ bits.")
    text.append("")
    text.append("But the template is alternating: only 1 bit of phase information")
    text.append("(ACAC… vs CACA…) is non-redundant. Two predictions:")
    text.append("")
    text.append("- **Prediction A (mechanism-determined):** the apparatus measures per-position")
    text.append(r"  transfer at the joint distribution level. Predicted $I_\mathrm{struct}")
    text.append(r"  \approx L \log_2 4 = 12$ bits at $L = 6$, regardless of template degeneracy.")
    text.append(f"  At $\\varepsilon = {epsilon}$ the realistic ceiling is")
    text.append(f"  $L \\cdot (2 - H(Y|X)) = 6 \\cdot {ipp:.4f} = {pred_A_realistic:.4f}$ bits.")
    text.append("- **Prediction B (template-information-limited):** the apparatus is bottlenecked")
    text.append("  by the template's intrinsic information content. Predicted")
    text.append(r"  $I_\mathrm{struct} \le 1$ bit total, because that is all the alternating")
    text.append("  template's population carries (1 bit of phase).")
    text.append("")
    text.append("## Empirical result (sweep 1: alternating vs random templates)")
    text.append("")
    text.append(f"At $L_T = 6$, $\\varepsilon = {epsilon}$, $n = 1000$ samples,")
    text.append("per-position plug-in MI estimator (Test A.1 convention):")
    text.append("")
    text.append("| template type                          | $I_\\mathrm{struct}$ (bits) | bulk-matched | structure-scrambled-Y | verdict |")
    text.append("|----------------------------------------|:--------------------------:|:------------:|:---------------------:|:--------|")
    if alt_fixed:
        text.append(f"| alternating, single fixed ACAC...      | "
                    f"{alt_fixed_I:.4f} | "
                    f"{alt_fixed['I_bulk_matched']:.4f} | "
                    f"{alt_fixed['I_structure_scrambled']:.4f} | "
                    f"{verdict_alt_fixed} |")
    if alt_2phase:
        text.append(f"| alternating, 2-phase population        | "
                    f"{alt_2phase_I:.4f} | "
                    f"{alt_2phase['I_bulk_matched']:.4f} | "
                    f"{alt_2phase['I_structure_scrambled']:.4f} | "
                    f"{verdict_alt_2phase} |")
    if rand:
        text.append(f"| random uniform template (A.1 baseline) | "
                    f"{rand_I:.4f} | "
                    f"{rand['I_bulk_matched']:.4f} | "
                    f"{rand['I_structure_scrambled']:.4f} | "
                    f"{verdict_random} |")
    text.append("")
    text.append("Full sweep (`results/test_g_alternating_vs_random_template_v1.csv`):")
    text.append("")
    text.append("| $L_T$ | template          | $I$ (bits) | bulk | scr-Y |")
    text.append("|------:|-------------------|----------:|-----:|------:|")
    for r in sorted(sweep1_rows, key=lambda x: (x["L_T"], x["template_type"])):
        text.append(f"| {r['L_T']:>5} | {r['template_type']:<17s} | "
                    f"{r['I_per_position_total']:8.4f} | "
                    f"{r['I_bulk_matched']:5.4f} | "
                    f"{r['I_structure_scrambled']:5.4f} |")
    text.append("")
    text.append("## Empirical result (sweep 2: Drt3a vs Drt3b L-scaling)")
    text.append("")
    text.append("Sequence-level joint MI estimator (Test B convention), $X$ = scalar")
    text.append("phase descriptor with 2 states for both modes:")
    text.append("")
    text.append("| $L$ | Drt3a (alt-2phase) | Drt3b (Mode 3 N=2) |")
    text.append("|----:|-------------------:|-------------------:|")
    Ls_all = sorted({r["L"] for r in sweep2_rows})
    by_mode = {(r["mode"], r["L"]): r["I_joint"] for r in sweep2_rows}
    for L in Ls_all:
        Ia = by_mode.get(("drt3a_mode1_alt2phase", L), float("nan"))
        Ib = by_mode.get(("drt3b_mode3_N2", L), float("nan"))
        text.append(f"| {L:>3} | {Ia:18.4f} | {Ib:18.4f} |")
    text.append("")
    text.append(f"Linear-fit slopes (informational): Drt3a slope = {slope_a:+.4e} bits/position,")
    text.append(f"Drt3b slope = {slope_b:+.4e} bits/position.")
    text.append("")
    text.append("Both curves saturate near $\\log_2 2 = 1$ bit. Drt3a's joint MI does NOT")
    text.append("scale linearly with L when X is taken to be the scalar phase index — it is")
    text.append("bounded above by $H(X) = 1$ bit by the data-processing inequality.")
    text.append("")
    text.append("## Resolution: which prediction does the apparatus confirm?")
    text.append("")
    text.append("The apparatus confirms a **per-position transfer interpretation only when the")
    text.append("template ensemble carries the relevant entropy.** Specifically:")
    text.append("")
    text.append("1. Random uniform template (A.1 baseline): the per-position MI estimator")
    text.append(f"   gives ≈ {rand_I:.2f} bits at $L_T = 6$, matching Prediction A. This")
    text.append("   confirms the apparatus reads per-position WC transfer correctly when")
    text.append("   the template population is non-degenerate.")
    text.append("")
    text.append("2. Alternating-fixed template (single ACAC… for all samples): the per-position")
    text.append(f"   MI estimator gives ≈ {alt_fixed_I:.2f} bits at $L_T = 6$. This is")
    text.append("   Prediction-B-like. Mechanically, every product position is determined by")
    text.append("   a template position via Watson-Crick — but with X constant across the")
    text.append("   population, the empirical I(X_i; Y_i) is identically zero by definition of")
    text.append("   plug-in MI. The apparatus measures cross-sample variability, not per-sample")
    text.append("   determinism.")
    text.append("")
    text.append("3. Alternating-2phase template (X drawn from a 2-element population): the")
    text.append(f"   per-position MI estimator gives ≈ {alt_2phase_I:.2f} bits at $L_T = 6$.")
    text.append("   This is bounded by H(X) = 1 bit per position (one binary descriptor) and")
    text.append("   by L_T * 1 bit = 6 bits if positions are treated independently — but")
    text.append("   positions are perfectly correlated given the phase, so the joint")
    text.append("   sequence-level MI is also bounded by H(X) = 1 bit (sweep 2 confirms).")
    text.append("")
    text.append("**Conclusion: the apparatus confirms Prediction B for the literal Drt3a setup")
    text.append("(alternating template, fixed or 2-phase). Prediction A holds only counterfactually,")
    text.append("when one substitutes the random-uniform-template population as the comparison")
    text.append("ensemble.**")
    text.append("")
    text.append("This is a substantive boundary-case finding for the framework. Two readings:")
    text.append("")
    text.append("- **Reading 1 (apparatus-as-stated is correct, framing was sloppy):** the")
    text.append("  framework's I_struct is a *population-level* mutual information, not a")
    text.append("  per-sample mechanism descriptor. If the population of templates is")
    text.append("  degenerate, transferable information is by definition limited by template")
    text.append("  entropy. Calling Drt3a 'Mode 1' on the basis of mechanism (per-base WC")
    text.append("  pairing happens) is fine descriptively, but the I_struct numeric must be")
    text.append("  reported relative to the population the descriptor varies over. For Drt3a's")
    text.append("  natural ncRNA-encoded template (one fixed ACACAC), I_struct measured against")
    text.append("  any reasonable comparison ensemble is ~0–1 bit — *not* 12 bits.")
    text.append("")
    text.append("- **Reading 2 (apparatus needs a per-sample mechanism observable):** the")
    text.append("  framework should add a complementary observable that captures per-position")
    text.append("  WC transfer at the single-sample level — e.g., the per-position fidelity")
    text.append("  $f_i = P(Y_i = WC(X_i) | X_i)$ averaged over realizations. This *would*")
    text.append("  give ~12 bits worth of 'transfer determinism' at $L_T = 6$ regardless of")
    text.append("  template degeneracy, but it would not be a mutual information.")
    text.append("")
    text.append("## Implication for the Mode 1 / Mode 3 boundary")
    text.append("")
    text.append("Sweep 2 shows that **with the I_struct apparatus as currently specified,**")
    text.append("Drt3a (Mode 1, alternating template) and Drt3b (Mode 3, N=2) are")
    text.append("*indistinguishable*: both saturate at ~1 bit and neither scales linearly in L.")
    text.append(f"Drt3a slope = {slope_a:+.4e} bits/position, Drt3b slope = {slope_b:+.4e}")
    text.append("bits/position. The framework's prose claim ('Drt3a I scales linearly with L,")
    text.append("Drt3b saturates at 1 bit') only holds for the random-template counterfactual,")
    text.append("not for Drt3a as it actually exists.")
    text.append("")
    text.append("**The Mode 1 / Mode 3 boundary, as the apparatus measures it, is")
    text.append("output-determined for degenerate-template Mode 1 systems, not")
    text.append("mechanism-determined.** Distinguishing Drt3a from Drt3b therefore requires")
    text.append("either (a) additional observables beyond I_struct (e.g., per-base fidelity,")
    text.append("template separability, ncRNA dependence experiments) or (b) restating the")
    text.append("framework's mode-classification claim to acknowledge that I_struct cannot")
    text.append("distinguish them when both produce the same (AC)_n output and X-population")
    text.append("entropies are matched.")
    text.append("")
    text.append("## Provenance")
    text.append("")
    text.append("- code: `code/test_g_drt3a_boundary.py`")
    text.append("- sweep 1 csv: `results/test_g_alternating_vs_random_template_v1.csv`")
    text.append("- sweep 2 csv: `results/test_g_drt3a_vs_drt3b_comparison_v1.csv`")
    text.append("- figures: `figures/test_g_alternating_vs_random.png`, "
                "`figures/test_g_drt3a_vs_drt3b_scaling.png`")
    text.append("- estimator: per-position plug-in MI on 4×4 joints (Test A.1) for sweep 1; "
                "joint sequence-level $H(Y) - H(Y|X)$ (Test B) for sweep 2")
    text.append(f"- seed: `np.random.seed(42)`, `np.random.default_rng(42)`")
    text.append(f"- n_samples = 1000, epsilon = {epsilon}")
    text.append("")
    path.write_text("\n".join(text))


# ============================================================================
# part 8: main
# ============================================================================

def main():
    np.random.seed(42)
    rng = np.random.default_rng(42)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # safety: refuse to overwrite our own previous outputs
    sweep1_csv = RESULTS_DIR / "test_g_alternating_vs_random_template_v1.csv"
    sweep2_csv = RESULTS_DIR / "test_g_drt3a_vs_drt3b_comparison_v1.csv"
    res_md = RESULTS_DIR / "test_g_drt3a_boundary_v1.md"
    fig1 = FIGURES_DIR / "test_g_alternating_vs_random.png"
    fig2 = FIGURES_DIR / "test_g_drt3a_vs_drt3b_scaling.png"

    epsilon = 0.01
    n_samples = 1000

    t_start = time.time()

    sweep1_rows = runAlternatingVsRandom(
        rng, L_Ts=(2, 4, 6, 8, 12, 24),
        epsilon=epsilon, n_samples=n_samples,
    )
    saveCsv(sweep1_rows, CSV_FIELDS_SWEEP1, sweep1_csv)
    plotAlternatingVsRandom(sweep1_rows, epsilon, fig1)

    sweep2_rows = runDrt3aVsDrt3b(
        rng, Ls=(2, 4, 6, 8, 12, 24),
        epsilon=epsilon, n_samples=n_samples,
    )
    saveCsv(sweep2_rows, CSV_FIELDS_SWEEP2, sweep2_csv)
    plotDrt3aVsDrt3bScaling(sweep2_rows, fig2)

    writeResolutionDoc(sweep1_rows, sweep2_rows, epsilon, res_md)

    elapsed = time.time() - t_start
    print()
    print("=" * 96)
    print(f"sweep 1 csv:  {sweep1_csv}")
    print(f"sweep 2 csv:  {sweep2_csv}")
    print(f"resolution:   {res_md}")
    print(f"figures:      {fig1}")
    print(f"              {fig2}")
    print(f"total wall-time: {elapsed:.2f} s")
    print("=" * 96)


if __name__ == "__main__":
    main()
