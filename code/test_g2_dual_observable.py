"""
test G2 — apparatus repair via dual observables (population MI + channel-as-X MI)
================================================================================

Templating Substrates Framework, Test G2 (apparatus repair)
-----------------------------------------------------------
Test G v1 found that the I_struct apparatus, computed with X = template
realization (= "I_struct^pop", population MI), gives ~0 bits for Drt3a's
actual fixed-ACACAC ncRNA template and ~1 bit for the 2-phase alternating
counterfactual. This breaks the v2 draft's central distinguishing claim
between Mode 1 (Drt3a) and Mode 3 (Drt3b N=2): both saturate at ~1 bit.

G v1's resolution document proposed two readings ("apparatus-as-stated"
vs "needs per-sample observable"). G2 reframes them as COMPLEMENTARY
observables, not competing repairs:

  Observable 1: I_struct^pop  (X = template realization)
                = H(Y) - H(Y|X = x), averaged over x in the template population.
                This is what Test G v1 measured. Limited by H(template ensemble).

  Observable 2: I_struct^chan (X = the channel identity itself)
                = sum_c P(C=c) * KL( P(Y|C=c) || P(Y) ),
                where C ranges over a fixed channel ensemble
                {Drt3a-WC, Drt3b-Mode3-N=2, Drt3b-Mode3-N=3, Drt3b-Mode3-N=5,
                 AbiK-uniform}. The per-channel contribution
                KL( P(Y|C=c_focal) || P(Y) ) is the focal channel's MI signature
                against the comparison ensemble.

The product (I_struct^pop, I_struct^chan, per-base fidelity) is the refined
apparatus signature. Mode classification uses the joint signature, not
either observable alone. This file computes both observables on the v1
sweep harness and reports them jointly.

This file is intentionally self-contained: simulators and estimators are
duplicated from test_g_drt3a_boundary.py (which duplicates from
test_a1_mode1_scaling.py and test_b_mode3_capacity.py). per repo CLAUDE.md
rule 4 (no shared utils module).
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

# fixed channel ensemble for I_chan computation (5 channels)
# each entry is (name, mode_kind, params)
#   mode_kind ∈ {"mode1_alt2phase", "mode3", "uniform_random"}
#   params: dict with kind-specific keys
CHANNEL_ENSEMBLE = [
    ("Drt3a-WC",       "mode1_alt2phase", {}),
    ("Drt3b-Mode3-N2", "mode3",           {"N": 2}),
    ("Drt3b-Mode3-N3", "mode3",           {"N": 3}),
    ("Drt3b-Mode3-N5", "mode3",           {"N": 5}),
    ("AbiK-uniform",   "uniform_random",  {}),
]


# ============================================================================
# part 1: simulators (duplicated from test_g_drt3a_boundary.py)
# ============================================================================

# -- Mode 1 (Drt3a) ----------------------------------------------------------

def simulateMode1Random(L, epsilon, n_samples, rng):
    """
    Mode 1 with uniformly random templates per sample (Test A.1 baseline).

    each of n_samples gets its own independently drawn template X of length L,
    and produces Y = WC(X) with per-position fidelity 1 - epsilon. returns
    (X, Y) both shape (n_samples, L), dtype int8.
    """
    X = rng.integers(0, 4, size=(n_samples, L), dtype=np.int8)
    Y_correct = WC_PAIR[X]
    is_error = rng.random(size=(n_samples, L)) < epsilon
    offset = rng.integers(0, 3, size=(n_samples, L), dtype=np.int8)
    Y_wrong = np.where(offset >= Y_correct, offset + 1, offset).astype(np.int8)
    Y = np.where(is_error, Y_wrong, Y_correct).astype(np.int8)
    return X, Y


def simulateMode1AlternatingFixed(L, epsilon, n_samples, rng):
    """
    Mode 1 with a single fixed alternating ACAC...AC template for all samples.
    X is identical across samples; per-sample mechanism is identical to Mode 1.
    """
    template = np.tile(np.array([0, 1], dtype=np.int8), (L + 1) // 2)[:L]
    X = np.broadcast_to(template, (n_samples, L)).copy()
    Y_correct = WC_PAIR[X]
    is_error = rng.random(size=(n_samples, L)) < epsilon
    offset = rng.integers(0, 3, size=(n_samples, L), dtype=np.int8)
    Y_wrong = np.where(offset >= Y_correct, offset + 1, offset).astype(np.int8)
    Y = np.where(is_error, Y_wrong, Y_correct).astype(np.int8)
    return X, Y


def simulateMode1Alternating2Phase(L, epsilon, n_samples, rng):
    """
    Mode 1 with alternating template drawn from a 2-element population.
    each sample picks phase ∈ {0, 1} uniformly:
      phase 0 → X = [A, C, A, C, ...]
      phase 1 → X = [C, A, C, A, ...]
    returns (phase, X, Y).
    """
    phase = rng.integers(0, 2, size=n_samples, dtype=np.int64)
    pos = np.arange(L, dtype=np.int64)[None, :]
    X = ((phase[:, None] + pos) % 2).astype(np.int8)
    Y_correct = WC_PAIR[X]
    is_error = rng.random(size=(n_samples, L)) < epsilon
    offset = rng.integers(0, 3, size=(n_samples, L), dtype=np.int8)
    Y_wrong = np.where(offset >= Y_correct, offset + 1, offset).astype(np.int8)
    Y = np.where(is_error, Y_wrong, Y_correct).astype(np.int8)
    return phase, X, Y


# -- Mode 3 (Drt3b) ----------------------------------------------------------

def simulateMode3(N, L, epsilon, n_samples, rng):
    """
    Mode 3 cyclic active-site channel.
    X = phi_0 ∈ {0, ..., N-1} uniform; intended Y_k = (phi_0 + k) mod N.
    P(Y_k = correct | X) = 1 - eps; P(Y_k = each wrong | X) = eps / (N-1).
    NOTE: Y values live in {0, ..., N-1}, NOT {0..3}. for 4-letter Y-space
    mapping we use simulateMode3As4Letter below.
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


def simulateMode3As4Letter(N, L, epsilon, n_samples, rng):
    """
    Mode 3 cyclic active-site channel mapped into the 4-letter {0,1,2,3} alphabet
    so Y rows live in the same space as Mode 1 / AbiK outputs and pooled
    distributional comparisons (I_chan) are well-defined.

    mapping: phase index k ∈ {0,...,N-1} → nucleotide via the lookup
        N=2: {0:A=0, 1:C=1}                     (poly(AC) — Drt3b output)
        N=3: {0:A=0, 1:C=1, 2:G=2}              (3-cycle output)
        N=5: {0:A=0, 1:C=1, 2:G=2, 3:T=3, 4:A=0}(5-cycle on 4-letter alphabet,
                                                 with one repeat — best
                                                 squeeze of N=5 into a 4-letter
                                                 product space)
    misincorporations occur in the {0..N-1} phase space and are then mapped
    to the 4-letter alphabet via the same lookup. this preserves the framework's
    Mode 3 mechanism (cyclic phase, fidelity 1-eps per position) while putting
    Y in the same alphabet as Mode 1 / AbiK so I_chan is computable.
    """
    # canonical lookups to 4-letter alphabet
    lookups = {
        2: np.array([0, 1], dtype=np.int8),
        3: np.array([0, 1, 2], dtype=np.int8),
        5: np.array([0, 1, 2, 3, 0], dtype=np.int8),
    }
    if N not in lookups:
        raise ValueError(f"no 4-letter lookup defined for N={N}")
    lookup = lookups[N]
    X_phase, Y_phase = simulateMode3(N=N, L=L, epsilon=epsilon,
                                     n_samples=n_samples, rng=rng)
    Y_4letter = lookup[Y_phase].astype(np.int8)
    return X_phase, Y_4letter


def simulateUniformRandom(L, epsilon, n_samples, rng):
    """
    AbiK-style channel: produces uniform random output independent of any
    template. epsilon is ignored (the channel is intrinsically maximally
    noisy). returns (X, Y) with X = None placeholder (no template) and
    Y shape (n_samples, L) dtype int8. for I_pop we treat X = "no template"
    (single class) so I_pop = 0 by construction, matching the framework's
    "AbiK is not a templating system" verdict.
    """
    Y = rng.integers(0, 4, size=(n_samples, L), dtype=np.int8)
    return None, Y


# ============================================================================
# part 2: estimators
# ============================================================================

# -- per-position plug-in MI (Test A.1 convention) ---------------------------

def findPerPositionMutualInformation(X, Y, alphabet=4):
    """
    plug-in I(X_i; Y_i) per position from empirical (alphabet x alphabet) joints,
    summed across positions. for (n_samples, L) X and Y.

    used here as the Observable 1 "I_struct^pop" so v1 numbers are reproduced
    exactly (modulo n_samples=5000 vs v1's 1000).

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


def findBulkMatchedControl(X, Y, alphabet, rng):
    """bulk-matched scramble: row-permute X. preserves X composition,
    breaks per-position alignment with Y."""
    perm = rng.permutation(X.shape[0])
    X_scrambled = X[perm]
    return findPerPositionMutualInformation(X_scrambled, Y, alphabet=alphabet)


# -- per-base fidelity (the per-sample mechanism observable, for Step 3) ----

def findPerBaseFidelity(X, Y):
    """
    per-position fidelity f_i = P(Y_i = WC(X_i)) over samples, then averaged
    over positions. only well-defined when WC pairing applies (Mode 1, where
    X is a 4-letter template). returns mean f and per-position array.
    """
    Y_correct = WC_PAIR[X.astype(np.int8)]
    match = (Y == Y_correct).astype(np.float64)
    f_per_pos = match.mean(axis=0)
    return float(f_per_pos.mean()), f_per_pos


# -- distributional helpers for I_chan (Observable 2) -----------------------

def encodeRowsToInts(Y, alphabet=4):
    """
    encode each row of Y (shape (n, L), values in {0..alphabet-1}) as a
    single int64. used as a hash key for empirical row-distribution counts.
    """
    n, L = Y.shape
    Yc = Y.astype(np.int64)
    powers = (alphabet ** np.arange(L, dtype=np.int64))[::-1]  # MSB first
    return Yc @ powers


def empiricalRowDistribution(Y, alphabet=4):
    """
    plug-in empirical distribution over Y-rows. returns dict {int64 row-code : prob}.
    """
    codes = encodeRowsToInts(Y, alphabet=alphabet)
    n = len(codes)
    counts = {}
    for c in codes:
        c_i = int(c)
        counts[c_i] = counts.get(c_i, 0) + 1
    return {k: v / n for k, v in counts.items()}


def klDivergence(p, q):
    """
    KL(p || q) in bits, where p and q are dicts {key: prob}. keys missing
    in q are treated as 0 → KL = +inf if p has positive mass on them
    (we clip via EPS_CLIP to avoid inf and to handle finite-sample bias).
    returns finite KL in bits.
    """
    kl = 0.0
    for k, p_k in p.items():
        if p_k <= 0:
            continue
        q_k = q.get(k, EPS_CLIP)
        if q_k <= 0:
            q_k = EPS_CLIP
        kl += p_k * np.log2(p_k / q_k)
    return float(kl)


def jsDivergence(p, q):
    """
    Jensen-Shannon divergence in bits. symmetric, bounded above by 1 bit.
    """
    keys = set(p.keys()) | set(q.keys())
    m = {k: 0.5 * (p.get(k, 0.0) + q.get(k, 0.0)) for k in keys}
    return 0.5 * klDivergence(p, m) + 0.5 * klDivergence(q, m)


def mixDistributions(dists, weights=None):
    """
    convex combination of a list of {key: prob} distributions.
    weights default to uniform.
    """
    n = len(dists)
    if weights is None:
        weights = [1.0 / n] * n
    keys = set()
    for d in dists:
        keys |= d.keys()
    out = {}
    for k in keys:
        out[k] = sum(w * d.get(k, 0.0) for w, d in zip(weights, dists))
    return out


# ============================================================================
# part 3: channel-ensemble simulator (5-channel ensemble at given L)
# ============================================================================

def simulateChannelEnsemble(L, epsilon, n_samples_per_channel, rng):
    """
    generate n_samples_per_channel rows of Y from each of the 5 channels
    in CHANNEL_ENSEMBLE. all Y rows are in the 4-letter alphabet so they
    live in a common Y-space at length L. returns a list of
    (channel_name, Y_rows) tuples, in the order of CHANNEL_ENSEMBLE.
    """
    out = []
    for name, kind, params in CHANNEL_ENSEMBLE:
        if kind == "mode1_alt2phase":
            _, _, Y = simulateMode1Alternating2Phase(L=L, epsilon=epsilon,
                                                    n_samples=n_samples_per_channel,
                                                    rng=rng)
        elif kind == "mode3":
            _, Y = simulateMode3As4Letter(N=params["N"], L=L, epsilon=epsilon,
                                          n_samples=n_samples_per_channel,
                                          rng=rng)
        elif kind == "uniform_random":
            _, Y = simulateUniformRandom(L=L, epsilon=epsilon,
                                         n_samples=n_samples_per_channel,
                                         rng=rng)
        else:
            raise ValueError(f"unknown channel kind {kind}")
        out.append((name, Y))
    return out


def findIChan(focal_Y, ensemble, alphabet=4):
    """
    compute the channel-as-X MI contribution for the focal channel against
    the channel ensemble.

    formally: I(C; Y) = sum_c P(C=c) * KL( P(Y|C=c) || P(Y) )
    we report the focal channel's per-channel summand, normalized as if all
    channels had equal mass:
        I_chan_focal = KL( P(Y|focal) || P(Y_ensemble_mean) )
    plus the total ensemble I(C;Y) for context.

    parameters
    ----------
    focal_Y : (n, L) array, the rows from the focal channel.
    ensemble: list of (name, Y) for all channels (uniform mixture).
    alphabet: alphabet size.

    returns dict with keys
        I_chan_focal_kl : KL(focal || mixture) in bits.
        I_chan_total    : full I(C;Y) for the ensemble in bits.
        focal_in_ens    : True if focal channel matches one of the ensemble
                          channels by row-distribution (for sanity).
    """
    p_focal = empiricalRowDistribution(focal_Y, alphabet=alphabet)
    dists = [empiricalRowDistribution(Y, alphabet=alphabet) for _, Y in ensemble]
    n_ch = len(ensemble)
    p_mix = mixDistributions(dists, weights=[1.0 / n_ch] * n_ch)
    kl_focal = klDivergence(p_focal, p_mix)
    # full I(C;Y)
    I_total = 0.0
    for d in dists:
        I_total += (1.0 / n_ch) * klDivergence(d, p_mix)
    return {
        "I_chan_focal_kl": kl_focal,
        "I_chan_total":    I_total,
    }


# ============================================================================
# part 4: periodicity and marginals
# ============================================================================

def findPeriodicityPeak(Y, max_lag=12):
    """
    autocorrelation of Y rows averaged over rows and over the alphabet:
        rho(lag) = mean over (i, k) of (1 if Y[i,k] == Y[i,k+lag] else 0)
    return (peak_lag, peak_value) over lag ∈ {2, 3, ..., min(max_lag, L-1)}.
    """
    n, L = Y.shape
    if L < 3:
        return 0, 0.0
    best_lag, best_val = 0, 0.0
    for lag in range(2, min(max_lag + 1, L)):
        match = (Y[:, :L - lag] == Y[:, lag:]).astype(np.float64).mean()
        if match > best_val:
            best_val = match
            best_lag = lag
    return best_lag, float(best_val)


def findMarginalFreqs(Y, alphabet=4):
    """marginal frequencies of each alphabet symbol in Y."""
    flat = Y.ravel().astype(np.int64)
    counts = np.bincount(flat, minlength=alphabet)
    return (counts / counts.sum()).tolist()


# ============================================================================
# part 5: STEP 1 — full v1 sweep with both observables
# ============================================================================

def runStep1Sweep(rng, L_Ts=(2, 4, 6, 8, 12, 24), epsilon=0.01, n_samples=5000):
    """
    for each (L_T, template_type) cell from v1, compute both observables
    plus periodicity and marginal stats. template_types: alt_fixed,
    alt_2phase, random.
    """
    rows = []
    print()
    print("=" * 110)
    print(f"Step 1 — v1 sweep with both observables (eps={epsilon}, n={n_samples})")
    print("=" * 110)
    print(f"{'L_T':>4s} {'tpl_type':>16s} "
          f"{'I_pop':>9s} {'I_chan':>9s} {'sep_pop':>9s} {'sep_chan':>9s} "
          f"{'period_pk':>10s} {'mA':>6s} {'mC':>6s} {'mG':>6s} {'mT':>6s} "
          f"{'time':>6s}")
    print("-" * 110)

    for L_T in L_Ts:
        # generate the channel ensemble at this L_T (5 channels x n_samples each)
        ensemble = simulateChannelEnsemble(L=L_T, epsilon=epsilon,
                                           n_samples_per_channel=n_samples,
                                           rng=rng)

        for tpl_type, sim_fn in [
            ("alt_fixed",   simulateMode1AlternatingFixed),
            ("alt_2phase",  lambda L, e, n, r: simulateMode1Alternating2Phase(L, e, n, r)[1:]),
            ("random",      simulateMode1Random),
        ]:
            t0 = time.time()
            X, Y = sim_fn(L_T, epsilon, n_samples, rng)

            # observable 1: I_struct^pop (per-position plug-in MI summed over L)
            I_pop, _, _ = findPerPositionMutualInformation(X, Y, alphabet=4)
            # bulk-matched control for I_pop
            I_pop_bulk, _, _ = findBulkMatchedControl(X, Y, alphabet=4, rng=rng)

            # observable 2: I_struct^chan against the 5-channel ensemble
            chan_stats = findIChan(Y, ensemble, alphabet=4)
            I_chan = chan_stats["I_chan_focal_kl"]

            # bulk-matched separation for I_chan: use a "no-templating-row"
            # control by row-shuffling Y across positions per row, then re-doing
            # KL against the same ensemble. crude but consistent
            Y_scr = np.empty_like(Y)
            for r_idx in range(Y.shape[0]):
                Y_scr[r_idx] = rng.permutation(Y[r_idx])
            chan_stats_scr = findIChan(Y_scr, ensemble, alphabet=4)
            I_chan_bulk = chan_stats_scr["I_chan_focal_kl"]

            sep_pop  = I_pop  - I_pop_bulk
            sep_chan = I_chan - I_chan_bulk

            period_lag, period_val = findPeriodicityPeak(Y)
            margins = findMarginalFreqs(Y, alphabet=4)

            elapsed = time.time() - t0
            rows.append({
                "template_type": tpl_type,
                "L_T":           L_T,
                "epsilon":       epsilon,
                "n_samples":     n_samples,
                "I_pop":         I_pop,
                "I_chan":        I_chan,
                "I_pop_bulk":    I_pop_bulk,
                "I_chan_bulk":   I_chan_bulk,
                "sep_pop":       sep_pop,
                "sep_chan":      sep_chan,
                "period_lag":    period_lag,
                "period_peak":   period_val,
                "margin_A":      margins[0],
                "margin_C":      margins[1],
                "margin_G":      margins[2],
                "margin_T":      margins[3],
            })
            print(f"{L_T:4d} {tpl_type:>16s} "
                  f"{I_pop:9.4f} {I_chan:9.4f} {sep_pop:9.4f} {sep_chan:9.4f} "
                  f"{period_val:10.4f} "
                  f"{margins[0]:6.3f} {margins[1]:6.3f} "
                  f"{margins[2]:6.3f} {margins[3]:6.3f} "
                  f"{elapsed:5.2f}s")
    return rows


# ============================================================================
# part 6: STEP 2 — 5×5 cross-mode separation matrix
# ============================================================================

def runStep2SeparationMatrix(rng, L=6, epsilon=0.01, n_samples=5000):
    """
    pairwise cross-mode separation for the 5 channels at fixed L.
    two parallel matrices:
      - JS divergence (symmetric, bounded by 1 bit)
      - I_chan contribution: KL(P(Y|c_i) || P(Y|c_j)) (asymmetric)
    """
    print()
    print("=" * 96)
    print(f"Step 2 — 5×5 cross-mode separation matrix (L={L}, eps={epsilon}, n={n_samples})")
    print("=" * 96)

    ensemble = simulateChannelEnsemble(L=L, epsilon=epsilon,
                                       n_samples_per_channel=n_samples,
                                       rng=rng)
    names = [n for n, _ in ensemble]
    dists = [empiricalRowDistribution(Y, alphabet=4) for _, Y in ensemble]
    n_ch = len(names)

    js_matrix = np.zeros((n_ch, n_ch))
    kl_matrix = np.zeros((n_ch, n_ch))
    for i in range(n_ch):
        for j in range(n_ch):
            if i == j:
                continue
            js_matrix[i, j] = jsDivergence(dists[i], dists[j])
            kl_matrix[i, j] = klDivergence(dists[i], dists[j])

    # tabular print
    print(f"{'JS divergence (bits)':<40s}")
    print(f"{'':<22s}" + "".join(f"{n:>16s}" for n in names))
    for i, n_i in enumerate(names):
        print(f"{n_i:<22s}" + "".join(f"{js_matrix[i, j]:16.4f}" for j in range(n_ch)))
    print()
    print(f"{'KL divergence i→j (bits)':<40s}")
    print(f"{'':<22s}" + "".join(f"{n:>16s}" for n in names))
    for i, n_i in enumerate(names):
        print(f"{n_i:<22s}" + "".join(f"{kl_matrix[i, j]:16.4f}" for j in range(n_ch)))

    # build flat row representation for csv
    rows = []
    for i, n_i in enumerate(names):
        for j, n_j in enumerate(names):
            rows.append({
                "channel_i":  n_i,
                "channel_j":  n_j,
                "L":          L,
                "epsilon":    epsilon,
                "n_samples":  n_samples,
                "js_divergence": js_matrix[i, j],
                "kl_divergence": kl_matrix[i, j],
            })
    return rows, names, js_matrix, kl_matrix


# ============================================================================
# part 7: STEP 3 — Drt3a under three template ensembles (E1/E2/E3)
# ============================================================================

def runStep3Drt3aThreeEnsembles(rng, L=6, epsilon=0.01, n_samples=5000):
    """
    full apparatus signature for Drt3a under:
      E1: single fixed ACACAC template
      E2: 2-phase {ACACAC, CACACA}
      E3: random uniform 6-nt template
    each reports (I_pop, I_chan, per-base fidelity, JS vs Drt3b).
    """
    print()
    print("=" * 110)
    print(f"Step 3 — Drt3a under E1/E2/E3 (L={L}, eps={epsilon}, n={n_samples})")
    print("=" * 110)
    print(f"{'ensemble':>10s} {'I_pop':>9s} {'I_chan':>9s} "
          f"{'f_per_base':>11s} {'JS_vs_Drt3b':>13s}")
    print("-" * 110)

    # ensemble for I_chan
    ensemble = simulateChannelEnsemble(L=L, epsilon=epsilon,
                                       n_samples_per_channel=n_samples,
                                       rng=rng)
    # Drt3b (Mode 3, N=2) reference distribution for JS
    drt3b_dist = empiricalRowDistribution(
        [Y for n, Y in ensemble if n == "Drt3b-Mode3-N2"][0], alphabet=4
    )

    rows = []

    # ---- E1: single fixed ACACAC ----
    X1, Y1 = simulateMode1AlternatingFixed(L, epsilon, n_samples, rng)
    I_pop_1, _, _ = findPerPositionMutualInformation(X1, Y1, alphabet=4)
    I_chan_1 = findIChan(Y1, ensemble, alphabet=4)["I_chan_focal_kl"]
    f_1, _ = findPerBaseFidelity(X1, Y1)
    js_1 = jsDivergence(empiricalRowDistribution(Y1, alphabet=4), drt3b_dist)
    rows.append({
        "ensemble":   "E1_fixed_ACACAC",
        "L":          L,
        "epsilon":    epsilon,
        "n_samples":  n_samples,
        "I_pop":      I_pop_1,
        "I_chan":     I_chan_1,
        "f_per_base": f_1,
        "JS_vs_Drt3b": js_1,
    })
    print(f"{'E1':>10s} {I_pop_1:9.4f} {I_chan_1:9.4f} {f_1:11.4f} {js_1:13.4f}")

    # ---- E2: 2-phase ----
    _, X2, Y2 = simulateMode1Alternating2Phase(L, epsilon, n_samples, rng)
    I_pop_2, _, _ = findPerPositionMutualInformation(X2, Y2, alphabet=4)
    I_chan_2 = findIChan(Y2, ensemble, alphabet=4)["I_chan_focal_kl"]
    f_2, _ = findPerBaseFidelity(X2, Y2)
    js_2 = jsDivergence(empiricalRowDistribution(Y2, alphabet=4), drt3b_dist)
    rows.append({
        "ensemble":   "E2_2phase_alt",
        "L":          L,
        "epsilon":    epsilon,
        "n_samples":  n_samples,
        "I_pop":      I_pop_2,
        "I_chan":     I_chan_2,
        "f_per_base": f_2,
        "JS_vs_Drt3b": js_2,
    })
    print(f"{'E2':>10s} {I_pop_2:9.4f} {I_chan_2:9.4f} {f_2:11.4f} {js_2:13.4f}")

    # ---- E3: random uniform 6-nt ----
    X3, Y3 = simulateMode1Random(L, epsilon, n_samples, rng)
    I_pop_3, _, _ = findPerPositionMutualInformation(X3, Y3, alphabet=4)
    I_chan_3 = findIChan(Y3, ensemble, alphabet=4)["I_chan_focal_kl"]
    f_3, _ = findPerBaseFidelity(X3, Y3)
    js_3 = jsDivergence(empiricalRowDistribution(Y3, alphabet=4), drt3b_dist)
    rows.append({
        "ensemble":   "E3_random_uniform",
        "L":          L,
        "epsilon":    epsilon,
        "n_samples":  n_samples,
        "I_pop":      I_pop_3,
        "I_chan":     I_chan_3,
        "f_per_base": f_3,
        "JS_vs_Drt3b": js_3,
    })
    print(f"{'E3':>10s} {I_pop_3:9.4f} {I_chan_3:9.4f} {f_3:11.4f} {js_3:13.4f}")

    return rows


# ============================================================================
# part 8: STEP 3b — Drt3b at N=2, N=3 and AbiK signature for the joint table
# ============================================================================

def runStep3bOtherCases(rng, L=6, epsilon=0.01, n_samples=5000):
    """
    compute the joint (I_pop, I_chan, f, JS_vs_Drt3b) signature for
    Drt3b-N=2, Drt3b-N=3, AbiK so we have 6 total test cases.
    NOTE: per-base fidelity for Mode 3 / AbiK is not WC; we report
    P(Y_k = (X_phase + k) mod N mapped to 4-letter) for Mode 3
    (i.e., per-position mechanism fidelity, not Watson-Crick fidelity)
    and for AbiK we report 0.25 (chance).
    """
    print()
    print("=" * 110)
    print(f"Step 3b — Drt3b N=2, N=3, AbiK signatures (L={L}, eps={epsilon}, n={n_samples})")
    print("=" * 110)

    ensemble = simulateChannelEnsemble(L=L, epsilon=epsilon,
                                       n_samples_per_channel=n_samples,
                                       rng=rng)
    drt3b_dist = empiricalRowDistribution(
        [Y for n, Y in ensemble if n == "Drt3b-Mode3-N2"][0], alphabet=4
    )

    rows = []

    # Drt3b N=2: I_pop = H(Y) - H(Y|X=phase) at row level using joint estimator
    # we use the per-position MI for Y in 4-letter space with X = phase (size 2)
    # but per-position MI requires X and Y of same shape. for row-level joint MI
    # we use H(Y) - H(Y|X) on row codes.
    X_b2_phase, Y_b2 = simulateMode3As4Letter(N=2, L=L, epsilon=epsilon,
                                              n_samples=n_samples, rng=rng)
    # row-level joint MI for Mode 3 (X = scalar phase ∈ {0..N-1})
    I_pop_b2 = jointMIRowsScalarX(X_b2_phase, Y_b2, n_states=2)
    I_chan_b2 = findIChan(Y_b2, ensemble, alphabet=4)["I_chan_focal_kl"]
    # per-position mechanism fidelity for Mode 3:
    # check Y[k] equals lookup[(phase+k) mod N]
    lookup_2 = np.array([0, 1], dtype=np.int8)
    pos = np.arange(L)[None, :]
    Y_pred_b2 = lookup_2[(X_b2_phase[:, None] + pos) % 2]
    f_b2 = float((Y_b2 == Y_pred_b2).mean())
    js_b2 = jsDivergence(empiricalRowDistribution(Y_b2, alphabet=4), drt3b_dist)

    rows.append({
        "ensemble":   "Drt3b_N2",
        "L":          L,
        "epsilon":    epsilon,
        "n_samples":  n_samples,
        "I_pop":      I_pop_b2,
        "I_chan":     I_chan_b2,
        "f_per_base": f_b2,
        "JS_vs_Drt3b": js_b2,
    })
    print(f"{'Drt3b_N2':>20s} I_pop={I_pop_b2:.4f}  I_chan={I_chan_b2:.4f}  "
          f"f={f_b2:.4f}  JS_vs_Drt3b={js_b2:.4f}")

    # Drt3b N=3
    X_b3_phase, Y_b3 = simulateMode3As4Letter(N=3, L=L, epsilon=epsilon,
                                              n_samples=n_samples, rng=rng)
    I_pop_b3 = jointMIRowsScalarX(X_b3_phase, Y_b3, n_states=3)
    I_chan_b3 = findIChan(Y_b3, ensemble, alphabet=4)["I_chan_focal_kl"]
    lookup_3 = np.array([0, 1, 2], dtype=np.int8)
    Y_pred_b3 = lookup_3[(X_b3_phase[:, None] + pos) % 3]
    f_b3 = float((Y_b3 == Y_pred_b3).mean())
    js_b3 = jsDivergence(empiricalRowDistribution(Y_b3, alphabet=4), drt3b_dist)
    rows.append({
        "ensemble":   "Drt3b_N3",
        "L":          L,
        "epsilon":    epsilon,
        "n_samples":  n_samples,
        "I_pop":      I_pop_b3,
        "I_chan":     I_chan_b3,
        "f_per_base": f_b3,
        "JS_vs_Drt3b": js_b3,
    })
    print(f"{'Drt3b_N3':>20s} I_pop={I_pop_b3:.4f}  I_chan={I_chan_b3:.4f}  "
          f"f={f_b3:.4f}  JS_vs_Drt3b={js_b3:.4f}")

    # AbiK: uniform random; I_pop = 0 (no template), I_chan against ensemble,
    # per-base fidelity = 0.25 (chance)
    _, Y_abik = simulateUniformRandom(L=L, epsilon=epsilon,
                                      n_samples=n_samples, rng=rng)
    I_pop_abik = 0.0  # no template descriptor; X = constant
    I_chan_abik = findIChan(Y_abik, ensemble, alphabet=4)["I_chan_focal_kl"]
    f_abik = 0.25
    js_abik = jsDivergence(empiricalRowDistribution(Y_abik, alphabet=4), drt3b_dist)
    rows.append({
        "ensemble":   "AbiK_uniform",
        "L":          L,
        "epsilon":    epsilon,
        "n_samples":  n_samples,
        "I_pop":      I_pop_abik,
        "I_chan":     I_chan_abik,
        "f_per_base": f_abik,
        "JS_vs_Drt3b": js_abik,
    })
    print(f"{'AbiK_uniform':>20s} I_pop={I_pop_abik:.4f}  I_chan={I_chan_abik:.4f}  "
          f"f={f_abik:.4f}  JS_vs_Drt3b={js_abik:.4f}")

    return rows


def jointMIRowsScalarX(X_scalar, Y, n_states):
    """
    sequence-level I(X;Y) = H(Y) - H(Y|X) for X scalar in {0..n_states-1},
    Y shape (n, L). plug-in entropy on row codes.
    """
    n_samples, L = Y.shape
    codes = encodeRowsToInts(Y, alphabet=4)
    counts_y = {}
    for c in codes:
        c_i = int(c)
        counts_y[c_i] = counts_y.get(c_i, 0) + 1
    p_y = np.array(list(counts_y.values()), dtype=np.float64) / n_samples
    H_y = float(-(p_y * np.log2(p_y)).sum())
    # H(Y|X)
    H_y_given_x = 0.0
    X_arr = np.asarray(X_scalar, dtype=np.int64)
    for x in range(n_states):
        idx = np.flatnonzero(X_arr == x)
        if len(idx) == 0:
            continue
        sub_codes = codes[idx]
        sub_counts = {}
        for c in sub_codes:
            c_i = int(c)
            sub_counts[c_i] = sub_counts.get(c_i, 0) + 1
        p_sub = np.array(list(sub_counts.values()), dtype=np.float64) / len(idx)
        H_sub = float(-(p_sub * np.log2(p_sub)).sum())
        H_y_given_x += (len(idx) / n_samples) * H_sub
    return H_y - H_y_given_x


# ============================================================================
# part 9: I/O
# ============================================================================

CSV_FIELDS_STEP1 = [
    "template_type", "L_T", "epsilon", "n_samples",
    "I_pop", "I_chan", "I_pop_bulk", "I_chan_bulk",
    "sep_pop", "sep_chan",
    "period_lag", "period_peak",
    "margin_A", "margin_C", "margin_G", "margin_T",
]

CSV_FIELDS_STEP2 = [
    "channel_i", "channel_j", "L", "epsilon", "n_samples",
    "js_divergence", "kl_divergence",
]

CSV_FIELDS_STEP3 = [
    "ensemble", "L", "epsilon", "n_samples",
    "I_pop", "I_chan", "f_per_base", "JS_vs_Drt3b",
]


def saveCsv(rows, fields, path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})


# ============================================================================
# part 10: figures
# ============================================================================

def plotJointSignature(step3_rows, path):
    """
    scatter (I_pop, I_chan) for all 6 test cases (E1, E2, E3, Drt3b-N2,
    Drt3b-N3, AbiK), color by mode-class.
    """
    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    colors = {
        "E1_fixed_ACACAC":   "#d62728",  # red
        "E2_2phase_alt":     "#ff7f0e",  # orange
        "E3_random_uniform": "#2ca02c",  # green
        "Drt3b_N2":          "#1f77b4",  # blue
        "Drt3b_N3":          "#9467bd",  # purple
        "AbiK_uniform":      "#7f7f7f",  # gray
    }
    markers = {
        "E1_fixed_ACACAC":   "o",
        "E2_2phase_alt":     "o",
        "E3_random_uniform": "o",
        "Drt3b_N2":          "s",
        "Drt3b_N3":          "s",
        "AbiK_uniform":      "X",
    }
    for r in step3_rows:
        ax.scatter(r["I_pop"], r["I_chan"], s=180,
                   color=colors.get(r["ensemble"], "k"),
                   marker=markers.get(r["ensemble"], "o"),
                   edgecolor="black", linewidth=1.0,
                   label=r["ensemble"])
        ax.annotate(r["ensemble"], (r["I_pop"], r["I_chan"]),
                    xytext=(7, 7), textcoords="offset points",
                    fontsize=9)
    ax.set_xlabel(r"$I_\mathrm{struct}^\mathrm{pop}$ (population MI, X = template realization) [bits]")
    ax.set_ylabel(r"$I_\mathrm{struct}^\mathrm{chan}$ (channel-as-X MI, vs 5-channel ensemble) [bits]")
    ax.set_title("Test G2 — joint dual-observable apparatus signature\n"
                 "(L=6, ε=0.01, n=5000)")
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="black", lw=0.5, alpha=0.5)
    ax.axvline(0, color="black", lw=0.5, alpha=0.5)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plotSeparationHeatmap(names, js_matrix, path):
    """heatmap of the 5x5 JS-divergence matrix."""
    fig, ax = plt.subplots(figsize=(7.5, 6.0))
    im = ax.imshow(js_matrix, cmap="viridis", vmin=0, vmax=1.0, aspect="equal")
    ax.set_xticks(range(len(names)))
    ax.set_yticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, ha="right")
    ax.set_yticklabels(names)
    for i in range(len(names)):
        for j in range(len(names)):
            ax.text(j, i, f"{js_matrix[i, j]:.3f}",
                    ha="center", va="center",
                    color="white" if js_matrix[i, j] < 0.5 else "black",
                    fontsize=9)
    ax.set_title("Test G2 — 5×5 cross-mode JS divergence (bits)\n"
                 "(L=6, ε=0.01, n=5000)")
    fig.colorbar(im, ax=ax, label="JS divergence (bits)")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ============================================================================
# part 11: STEP 4 — apparatus decision document
# ============================================================================

def writeApparatusDecisionDoc(step1_rows, step2_rows, step3a_rows, step3b_rows,
                               names, js_matrix, kl_matrix, path):
    """
    write the v3 apparatus-decision markdown.

    structure:
      1. joint signature for the 6 test cases (E1, E2, E3, Drt3b-N=2,
         Drt3b-N=3, AbiK)
      2. which observable separates which pair
      3. v3 Methods recommendation
      4. v3 Results recommendation
    """
    # gather joint signature rows
    joint_rows = step3a_rows + step3b_rows  # 6 rows total
    # build a lookup
    sig = {r["ensemble"]: r for r in joint_rows}

    # decision logic: which observables distinguish which pairs?
    pairs_of_interest = [
        ("E1_fixed_ACACAC", "Drt3b_N2",
         "Drt3a as it actually exists vs Drt3b N=2 (the central v3 boundary)"),
        ("E2_2phase_alt", "Drt3b_N2",
         "Drt3a-2phase counterfactual vs Drt3b N=2"),
        ("E3_random_uniform", "Drt3b_N2",
         "Drt3a-random counterfactual vs Drt3b N=2"),
        ("E1_fixed_ACACAC", "AbiK_uniform",
         "Drt3a vs AbiK"),
        ("Drt3b_N2", "Drt3b_N3",
         "Drt3b N=2 vs N=3"),
    ]

    def diff(a, b, key):
        return abs(sig[a][key] - sig[b][key])

    text = []
    text.append("# Test G2 — apparatus repair decision (v1)")
    text.append("")
    text.append("## Background")
    text.append("")
    text.append("Test G v1 found that the I_struct apparatus, computed with X = template")
    text.append("realization (= I_struct^pop, population mutual information), gives ~0 bits")
    text.append("for Drt3a's actual fixed-ACACAC ncRNA template and ~1 bit for the 2-phase")
    text.append("alternating counterfactual. Both Drt3a and Drt3b-N=2 saturate at ~1 bit by")
    text.append("the joint sequence-level estimator (sweep 2 of v1). The v2 draft's central")
    text.append("Mode 1 / Mode 3 distinguishing claim — that I_struct grows linearly with L for")
    text.append("Drt3a and saturates at log2(N) for Drt3b — does not survive when applied")
    text.append("to Drt3a as it actually exists.")
    text.append("")
    text.append("Test G v1's resolution document offered two interpretive readings. Test G2")
    text.append("reframes them as **complementary observables**, not competing repairs:")
    text.append("")
    text.append("- **Observable 1: I_struct^pop** — population MI with X = template realization.")
    text.append("  $I_\\mathrm{struct}^\\mathrm{pop} = H(Y) - H(Y \\mid X = x)$ averaged over the")
    text.append("  template population. Bounded above by H(template ensemble). Equal to the v1")
    text.append("  measurement.")
    text.append("- **Observable 2: I_struct^chan** — channel-as-X MI with X = the channel identity.")
    text.append("  Computed against a fixed 5-channel ensemble {Drt3a-WC, Drt3b-Mode3-N=2,")
    text.append("  Drt3b-Mode3-N=3, Drt3b-Mode3-N=5, AbiK-uniform}. The focal-channel summand is")
    text.append("  $\\mathrm{KL}(P(Y \\mid C = c_\\mathrm{focal}) \\,\\|\\, \\bar P(Y))$, where")
    text.append("  $\\bar P(Y)$ is the equiprobable mixture of the 5 channels' empirical product")
    text.append("  distributions.")
    text.append("")
    text.append("The product (I_struct^pop, I_struct^chan, per-base fidelity) is the refined")
    text.append("apparatus signature.")
    text.append("")
    text.append("## Joint apparatus signature for the 6 test cases (L=6, ε=0.01, n=5000)")
    text.append("")
    text.append("| Test case            | I_pop (bits) | I_chan (bits) | per-base fidelity | JS vs Drt3b-N=2 |")
    text.append("|----------------------|-------------:|--------------:|------------------:|----------------:|")
    for r in joint_rows:
        text.append(f"| {r['ensemble']:<20s} | "
                    f"{r['I_pop']:12.4f} | {r['I_chan']:13.4f} | "
                    f"{r['f_per_base']:17.4f} | {r['JS_vs_Drt3b']:15.4f} |")
    text.append("")
    text.append("Notation:")
    text.append("- E1 = Drt3a with the single fixed ACACAC template (the actual biological case)")
    text.append("- E2 = Drt3a with a 2-phase {ACACAC, CACACA} template population")
    text.append("- E3 = Drt3a with a random uniform 6-nt template (counterfactual baseline)")
    text.append("- per-base fidelity for Mode 3 entries is mechanism fidelity")
    text.append("  $P(Y_k = (\\phi_0 + k) \\bmod N)$, not Watson-Crick fidelity")
    text.append("")
    text.append("## Pairwise observable separation")
    text.append("")
    text.append("| Pair | ΔI_pop | ΔI_chan | Δf | ΔJS | Which observable separates? |")
    text.append("|------|-------:|--------:|---:|----:|-----------------------------|")
    for a, b, label in pairs_of_interest:
        d_pop = diff(a, b, "I_pop")
        d_chan = diff(a, b, "I_chan")
        d_f = diff(a, b, "f_per_base")
        d_js = diff(a, b, "JS_vs_Drt3b")
        # separation threshold: ≥ 0.5 bits = "clear", 0.1 - 0.5 = "weak", <0.1 = "no"
        sep_obs = []
        if d_pop >= 0.5:
            sep_obs.append("I_pop")
        if d_chan >= 0.5:
            sep_obs.append("I_chan")
        if d_f >= 0.1:
            sep_obs.append("fidelity")
        if d_js >= 0.1:
            sep_obs.append("JS")
        sep_str = ", ".join(sep_obs) if sep_obs else "**none above threshold**"
        text.append(f"| {a} vs {b} | {d_pop:.3f} | {d_chan:.3f} | {d_f:.3f} | {d_js:.3f} | {sep_str} |")
    text.append("")
    text.append("Threshold conventions: ΔI ≥ 0.5 bits = clear separation;")
    text.append("Δf ≥ 0.1 or ΔJS ≥ 0.1 = clear distributional separation.")
    text.append("")

    # ---- explicit observable-separates-which-pair statement ----
    e1_b2_pop = diff("E1_fixed_ACACAC", "Drt3b_N2", "I_pop")
    e1_b2_chan = diff("E1_fixed_ACACAC", "Drt3b_N2", "I_chan")
    e1_b2_f = diff("E1_fixed_ACACAC", "Drt3b_N2", "f_per_base")
    e1_b2_js = diff("E1_fixed_ACACAC", "Drt3b_N2", "JS_vs_Drt3b")

    text.append("## The Drt3a-vs-Drt3b decision")
    text.append("")
    text.append("The central v3 question: does the dual-observable apparatus distinguish Drt3a")
    text.append("(as it actually exists, E1: fixed ACACAC) from Drt3b-N=2 (Mode 3 alternating)?")
    text.append("")
    text.append(f"- ΔI_pop  = {e1_b2_pop:.4f} bits")
    text.append(f"- ΔI_chan = {e1_b2_chan:.4f} bits")
    text.append(f"- Δf (per-base fidelity) = {e1_b2_f:.4f}")
    text.append(f"- ΔJS(output dist vs Drt3b) = {e1_b2_js:.4f}")
    text.append("")
    if e1_b2_chan >= 0.5:
        text.append("**Verdict:** the channel-as-X observable I_chan **does** distinguish Drt3a-E1")
        text.append("from Drt3b-N=2. The biological Drt3a templates a *random-template-class* output")
        text.append("(Y rows distributed similarly to the WC channel applied to a non-degenerate")
        text.append("template), not a phase-restricted alternating-output distribution. The 5-channel")
        text.append("ensemble's empirical row distributions for Drt3a-WC and Drt3b-Mode3-N=2 differ")
        text.append("substantially in the row-distribution sense, and that difference is what I_chan")
        text.append("measures.")
    elif e1_b2_chan >= 0.1:
        text.append("**Verdict (weak):** I_chan separates Drt3a-E1 from Drt3b-N=2 at a weak level")
        text.append(f"(ΔI_chan ≈ {e1_b2_chan:.2f} bits). The channel ensembles have measurably")
        text.append("different output distributions, but the gap is not as large as the Mode 1 vs")
        text.append("Mode 3 mechanistic distinction would naively suggest. Consider per-base")
        text.append("fidelity as a third-leg observable.")
    else:
        text.append("**Verdict (negative):** I_chan does NOT distinguish Drt3a-E1 from Drt3b-N=2")
        text.append(f"(ΔI_chan ≈ {e1_b2_chan:.4f} bits, below threshold). The dual-observable")
        text.append("apparatus is insufficient. A third observable — per-base fidelity")
        text.append("$f_i = P(Y_i = \\mathrm{WC}(X_i) \\mid X_i)$ for Drt3a vs the channel-fidelity")
        text.append("$P(Y_k = (\\phi + k) \\bmod N)$ for Drt3b — is required. This is the orchestration")
        text.append("plan's flagged-negative outcome.")
    text.append("")

    # ---- Drt3a vs Drt3b global statement ----
    drt3a_pop = sig["E1_fixed_ACACAC"]["I_pop"]
    drt3a_chan = sig["E1_fixed_ACACAC"]["I_chan"]
    drt3b_pop = sig["Drt3b_N2"]["I_pop"]
    drt3b_chan = sig["Drt3b_N2"]["I_chan"]
    text.append("Direct numerical comparison:")
    text.append("")
    text.append(f"- Drt3a-E1: (I_pop, I_chan) = ({drt3a_pop:.4f}, {drt3a_chan:.4f}) bits")
    text.append(f"- Drt3b-N2: (I_pop, I_chan) = ({drt3b_pop:.4f}, {drt3b_chan:.4f}) bits")
    text.append("")

    # ---- v3 Methods recommendation ----
    text.append("## Recommendation for the v3 Methods section")
    text.append("")
    text.append("The v3 Methods apparatus paragraph should:")
    text.append("")
    text.append("1. Replace the single ambiguous I_struct(X;Y) of v2 with **two explicitly-defined**")
    text.append("   observables: I_struct^pop and I_struct^chan, with formulae and stated domains")
    text.append("   of validity.")
    text.append("2. State that the joint signature (I_struct^pop, I_struct^chan) is the primary")
    text.append("   apparatus output. Per-base fidelity f is a supplementary, non-MI observable")
    text.append("   reported alongside.")
    text.append("3. Note that I_struct^pop is bounded above by H(template ensemble) and is")
    text.append("   identically zero for a fixed-template population (Drt3a-E1 case). I_struct^chan")
    text.append("   is bounded above by log2(|channel ensemble|) and is non-zero whenever the")
    text.append("   focal channel's row distribution differs from the ensemble mixture.")
    text.append("4. Define the channel ensemble explicitly: in this paper, the canonical ensemble")
    text.append("   is {Drt3a-WC, Drt3b-Mode3-N=2, Drt3b-Mode3-N=3, Drt3b-Mode3-N=5, AbiK-uniform}.")
    text.append("   Other choices are valid but must be stated.")
    text.append("")
    text.append("See `test_g2_v3_methods_paragraph.md` for a drop-in draft.")
    text.append("")

    # ---- v3 Results recommendation ----
    text.append("## Recommendation for the v3 Results section (Mode 1 vs Mode 3)")
    text.append("")
    if e1_b2_chan >= 0.5:
        text.append("The Mode 1 vs Mode 3 distinguishing claim should be **carried by I_struct^chan**,")
        text.append("not by I_struct^pop. v3 should state:")
        text.append("")
        text.append("> Drt3a (Mode 1 with biologically-encoded ACACAC template) and Drt3b (Mode 3")
        text.append("> with N=2 cyclic active site) are **distinguishable by their output row")
        text.append("> distributions** (I_struct^chan against the 5-channel ensemble:")
        text.append(f"> {drt3a_chan:.2f} bits for Drt3a vs {drt3b_chan:.2f} bits for Drt3b-N=2),")
        text.append("> not by per-template-realization population MI. The L-scaling claim of v2")
        text.append("> survives only when applied to the random-template counterfactual (E3); for")
        text.append("> the actual biological systems, the framework distinguishes them by")
        text.append("> product-distribution divergence and per-base fidelity, not by linear-in-L")
        text.append("> growth of I_struct^pop.")
    else:
        text.append("The dual-observable apparatus does not cleanly separate Drt3a-E1 from Drt3b-N=2")
        text.append(f"(ΔI_chan ≈ {e1_b2_chan:.4f} bits). v3 must:")
        text.append("")
        text.append("1. **Drop the unconditional Mode 1 vs Mode 3 distinguishing claim** for the")
        text.append("   biological Drt3a (E1) case. The L-scaling claim survives only as a")
        text.append("   counterfactual statement about the random-template ensemble (E3).")
        text.append("2. **Add a third observable** — per-base fidelity $f_i$ — as the primary")
        text.append("   distinguishing measurement for biological systems on identical scaffolds")
        text.append("   (Drt3a vs Drt3b vs AbiK). f is mechanism-determined and is non-zero exactly")
        text.append("   when WC pairing happens, regardless of template degeneracy.")
        text.append("3. State that the 'Mode 1 vs Mode 3' distinction is a **mechanism claim**, not")
        text.append("   an MI-apparatus claim, when the template is degenerate.")
    text.append("")

    # ---- 5x5 cross-mode separation matrix ----
    text.append("## 5×5 cross-mode separation matrix (JS divergence in bits)")
    text.append("")
    header = "| | " + " | ".join(names) + " |"
    sep_line = "|---|" + "|".join(["---:"] * len(names)) + "|"
    text.append(header)
    text.append(sep_line)
    for i, n_i in enumerate(names):
        row = f"| **{n_i}** | " + " | ".join(f"{js_matrix[i,j]:.4f}" for j in range(len(names))) + " |"
        text.append(row)
    text.append("")
    text.append("KL divergence (i → j) in bits:")
    text.append("")
    text.append(header)
    text.append(sep_line)
    for i, n_i in enumerate(names):
        row = f"| **{n_i}** | " + " | ".join(f"{kl_matrix[i,j]:.4f}" for j in range(len(names))) + " |"
        text.append(row)
    text.append("")

    text.append("## Provenance")
    text.append("")
    text.append("- code: `code/test_g2_dual_observable.py`")
    text.append("- step 1 csv: `results/test_g2_dual_observable_v1.csv`")
    text.append("- step 2 csv: `results/test_g2_mode_separation_matrix_v1.csv`")
    text.append("- step 3 csv: `results/test_g2_drt3a_three_ensembles_v1.csv`")
    text.append("- figures: `figures/test_g2_dual_observable_signature.png`,")
    text.append("            `figures/test_g2_mode_separation_heatmap.png`")
    text.append("- estimators: per-position plug-in MI for I_pop; KL/JS over empirical row")
    text.append("  distributions for I_chan and the cross-mode matrix.")
    text.append("- seed: `np.random.seed(42)`, `np.random.default_rng(42)`")
    text.append("- n_samples = 5000 per cell (5× v1's 1000), epsilon = 0.01")
    text.append("")

    path.write_text("\n".join(text))


# ============================================================================
# part 12: STEP 5 — v3 Methods paragraph draft
# ============================================================================

def writeMethodsParagraph(path):
    """drop-in v3 Methods paragraph defining both observables explicitly."""
    text = []
    text.append("# v3 Methods — apparatus paragraph (drop-in)")
    text.append("")
    text.append("(this paragraph replaces the single-observable apparatus definition in v2")
    text.append("Methods/Results lines 35–46 of `templating_substrates_draft_v2.md`.)")
    text.append("")
    text.append("## Apparatus definition")
    text.append("")
    text.append("We define a templating event as the tuple")
    text.append("$(X, O, S, \\Delta G; \\phi_X, \\phi_Y) \\to Y$, with notation as in v2. The framework's")
    text.append("apparatus uses two complementary mutual-information observables, evaluated")
    text.append("descriptor-relatively against the bulk-matched and structure-scrambled controls of v2.")
    text.append("")
    text.append("**Observable 1: population mutual information.** With $X$ taken as the *realization*")
    text.append("of the template (a specific 1D sequence, conformer, or module sequence drawn from a")
    text.append("template population $\\mathcal{X}$ with distribution $\\Pi(X)$), the population MI is")
    text.append("")
    text.append("$$ I_\\text{struct}^\\text{pop}(\\mathcal{X}) \\;=\\; H(\\phi_Y(Y)) \\;-\\; \\mathbb{E}_{X \\sim \\Pi}\\left[H(\\phi_Y(Y) \\mid X)\\right]. $$")
    text.append("")
    text.append("$I_\\text{struct}^\\text{pop}$ is bounded above by the data-processing inequality at")
    text.append("$H(\\Pi)$, the entropy of the template population. It is identically zero for a fixed")
    text.append("template ($|\\mathcal{X}| = 1$, e.g., a single ncRNA sequence in a cell) regardless of")
    text.append("the per-sample mechanism. $I_\\text{struct}^\\text{pop}$ measures cross-realization")
    text.append("transferable information; it is the *informational capacity of the system as an")
    text.append("encoder*.")
    text.append("")
    text.append("**Observable 2: channel-as-X mutual information.** With $X$ taken as the *identity")
    text.append("of the channel itself* (the parameterized templating apparatus, including substrate-")
    text.append("alphabet, mechanism class, and channel parameters), the channel-as-X MI against a")
    text.append("comparison ensemble $\\mathcal{C} = \\{c_1, \\ldots, c_K\\}$ of channels is")
    text.append("")
    text.append("$$ I_\\text{struct}^\\text{chan}(\\mathcal{C}) \\;=\\; \\sum_{c \\in \\mathcal{C}} P(C=c)\\, D_\\text{KL}\\!\\left(P(\\phi_Y(Y) \\mid C=c) \\,\\|\\, \\bar P(\\phi_Y(Y))\\right), $$")
    text.append("")
    text.append("where $\\bar P(\\phi_Y(Y)) = \\sum_c P(C=c) P(\\phi_Y(Y) \\mid C=c)$ is the mixture row")
    text.append("distribution. The focal-channel summand,")
    text.append("$D_\\text{KL}(P(\\phi_Y(Y) \\mid C=c_\\text{focal}) \\,\\|\\, \\bar P(\\phi_Y(Y)))$, is")
    text.append("the focal channel's contribution and is the natural per-channel signature.")
    text.append("$I_\\text{struct}^\\text{chan}$ is bounded above by $\\log_2 K$. It is non-zero")
    text.append("whenever the focal channel's product distribution differs from the ensemble")
    text.append("mixture, regardless of whether the focal channel's template population is degenerate.")
    text.append("$I_\\text{struct}^\\text{chan}$ measures *how distinguishable the channel is from")
    text.append("its alternatives at the product-distribution level*.")
    text.append("")
    text.append("**Conditions of validity.** $I_\\text{struct}^\\text{pop}$ requires a template")
    text.append("population with $|\\mathcal{X}| > 1$ to be informative; for a fixed template it is")
    text.append("zero by construction. $I_\\text{struct}^\\text{chan}$ requires a stated comparison")
    text.append("ensemble $\\mathcal{C}$ that includes the focal channel and at least one alternative;")
    text.append("its numerical value depends on the ensemble choice and that choice must be reported")
    text.append("alongside the value.")
    text.append("")
    text.append("**Joint signature.** The refined apparatus output is the triple")
    text.append("$(I_\\text{struct}^\\text{pop},\\; I_\\text{struct}^\\text{chan},\\; \\bar f)$, where")
    text.append("$\\bar f = \\mathbb{E}[\\mathbf{1}\\{Y_i = \\text{WC}(X_i)\\}]$ is the per-base mechanism")
    text.append("fidelity averaged over positions and realizations (or the appropriate channel-")
    text.append("specific per-position correctness probability for non-WC channels). $\\bar f$ is not")
    text.append("a mutual information but is a non-negotiable mechanism observable: it distinguishes")
    text.append("Drt3a's per-base WC pairing (≈ 0.99) from AbiK's chance-level output (= 0.25)")
    text.append("regardless of template degeneracy or channel-ensemble choice.")
    text.append("")
    text.append("**Mode classification with the joint signature.**")
    text.append("- Mode 1 (Drt3a): $\\bar f \\approx 1 - \\varepsilon$ regardless of template ensemble;")
    text.append("  $I_\\text{struct}^\\text{pop}$ scales linearly with $L$ when $\\mathcal{X}$ is")
    text.append("  non-degenerate (E3 case) and is zero/bounded when $\\mathcal{X}$ is degenerate")
    text.append("  (E1/E2 cases); $I_\\text{struct}^\\text{chan}$ depends on the ensemble.")
    text.append("- Mode 3 (Drt3b N=2): $\\bar f \\approx 1 - \\varepsilon$ as a *channel*-fidelity;")
    text.append("  $I_\\text{struct}^\\text{pop}$ saturates at $\\log_2 N$ regardless of $L$;")
    text.append("  $I_\\text{struct}^\\text{chan}$ separates from Mode 1 when the channel ensemble")
    text.append("  contains both.")
    text.append("- Random / non-templating (AbiK): $\\bar f = 1/|\\mathcal{A}|$;")
    text.append("  $I_\\text{struct}^\\text{pop} = 0$; $I_\\text{struct}^\\text{chan}$ is the")
    text.append("  divergence of the uniform distribution from the ensemble mixture.")
    text.append("")
    text.append("**Controls.** The bulk-matched and structure-scrambled controls of v2 apply unchanged")
    text.append("to $I_\\text{struct}^\\text{pop}$. For $I_\\text{struct}^\\text{chan}$ the analogous")
    text.append("control is $D_\\text{KL}(P(\\phi_Y(Y_\\text{scrambled})) \\,\\|\\, \\bar P(\\phi_Y(Y)))$,")
    text.append("which probes whether the channel's signature against the ensemble survives positional")
    text.append("scrambling of the focal channel's outputs.")
    text.append("")
    path.write_text("\n".join(text))


# ============================================================================
# part 13: main
# ============================================================================

def main():
    np.random.seed(42)
    rng = np.random.default_rng(42)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # output paths
    step1_csv = RESULTS_DIR / "test_g2_dual_observable_v1.csv"
    step2_csv = RESULTS_DIR / "test_g2_mode_separation_matrix_v1.csv"
    step3_csv = RESULTS_DIR / "test_g2_drt3a_three_ensembles_v1.csv"
    decision_md = RESULTS_DIR / "test_g2_apparatus_decision_v1.md"
    methods_md  = RESULTS_DIR / "test_g2_v3_methods_paragraph.md"
    fig_signature = FIGURES_DIR / "test_g2_dual_observable_signature.png"
    fig_heatmap   = FIGURES_DIR / "test_g2_mode_separation_heatmap.png"

    # safety: do NOT overwrite Test G v1 outputs (the dispatch is explicit)
    v1_files = [
        RESULTS_DIR / "test_g_drt3a_boundary_v1.md",
        RESULTS_DIR / "test_g_alternating_vs_random_template_v1.csv",
        RESULTS_DIR / "test_g_drt3a_vs_drt3b_comparison_v1.csv",
    ]
    for p in v1_files:
        if not p.exists():
            print(f"WARNING: Test G v1 output {p} not found — G2 should supplement, not replace.")

    epsilon = 0.01
    n_samples = 5000

    t_start = time.time()

    # Step 1: full v1 sweep with both observables
    step1_rows = runStep1Sweep(rng,
                               L_Ts=(2, 4, 6, 8, 12, 24),
                               epsilon=epsilon, n_samples=n_samples)
    saveCsv(step1_rows, CSV_FIELDS_STEP1, step1_csv)

    # Step 2: 5×5 cross-mode separation matrix at L=6
    step2_rows, names, js_matrix, kl_matrix = runStep2SeparationMatrix(
        rng, L=6, epsilon=epsilon, n_samples=n_samples
    )
    saveCsv(step2_rows, CSV_FIELDS_STEP2, step2_csv)

    # Step 3: Drt3a under E1/E2/E3 + Step 3b: Drt3b-N2/N3 + AbiK
    step3a_rows = runStep3Drt3aThreeEnsembles(rng, L=6, epsilon=epsilon,
                                              n_samples=n_samples)
    step3b_rows = runStep3bOtherCases(rng, L=6, epsilon=epsilon,
                                      n_samples=n_samples)
    saveCsv(step3a_rows + step3b_rows, CSV_FIELDS_STEP3, step3_csv)

    # figures
    plotJointSignature(step3a_rows + step3b_rows, fig_signature)
    plotSeparationHeatmap(names, js_matrix, fig_heatmap)

    # Step 4: apparatus decision document
    writeApparatusDecisionDoc(step1_rows, step2_rows, step3a_rows, step3b_rows,
                              names, js_matrix, kl_matrix, decision_md)

    # Step 5: drop-in v3 Methods paragraph
    writeMethodsParagraph(methods_md)

    elapsed = time.time() - t_start
    print()
    print("=" * 96)
    print(f"step 1 csv:        {step1_csv}")
    print(f"step 2 csv:        {step2_csv}")
    print(f"step 3 csv:        {step3_csv}")
    print(f"decision md:       {decision_md}")
    print(f"methods md:        {methods_md}")
    print(f"signature figure:  {fig_signature}")
    print(f"heatmap figure:    {fig_heatmap}")
    print(f"total wall-time:   {elapsed:.2f} s")
    print("=" * 96)


if __name__ == "__main__":
    main()
