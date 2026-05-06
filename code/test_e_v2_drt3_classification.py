"""
Test E v2 -- corrected E26Q channel parameters from Deng et al. 2026
(P(dA)=0.80, P(dG)=0.20 at state_A in E26Q variant). Replaces the placeholder
used in v1.

Apply the framework's diagnostic apparatus to three biological systems on the
Drt3b protein scaffold:
  1. Drt3b WT       -- Mode 3 with N=2, predicted I_struct = log_2(2) = 1 bit
  2. Drt3b E26Q     -- degraded Mode 3 (one gate's selectivity broken)
  3. AbiK           -- same fold, no selectivity, predicted I_struct ~ 0

Model: each position is a Markov chain step. Phase X in {0, 1} drawn uniformly
(initial conformational state of the cycle). Y[i] in {0,1,2,3} = {A,C,G,T}.

For Drt3b WT and E26Q the active site cycles deterministically between two
states: state_A and state_C. Phase 0 means state_A first; phase 1 means state_C
first.

Per-state channel definitions (as biological calibration):
  WT       state_A: P(A)=0.99, others=0.0033 each
  WT       state_C: P(C)=0.99, others=0.0033 each
  E26Q     state_A: P(A)=0.80, P(G)=0.20, P(C)=P(T)=0.0001  (Deng et al. 2026)
  E26Q     state_C: P(C)=0.99, others=0.0033 each
  AbiK            : every position uniform 0.25 each, X has no causal effect
"""
import csv
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = PROJECT_ROOT / "figures"
sys.path.insert(0, str(SCRIPT_DIR))


# ----------------------------------------------------------------------------
# Per-state channel matrices
# ----------------------------------------------------------------------------
WT_STATE_A   = np.array([0.99, 0.0033, 0.0033, 0.0034])
WT_STATE_C   = np.array([0.0033, 0.99, 0.0033, 0.0034])
E26Q_STATE_A = np.array([0.80, 0.0001, 0.20, 0.0001])
E26Q_STATE_C = np.array([0.0033, 0.99, 0.0033, 0.0034])
ABIK_STATE  = np.array([0.25, 0.25, 0.25, 0.25])


def _draw_from(p, size, rng):
    """vectorized categorical draw from probability vector p, returning shape `size`."""
    cum = np.cumsum(p)
    u = rng.random(size=size)
    return np.searchsorted(cum, u).astype(np.int16)


def simulate_two_state_cycle(L, n_samples, state_A_p, state_C_p, rng):
    """generic two-state cyclic active-site simulator."""
    X = rng.integers(0, 2, size=n_samples, dtype=np.int16)  # phase
    pos = np.arange(L, dtype=np.int64)[None, :]
    intended = (X[:, None].astype(np.int64) + pos) % 2  # 0 = state_A, 1 = state_C
    Y = np.empty((n_samples, L), dtype=np.int16)
    # mass-draw for each state, then assemble
    yA = _draw_from(state_A_p, (n_samples, L), rng)
    yC = _draw_from(state_C_p, (n_samples, L), rng)
    Y = np.where(intended == 0, yA, yC).astype(np.int16)
    return X, Y


def simulate_drt3b_wt(L, n_samples, rng):
    return simulate_two_state_cycle(L, n_samples, WT_STATE_A, WT_STATE_C, rng)


def simulate_drt3b_e26q(L, n_samples, rng):
    return simulate_two_state_cycle(L, n_samples, E26Q_STATE_A, E26Q_STATE_C, rng)


def simulate_abik(L, n_samples, rng):
    X = rng.integers(0, 2, size=n_samples, dtype=np.int16)  # phase exists but unused
    Y = _draw_from(ABIK_STATE, (n_samples, L), rng)
    return X, Y


# ----------------------------------------------------------------------------
# Joint MI estimator (X-alphabet 2, Y-alphabet 4)
# ----------------------------------------------------------------------------
EPS_CLIP = 1e-12

def _entropy_from_counts(counts):
    total = counts.sum()
    if total == 0:
        return 0.0
    p = counts.astype(np.float64) / total
    nz = p[p > 0]
    return float(-(nz * np.log2(nz)).sum())


def estimate_joint_mi(X, Y, x_alphabet=2, y_alphabet=4, block_length=None):
    """plug-in joint MI: I(X; Y[0:k]) = H(Y[0:k]) - H(Y[0:k] | X)."""
    n_samples, L = Y.shape
    if block_length is None:
        block_length = min(L, 8)
    k = min(block_length, L)
    if k < 1:
        return 0.0, 0, 0
    # require y_alphabet^k * x_alphabet <= n/10 for stable plug-in
    while y_alphabet ** k * x_alphabet > n_samples / 10 and k > 1:
        k -= 1
    powers = y_alphabet ** np.arange(k, dtype=np.int64)
    y_codes = (Y[:, :k].astype(np.int64) * powers[None, :]).sum(axis=1)
    counts_y = np.bincount(y_codes, minlength=y_alphabet ** k)
    H_y = _entropy_from_counts(counts_y)
    H_y_cond = 0.0
    counts_x = np.bincount(X.astype(np.int64), minlength=x_alphabet)
    for x in range(x_alphabet):
        nx = counts_x[x]
        if nx == 0:
            continue
        mask = X == x
        cyx = np.bincount(y_codes[mask], minlength=y_alphabet ** k)
        H_y_cond += (nx / n_samples) * _entropy_from_counts(cyx)
    return H_y - H_y_cond, k, H_y


def estimate_periodicity(Y, max_lag):
    """symbol-equality autocorrelation."""
    _, L = Y.shape
    autocorr = np.zeros(max_lag, dtype=np.float64)
    for lag in range(1, max_lag + 1):
        if lag >= L:
            autocorr[lag - 1] = 0.0
            continue
        eq = (Y[:, : L - lag] == Y[:, lag:])
        autocorr[lag - 1] = float(eq.mean())
    return autocorr


# ----------------------------------------------------------------------------
# Bulk-matched control: marginal-matched, X uncorrelated with Y_bulk
# ----------------------------------------------------------------------------
def simulate_bulk_matched(L, n_samples, target_marginal, x_alphabet, rng):
    X_bulk = rng.integers(0, x_alphabet, size=n_samples, dtype=np.int16)
    Y_bulk = _draw_from(target_marginal, (n_samples, L), rng)
    return X_bulk, Y_bulk


def empirical_marginal(Y, alphabet_size=4):
    counts = np.bincount(Y.ravel(), minlength=alphabet_size)
    return counts / counts.sum()


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------
def run_experiment():
    np.random.seed(42)
    rng = np.random.default_rng(42)

    Ls = [2, 4, 10, 20, 100, 500]
    n_samples = 100_000

    systems = [
        ("Drt3b_WT",    simulate_drt3b_wt),
        ("Drt3b_E26Q",  simulate_drt3b_e26q),
        ("AbiK",        simulate_abik),
    ]

    progress_path = RESULTS_DIR / "test_e_v2_progress.txt"
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    print(f"{'system':<12s} {'L':>5s} {'I_emp':>8s} {'I_bulk':>8s} {'ratio':>9s} "
          f"{'peak_lag':>9s} {'peak_val':>9s} {'A':>5s} {'C':>5s} {'G':>5s} {'T':>5s}")
    print("-" * 95)

    with open(progress_path, "a", buffering=1) as pf:
        pf.write(f"\n# === run started at {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")
        for system_label, sim_fn in systems:
            for L in Ls:
                t0 = time.time()
                X, Y = sim_fn(L, n_samples, rng)
                target_marg = empirical_marginal(Y)
                Xb, Yb = simulate_bulk_matched(L, n_samples, target_marg, 2, rng)
                I_emp, _, _ = estimate_joint_mi(X, Y, x_alphabet=2, y_alphabet=4)
                I_bulk, _, _ = estimate_joint_mi(Xb, Yb, x_alphabet=2, y_alphabet=4)
                ratio = I_emp / max(I_bulk, EPS_CLIP)
                max_lag = min(8, L - 1) if L > 1 else 1
                if max_lag >= 1:
                    autocorr = estimate_periodicity(Y, max_lag=max_lag)
                    peak_idx = int(np.argmax(autocorr))
                    peak_lag = peak_idx + 1
                    peak_val = float(autocorr[peak_idx])
                else:
                    peak_lag = 0
                    peak_val = 0.0
                elapsed = time.time() - t0
                rows.append({
                    "system_label": system_label,
                    "L": L,
                    "n_samples": n_samples,
                    "I_struct_empirical": I_emp,
                    "I_struct_bulk_matched": I_bulk,
                    "separation_ratio": ratio,
                    "periodicity_peak_lag": peak_lag,
                    "periodicity_peak_value": peak_val,
                    "mean_marginal_A": float(target_marg[0]),
                    "mean_marginal_C": float(target_marg[1]),
                    "mean_marginal_G": float(target_marg[2]),
                    "mean_marginal_T": float(target_marg[3]),
                })
                line = (f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {system_label}, L={L}: "
                        f"I_emp={I_emp:.4f}, I_bulk={I_bulk:.4f}, ratio={ratio:.1f}x, "
                        f"peak_lag={peak_lag}, peak_val={peak_val:.3f}, "
                        f"marg_A={target_marg[0]:.3f}, elapsed={elapsed:.1f}s")
                pf.write(line + "\n")
                print(f"{system_label:<12s} {L:>5d} {I_emp:>8.4f} {I_bulk:>8.4f} "
                      f"{ratio:>9.1f} {peak_lag:>9d} {peak_val:>9.4f} "
                      f"{target_marg[0]:>5.3f} {target_marg[1]:>5.3f} "
                      f"{target_marg[2]:>5.3f} {target_marg[3]:>5.3f}")

    return rows


def save_csv(rows, path):
    fields = ["system_label", "L", "n_samples", "I_struct_empirical",
              "I_struct_bulk_matched", "separation_ratio",
              "periodicity_peak_lag", "periodicity_peak_value",
              "mean_marginal_A", "mean_marginal_C", "mean_marginal_G", "mean_marginal_T"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})


def plot_three_systems(rows, path):
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for system_label, color in [("Drt3b_WT", "#2563eb"),
                                 ("Drt3b_E26Q", "#f59e0b"),
                                 ("AbiK", "#9ca3af")]:
        pts = sorted([(r["L"], r["I_struct_empirical"]) for r in rows
                      if r["system_label"] == system_label])
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        ax.plot(xs, ys, "-o", color=color, markersize=7, linewidth=1.8,
                label=system_label)
    ax.axhline(1.0, color="black", linestyle="--", alpha=0.5, label="log_2(2) = 1 bit")
    ax.axhline(0.0, color="grey", linestyle=":", alpha=0.5)
    ax.set_xscale("log")
    ax.set_xlabel("polymer length L")
    ax.set_ylabel(r"$I_{\mathrm{struct}}(X; Y)$  [bits]")
    ax.set_title("Test E -- three biological systems on the Drt3b scaffold")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_periodicity(rows, path, L_target=100):
    """re-derive autocorr at L_target for each system."""
    rng = np.random.default_rng(7)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    for ax, (label, fn, color) in zip(axes, [
        ("Drt3b_WT", simulate_drt3b_wt, "#2563eb"),
        ("Drt3b_E26Q", simulate_drt3b_e26q, "#f59e0b"),
        ("AbiK", simulate_abik, "#9ca3af"),
    ]):
        _, Y = fn(L_target, 50000, rng)
        autocorr = estimate_periodicity(Y, max_lag=10)
        ax.plot(np.arange(1, 11), autocorr, "-o", color=color, markersize=6)
        ax.axhline(0.25, color="grey", linestyle=":", label="chance = 0.25")
        ax.axvline(2, color="red", linestyle="--", alpha=0.5, label="lag = 2")
        ax.set_xlabel("lag")
        ax.set_ylabel("P(Y_i = Y_{i+lag})")
        ax.set_title(f"{label} (L={L_target})")
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def evaluate_pass(rows):
    """PASS criteria as specified in the task."""
    crit_wt = []
    crit_e26q = []
    crit_abik = []
    wt_at_L = {r["L"]: r for r in rows if r["system_label"] == "Drt3b_WT"}
    e26q_at_L = {r["L"]: r for r in rows if r["system_label"] == "Drt3b_E26Q"}
    abik_at_L = {r["L"]: r for r in rows if r["system_label"] == "AbiK"}
    for L in [4, 10, 20, 100, 500]:
        if L in wt_at_L:
            r = wt_at_L[L]
            ok = (0.95 <= r["I_struct_empirical"] <= 1.0001) and \
                 (r["periodicity_peak_lag"] % 2 == 0) and \
                 (r["periodicity_peak_value"] > 0.95) and \
                 (r["separation_ratio"] > 100)
            crit_wt.append((L, r, ok))
        if L in e26q_at_L:
            r = e26q_at_L[L]
            wt_I = wt_at_L[L]["I_struct_empirical"]
            ok = (r["I_struct_empirical"] < wt_I) and \
                 (r["I_struct_empirical"] > 0.1) and \
                 (r["separation_ratio"] > 5)
            crit_e26q.append((L, r, ok, wt_I))
        if L in abik_at_L:
            r = abik_at_L[L]
            ok = (r["I_struct_empirical"] < 0.05) and (r["separation_ratio"] < 2)
            crit_abik.append((L, r, ok))
    c1 = all(t[-1] for t in crit_wt)
    c2 = all(t[2] for t in crit_e26q)
    c3 = all(t[-1] for t in crit_abik)
    return c1, c2, c3, crit_wt, crit_e26q, crit_abik


README_TEXT = """\
# Test E v2 -- Drt3 Biological Anchor Classification (corrected E26Q)

Test E v2 -- corrected E26Q channel parameters from Deng et al. 2026
(P(dA)=0.80, P(dG)=0.20 at state_A in E26Q variant). Replaces the placeholder
used in v1.

## Biological context

The Drt3 system (Sharma et al. 2026, Science) is a recently-discovered
bacterial anti-phage system that produces alternating poly(GT/AC) double-
stranded DNA. The complex has two reverse transcriptases:

- Drt3a uses a conserved ACACAC region of an ncRNA as template (Mode 1).
- Drt3b synthesizes the complementary poly(AC) strand WITHOUT a nucleic acid
  template, using two amino acid residues to enforce alternation (Mode 3
  with N=2).

Sharma et al. report that the E26Q point mutation breaks one gate's
selectivity (state_A accepts dG when dGTP is available). The same protein
fold appears in AbiK, where different active-site residues yield random
DNA output -- demonstrating that "a handful of residues separates random
from sequence-specific synthesis on an identical scaffold."

## Three systems compared

  Drt3b WT      -- N=2 cyclic active site, predicted I_struct = 1 bit
  Drt3b E26Q    -- degraded N=2, one gate broken, predicted I_struct < 1 bit
  AbiK          -- same fold, no selectivity, predicted I_struct ~ 0

## Channel parameterization

  WT       state_A: P(A)=0.99, others=0.0033 each
  WT       state_C: P(C)=0.99, others=0.0033 each
  E26Q     state_A: P(A)=0.80, P(G)=0.20, P(C)=P(T)=0.0001  (Deng et al. 2026)
  E26Q     state_C: P(C)=0.99, others=0.0033 each
  AbiK            : every position uniform 0.25 each, X has no causal effect

## PASS criterion

PASS if framework apparatus correctly classifies all three systems:
  1. Drt3b WT: I_struct in [0.95, 1.0] bits, periodicity peak at lag % 2 == 0
     with value > 0.95, separation ratio > 100 vs bulk-matched control.
  2. Drt3b E26Q: I_struct < WT's value (degraded), still > 0.1 bits, separation > 5.
  3. AbiK: I_struct < 0.05 bits, separation ratio < 2.

If all three pass, the framework's biological anchor is established: the
recently-discovered Drt3 system fits exactly the Mode 3 case the framework
predicted, and the apparatus distinguishes templating from biasing-participants
in genuine biological cases.
"""


def main():
    rows = run_experiment()
    save_csv(rows, RESULTS_DIR / "test_e_v2_results.csv")
    plot_three_systems(rows, FIGURES_DIR / "test_e_v2_three_systems.png")
    plot_periodicity(rows, FIGURES_DIR / "test_e_v2_periodicity.png")
    (RESULTS_DIR / "test_e_v2_README.md").write_text(README_TEXT)

    c1, c2, c3, crit_wt, crit_e26q, crit_abik = evaluate_pass(rows)

    print("\n" + "=" * 90)
    print("Drt3b WT classification:")
    for L, r, ok in crit_wt:
        flag = "OK" if ok else "FAIL"
        print(f"  [{flag}] L={L:4d}: I={r['I_struct_empirical']:.4f}, "
              f"peak_lag={r['periodicity_peak_lag']}, "
              f"peak_val={r['periodicity_peak_value']:.4f}, "
              f"sep={r['separation_ratio']:.1f}x")
    print(f"  -> {'PASS' if c1 else 'FAIL'}")
    print()
    print("Drt3b E26Q classification:")
    for L, r, ok, wt_I in crit_e26q:
        flag = "OK" if ok else "FAIL"
        print(f"  [{flag}] L={L:4d}: I={r['I_struct_empirical']:.4f} "
              f"(WT was {wt_I:.4f}), sep={r['separation_ratio']:.1f}x")
    print(f"  -> {'PASS' if c2 else 'FAIL'}")
    print()
    print("AbiK classification:")
    for L, r, ok in crit_abik:
        flag = "OK" if ok else "FAIL"
        print(f"  [{flag}] L={L:4d}: I={r['I_struct_empirical']:.4f}, "
              f"sep={r['separation_ratio']:.1f}x")
    print(f"  -> {'PASS' if c3 else 'FAIL'}")
    print()
    print("=" * 90)
    overall = c1 and c2 and c3
    print(f"OVERALL: {'PASS' if overall else 'FAIL'}")
    print("=" * 90)


if __name__ == "__main__":
    main()
