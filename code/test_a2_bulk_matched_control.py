"""
test A.2 — bulk-matched biasing-participant control discrimination

Templating Substrates Framework, Test A.2
-----------------------------------------
Purpose: Verify that the framework's per-position mutual information
measurement correctly distinguishes:

  (a) a Mode 1 templating system, where each template position determines
      its own product position (positional information transfer)
  (b) a biasing-participant control system, where the SAME composition
      bias applies at every position regardless of "template" (bulk-property
      selectivity, no positional information)

Concretely:
  System A (Mode 1)              : X uniform, Y = WC(X) with error rate eps
  System B (biasing-participant) : X uniform "template", Y drawn iid from
                                   a fixed biased distribution P_bias.
                                   Y has no causal dependence on X.

The framework's prediction:
  - System A per-position MI ≈ I_per_pos = 2 - H(Y|X)        (Mode 1 signal)
  - System B per-position MI ≈ 0                              (no positional info)
  - Marginal H(Y_A) ≈ 2 bits  (uniform under uniform X + WC channel)
  - Marginal H(Y_B) < 2 bits  (biased composition, despite zero MI)

This validates the framework's "ask about information position-by-position"
methodology: the SAME composition signal that fools a bulk analyst is
correctly read as zero positional information by per-position MI.

PASS criterion:
  - System A per-position MI > 1.5 bits           (close to theoretical 1.634)
  - System B per-position MI < 0.05 bits          (close to zero)
  - System A marginal H(Y)   > 1.95 bits          (close to 2)
  - System B marginal H(Y)   < 1.4 bits           (close to ~1.36, well below 2)

This file is intentionally self-contained: the simulator and MI estimator
are duplicated from test_a1 so changes to A.1 cannot silently shift A.2's
diagnostic.
"""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless backend
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

# ---- experiment configuration (frozen at spec) -----------------------------
# the biasing-participant's fixed product distribution: heavy A bias
P_BIAS = np.array([0.7, 0.1, 0.1, 0.1])
EPSILON_MODE1 = 0.05
L = 100
N_SAMPLES = 5000


# ----------------------------------------------------------------------------
# simulators
# ----------------------------------------------------------------------------
def simulate_mode1(L, epsilon, n_samples, rng):
    """
    Mode 1 (Watson-Crick) channel.

    X uniform over {0,1,2,3}; Y = WC(X) with prob 1-eps, uniform wrong with eps.
    returns X, Y both shape (n_samples, L), dtype int8.
    """
    X = rng.integers(0, 4, size=(n_samples, L), dtype=np.int8)

    Y_correct = WC_PAIR[X]
    is_error = rng.random(size=(n_samples, L)) < epsilon

    # for error positions: pick uniformly among the 3 wrong nucleotides
    # trick: pick offset in {0,1,2}; if offset >= correct, shift up by 1 to
    # skip the correct value, leaving exactly the 3 wrong choices
    offset = rng.integers(0, 3, size=(n_samples, L), dtype=np.int8)
    Y_wrong = np.where(offset >= Y_correct, offset + 1, offset).astype(np.int8)

    Y = np.where(is_error, Y_wrong, Y_correct).astype(np.int8)
    return X, Y


def simulate_biasing_participant(L, P_bias, n_samples, rng):
    """
    biasing-participant control.

    X_bulk uniform "template" (carries no information forward).
    Y drawn iid from fixed P_bias at every position; no dependence on X_bulk.
    returns X_bulk, Y both shape (n_samples, L), dtype int8.
    """
    P = np.asarray(P_bias, dtype=float)
    X_bulk = rng.integers(0, 4, size=(n_samples, L), dtype=np.int8)
    Y = rng.choice(4, size=(n_samples, L), p=P).astype(np.int8)
    return X_bulk, Y


# ----------------------------------------------------------------------------
# per-position plug-in MI estimator (same as test_a1; duplicated for standalone)
# ----------------------------------------------------------------------------
def estimate_mutual_information(X, Y):
    """
    plug-in estimator of I(X;Y) decomposed across positions.

    for each column i, build the empirical 4x4 joint over (x_i, y_i) and
    compute I_i = sum P(x,y) log2 [ P(x,y) / (P(x) P(y)) ].
    returns (I_total, I_per_position_mean, I_per_position_array).
    """
    n_samples, L = X.shape
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
        I_i = np.sum(joint[mask] * np.log2(joint[mask] / denom[mask]))
        I_per_pos[i] = I_i
        I_total += I_i
    return I_total, I_total / L, I_per_pos


# ----------------------------------------------------------------------------
# composition-level (marginal) entropy
# ----------------------------------------------------------------------------
def marginal_entropy(Y):
    """
    entropy of the marginal Y distribution (averaging over positions and
    samples). composition-level measurement: insensitive to positional
    structure, only sees overall nucleotide frequencies.
    """
    counts = np.bincount(Y.ravel(), minlength=4)
    p = counts / counts.sum()
    p_safe = np.where(p > 0, p, EPS_CLIP)
    return float(-np.sum(p * np.log2(p_safe)))


def theoretical_per_position_info(epsilon):
    """closed form: I_per_position = 2 - H(Y|X) for Mode 1 with uniform X."""
    if epsilon <= 0.0:
        return 2.0
    H_y_given_x = (
        -(1 - epsilon) * np.log2(max(1.0 - epsilon, EPS_CLIP))
        - epsilon * np.log2(max(epsilon / 3.0, EPS_CLIP))
    )
    return 2.0 - H_y_given_x


# ----------------------------------------------------------------------------
# experiment driver
# ----------------------------------------------------------------------------
def run_experiment():
    np.random.seed(42)
    rng = np.random.default_rng(42)

    # ---- expected reference numbers (closed form) --------------------------
    expected_ipp_mode1 = theoretical_per_position_info(EPSILON_MODE1)
    expected_marginal_H_mode1 = 2.0  # uniform X + symmetric WC channel ⇒ uniform Y

    P = P_BIAS / P_BIAS.sum()
    expected_marginal_H_bias = float(
        -np.sum(P * np.log2(np.where(P > 0, P, EPS_CLIP)))
    )
    expected_ipp_bias = 0.0

    # ---- System A: Mode 1 ---------------------------------------------------
    Xa, Ya = simulate_mode1(L, EPSILON_MODE1, N_SAMPLES, rng)
    I_total_a, ipp_a, _ = estimate_mutual_information(Xa, Ya)
    H_a = marginal_entropy(Ya)

    # ---- System B: biasing-participant -------------------------------------
    Xb, Yb = simulate_biasing_participant(L, P_BIAS, N_SAMPLES, rng)
    I_total_b, ipp_b, _ = estimate_mutual_information(Xb, Yb)
    H_b = marginal_entropy(Yb)

    rows = [
        {
            "system": "mode1",
            "L": L,
            "n_samples": N_SAMPLES,
            "I_per_position_mean": ipp_a,
            "I_total": I_total_a,
            "marginal_H_Y": H_a,
            "expected_I_per_position": expected_ipp_mode1,
            "expected_marginal_H_Y": expected_marginal_H_mode1,
        },
        {
            "system": "biasing_participant",
            "L": L,
            "n_samples": N_SAMPLES,
            "I_per_position_mean": ipp_b,
            "I_total": I_total_b,
            "marginal_H_Y": H_b,
            "expected_I_per_position": expected_ipp_bias,
            "expected_marginal_H_Y": expected_marginal_H_bias,
        },
    ]
    return rows


# ----------------------------------------------------------------------------
# I/O — CSV
# ----------------------------------------------------------------------------
def save_csv(rows, path):
    fieldnames = [
        "system",
        "L",
        "n_samples",
        "I_per_position_mean",
        "I_total",
        "marginal_H_Y",
        "expected_I_per_position",
        "expected_marginal_H_Y",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


# ----------------------------------------------------------------------------
# plot: 2-panel comparison
# ----------------------------------------------------------------------------
def plot_comparison(rows, path):
    """
    left panel : per-position MI (Mode 1 vs biasing-participant)
    right panel: marginal H(Y)  (Mode 1 vs biasing-participant)
    """
    by_system = {r["system"]: r for r in rows}
    labels = ["Mode 1\n(positional templating)",
              "biasing-participant\n(bulk-property selectivity)"]
    colors = ["C0", "C3"]

    ipp_vals = [
        by_system["mode1"]["I_per_position_mean"],
        by_system["biasing_participant"]["I_per_position_mean"],
    ]
    H_vals = [
        by_system["mode1"]["marginal_H_Y"],
        by_system["biasing_participant"]["marginal_H_Y"],
    ]

    expected_ipp = [
        by_system["mode1"]["expected_I_per_position"],
        by_system["biasing_participant"]["expected_I_per_position"],
    ]
    expected_H = [
        by_system["mode1"]["expected_marginal_H_Y"],
        by_system["biasing_participant"]["expected_marginal_H_Y"],
    ]

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(13.0, 5.5))

    # ----- left: per-position MI ----------------------------------------
    x = np.arange(2)
    bars_l = ax_left.bar(x, ipp_vals, color=colors, edgecolor="black",
                         width=0.55, label="empirical")
    # theoretical reference markers (horizontal ticks at expected values)
    for xi, exp_v in zip(x, expected_ipp):
        ax_left.hlines(exp_v, xi - 0.32, xi + 0.32, colors="black",
                       linestyles="--", linewidth=1.4)

    # value labels above bars
    for xi, v in zip(x, ipp_vals):
        ax_left.text(xi, v + 0.05, f"{v:.3f} bits",
                     ha="center", va="bottom", fontsize=10)

    # PASS thresholds (annotated as faint reference lines)
    ax_left.axhline(1.5, color="green", linestyle=":", linewidth=1.0, alpha=0.6)
    ax_left.text(1.45, 1.52, "PASS floor for Mode 1 (>1.5)",
                 color="green", fontsize=8, ha="right", va="bottom")
    ax_left.axhline(0.05, color="red", linestyle=":", linewidth=1.0, alpha=0.6)
    ax_left.text(1.45, 0.07, "PASS ceiling for biasing-participant (<0.05)",
                 color="red", fontsize=8, ha="right", va="bottom")

    ax_left.set_xticks(x)
    ax_left.set_xticklabels(labels, fontsize=9)
    ax_left.set_ylabel("per-position mutual information  $I(X_i;Y_i)$  [bits]")
    ax_left.set_ylim(-0.05, 2.1)
    ax_left.set_title("Per-position MI distinguishes positional templating\n"
                      "from bulk-property selectivity", fontsize=11)
    ax_left.grid(True, axis="y", alpha=0.3)

    # interpretive annotation
    ax_left.text(
        0.02, 0.98,
        "framework reads:\n"
        "  Mode 1            → positional info present\n"
        "  biasing-participant → no positional info",
        transform=ax_left.transAxes, fontsize=8.5, va="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                  edgecolor="gray", alpha=0.85),
    )

    # ----- right: marginal H(Y) -----------------------------------------
    bars_r = ax_right.bar(x, H_vals, color=colors, edgecolor="black",
                          width=0.55, label="empirical")
    for xi, exp_v in zip(x, expected_H):
        ax_right.hlines(exp_v, xi - 0.32, xi + 0.32, colors="black",
                        linestyles="--", linewidth=1.4)

    for xi, v in zip(x, H_vals):
        ax_right.text(xi, v + 0.04, f"{v:.3f} bits",
                      ha="center", va="bottom", fontsize=10)

    # H_uniform = 2 bits reference
    ax_right.axhline(2.0, color="black", linestyle=":", linewidth=1.0, alpha=0.5)
    ax_right.text(1.45, 2.02, "uniform H(Y) = 2 bits",
                  color="black", fontsize=8, ha="right", va="bottom")

    ax_right.set_xticks(x)
    ax_right.set_xticklabels(labels, fontsize=9)
    ax_right.set_ylabel("marginal output entropy  $H(Y)$  [bits]")
    ax_right.set_ylim(0, 2.25)
    ax_right.set_title("Composition-level signal IS detectable —\n"
                       "but not via per-position MI", fontsize=11)
    ax_right.grid(True, axis="y", alpha=0.3)

    ax_right.text(
        0.02, 0.98,
        "composition reads:\n"
        "  Mode 1            → uniform output\n"
        "  biasing-participant → biased output\n"
        "(both detect the bias; only per-position MI\n"
        " correctly attributes it to NON-templating)",
        transform=ax_right.transAxes, fontsize=8.5, va="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                  edgecolor="gray", alpha=0.85),
    )

    fig.suptitle(
        f"Test A.2 — bulk-matched control discrimination "
        f"(L={L}, n_samples={N_SAMPLES})",
        fontsize=12, y=1.00,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# pass / fail logic
# ----------------------------------------------------------------------------
def evaluate_pass(rows):
    """
    PASS if all four hold:
      Mode 1 per-position MI            > 1.5
      biasing-participant per-position MI < 0.05
      Mode 1 marginal H(Y)              > 1.95
      biasing-participant marginal H(Y)   < 1.4
    """
    by_system = {r["system"]: r for r in rows}
    a = by_system["mode1"]
    b = by_system["biasing_participant"]

    checks = [
        ("Mode 1 per-position MI > 1.5",
         a["I_per_position_mean"] > 1.5,
         f"got {a['I_per_position_mean']:.4f}"),
        ("biasing-participant per-position MI < 0.05",
         b["I_per_position_mean"] < 0.05,
         f"got {b['I_per_position_mean']:.4f}"),
        ("Mode 1 marginal H(Y) > 1.95",
         a["marginal_H_Y"] > 1.95,
         f"got {a['marginal_H_Y']:.4f}"),
        ("biasing-participant marginal H(Y) < 1.4",
         b["marginal_H_Y"] < 1.4,
         f"got {b['marginal_H_Y']:.4f}"),
    ]
    return checks


# ----------------------------------------------------------------------------
# README
# ----------------------------------------------------------------------------
README_TEXT = """\
# Test A.2 — Bulk-matched Control Discrimination

## What this test does

Verifies that the framework's per-position mutual information measurement
correctly distinguishes:

(a) a Mode 1 templating system, where each template position determines
    its own product position (positional information transfer)
(b) a biasing-participant control system, where the SAME composition bias
    applies at every position regardless of "template" (bulk-property
    selectivity, no positional information)

This is the empirical test of the framework's bulk-matched control
methodology.

## Two systems compared

System A — Mode 1 (genuine templating):
    X drawn iid uniformly over {A,C,G,T}
    Y produced by Watson-Crick channel with misincorporation rate ε
        P(Y_i = wc[X_i] | X_i)        = 1 - ε
        P(Y_i = each wrong | X_i)     = ε / 3
    fixed parameters: ε = 0.05, L = 100, n_samples = 5000

System B — biasing-participant control:
    X_bulk drawn iid uniformly (carries no information forward)
    Y drawn iid from fixed P_bias = [0.7, 0.1, 0.1, 0.1] at every position
    Y has NO causal dependence on X_bulk
    fixed parameters: same L = 100, n_samples = 5000

## Framework prediction

Per-position MI:
    System A : I_per_pos ≈ 2 - H(Y|X) ≈ 1.634 bits   (Mode 1 signal)
    System B : I_per_pos ≈ 0                         (no positional info)

Marginal H(Y):
    System A : H(Y) ≈ 2 bits         (uniform output: uniform X + symm. WC)
    System B : H(Y) ≈ 1.357 bits     (biased output, set by P_bias entropy)

## The methodological point

Both systems differ in output composition: System B has a non-uniform Y
distribution (heavy A bias), System A has uniform Y. A bulk-property
analyst measuring composition alone would see "biased output" and might
naively conclude that the biasing-participant is templating. The framework
avoids this by asking about information POSITION-BY-POSITION:

    I(X_i ; Y_i)  →  zero for biasing-participant, nonzero for Mode 1

The same composition signal that fools a bulk analyst is correctly read
as zero positional information by per-position MI.

## Outputs

- `results/test_a2_results.csv` — one row per system (Mode 1 + biasing-participant),
  with columns: system, L, n_samples, I_per_position_mean, I_total,
  marginal_H_Y, expected_I_per_position, expected_marginal_H_Y
- `figures/test_a2_comparison.png` — two-panel figure:
    left  : per-position MI bars (Mode 1 vs biasing-participant)
    right : marginal H(Y) bars   (Mode 1 vs biasing-participant)
  Both panels include theoretical reference markers (dashed horizontal lines).

## PASS / FAIL

PASS if all four hold:
    Mode 1 per-position MI              > 1.5 bits   (close to theoretical 1.634)
    biasing-participant per-position MI < 0.05 bits  (close to zero, accounting
                                                      for finite-sample bias)
    Mode 1 marginal H(Y)                > 1.95 bits  (close to 2)
    biasing-participant marginal H(Y)   < 1.4 bits   (close to ~1.36, well
                                                      below 2)

## What a PASS means

The per-position MI measurement correctly identifies positional templating
(Mode 1) and correctly excludes biasing-participants (which produce
composition bias but no positional information). This validates the
bulk-matched control methodology — the framework can distinguish
"templating-like surface effects" (such as montmorillonite-style mineral
surfaces, generic catalysts, β-cyclodextrin selectivity) from genuine
templating systems.

## What a FAIL means

Either:
- The per-position MI estimator is conflating composition with positional
  structure (would need re-design)
- The biasing-participant simulation is not actually composition-biasing-only
  (implementation bug)
- The framework's specification of bulk-matched control is incoherent
  (substantive framework issue)

## Reproducibility

`np.random.seed(42)` and an explicitly seeded `default_rng(42)` are set
at the start of `run_experiment()`. Re-running the script reproduces
the table exactly.
"""


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main():
    rows = run_experiment()

    csv_path = RESULTS_DIR / "test_a2_results.csv"
    save_csv(rows, csv_path)

    fig_path = FIGURES_DIR / "test_a2_comparison.png"
    plot_comparison(rows, fig_path)

    readme_path = RESULTS_DIR / "test_a2_README.md"
    readme_path.write_text(README_TEXT)

    # ---- summary table to stdout ------------------------------------------
    a = next(r for r in rows if r["system"] == "mode1")
    b = next(r for r in rows if r["system"] == "biasing_participant")

    print("=" * 80)
    print("Test A.2 — bulk-matched control discrimination")
    print("=" * 80)
    print(f"L = {L},  n_samples = {N_SAMPLES}")
    print(f"P_bias (biasing-participant Y distribution): "
          f"[{', '.join(f'{p:.2f}' for p in P_BIAS)}]")
    print(f"Mode 1 misincorporation rate: ε = {EPSILON_MODE1}")
    print()
    print(f"{'metric':<40s} {'Mode 1':>14s} {'biasing-part.':>16s}")
    print("-" * 72)
    print(f"{'per-position MI (empirical) [bits]':<40s} "
          f"{a['I_per_position_mean']:>14.4f} {b['I_per_position_mean']:>16.4f}")
    print(f"{'per-position MI (expected)  [bits]':<40s} "
          f"{a['expected_I_per_position']:>14.4f} {b['expected_I_per_position']:>16.4f}")
    print(f"{'I_total (= L · per-position) [bits]':<40s} "
          f"{a['I_total']:>14.3f} {b['I_total']:>16.3f}")
    print(f"{'marginal H(Y) (empirical)  [bits]':<40s} "
          f"{a['marginal_H_Y']:>14.4f} {b['marginal_H_Y']:>16.4f}")
    print(f"{'marginal H(Y) (expected)   [bits]':<40s} "
          f"{a['expected_marginal_H_Y']:>14.4f} {b['expected_marginal_H_Y']:>16.4f}")
    print()

    # ---- PASS / FAIL ------------------------------------------------------
    checks = evaluate_pass(rows)
    print("PASS criteria:")
    print("-" * 72)
    all_pass = True
    for desc, ok, detail in checks:
        marker = "PASS" if ok else "FAIL"
        print(f"  [{marker}]  {desc:<55s}  ({detail})")
        if not ok:
            all_pass = False
    print("-" * 72)
    print()
    print(f"results saved to:  {csv_path}")
    print(f"figure  saved to:  {fig_path}")
    print(f"readme  saved to:  {readme_path}")
    print()
    if all_pass:
        print("OVERALL: PASS  — the framework's per-position MI correctly "
              "discriminates positional templating from bulk-property "
              "selectivity.")
    else:
        print("OVERALL: FAIL  — one or more PASS criteria violated. See above.")


if __name__ == "__main__":
    main()
