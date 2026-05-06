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
