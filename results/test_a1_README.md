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
