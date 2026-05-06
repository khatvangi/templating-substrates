# Test A.2 — Bulk-matched Control Discrimination

## What this test does

Verifies that the framework's diagnostic apparatus correctly distinguishes a
genuine Mode 1 template from a *biasing participant*: a system that produces
the SAME single-nucleotide output composition (bulk-matched), but with no
positional structure between X and Y.

## Why this matters

Output composition alone (the "bulk" statistic) cannot tell you whether bias
came from positional templating or from uniform compositional bias. The
framework's claim is that I_struct(X;Y), evaluated position by position,
collapses to ~ 0 for the bulk-matched control while staying L * (2 - H_eps) for
genuine templating. A.2 stresses this discrimination directly.

## Two systems compared

System 1 -- biased Mode 1 (genuine templating):
    X drawn iid from biased distribution pi over {A,C,G,T}
    Y produced by WC channel with misincorporation rate eps
    Output marginal Q(a) = sum_x pi(x) * P_chan(Y=a | X=x)

System 2 -- bulk-matched biasing participant:
    X_bulk drawn iid from pi (carries no information forward)
    Y_bulk drawn iid from Q at each position, INDEPENDENT of X_bulk
    Single-nucleotide composition matches System 1 exactly (in expectation),
    but P(y_i | x_i) = P(y_i), so I(x_i; y_i) -> 0 (up to finite-sample bias)

## Bulk-matching logic

For each pi and eps we compute Q analytically (no sampling) using the
involutivity of the Watson-Crick map (wc[wc[a]] = a):

    Q(a) = pi(wc[a]) * (1 - eps) + (1 - pi(wc[a])) * (eps / 3)

So Q is exactly the System 1 output marginal (in expectation), and
verify_bulk_match confirms the empirical Y compositions agree to the
sampling-error floor (~ 0.001-0.005 per nucleotide at our sample sizes).

## Sweeps

3 template biases:
    pi_uniform = (0.25, 0.25, 0.25, 0.25)
    pi_AT_skew = (0.40, 0.10, 0.10, 0.40)
    pi_GC_skew = (0.10, 0.40, 0.40, 0.10)
4 misincorporation rates: eps in {0.01, 0.05, 0.10, 0.25}
3 lengths:                L in {25, 100, 500}
n_samples = 5000 per cell. Total: 36 cells.

## Outputs

- results/test_a2_results.csv -- one row per (pi, eps, L) cell
- results/test_a2_progress.txt -- appended in real time during the run
- figures/test_a2_separation.png -- bar chart, log y-axis: I_sys1 vs I_bulk
- figures/test_a2_marginal_check.png -- scatter of I_bulk per position vs
  marginal_match_max_dev. I_bulk per position should hover near zero
  (well below the 0.05-bit PASS threshold) regardless of how good the
  marginal match is.

## Note on the uniform-pi edge case

For pi_uniform the System 1 output marginal is itself uniform at every eps
(uniform template + symmetric WC channel). So Q is uniform, both systems
have uniform output composition, and the "bulk-match" is automatic. In
this edge case System 1 still gives I = L * (2 - H_eps) (just like Test A.1)
and System 2 gives I ~ 0, but the discrimination test isn't really being
exercised because there is nothing for an analyst to be fooled by. The
interesting cells are pi_AT_skew and pi_GC_skew -- the PASS criterion
excludes pi_uniform for that reason.

## PASS / FAIL

PASS if for every (pi != pi_uniform, every eps, L >= 100):
    I_bulk_per_position    < 0.05 bits        (apparatus rules out the impostor)
    separation_ratio       > 20x              (strong discrimination)
    marginal_match_max_dev < 0.015            (bulk-match is honest)

FAIL otherwise -- the script's stdout localizes which cell(s) failed and which
of the three sub-criteria each failure tripped.

## What a PASS means

I_struct(X; Y) at the position level discriminates positional templating
from same-composition uniform bias. The diagnostic apparatus works as
specified by the framework and is ready for tougher cases (e.g., Test B
on Mode 3 capacity prediction).

## What a FAIL means

Either the bulk-matching is dishonest (mdev too large -- Q computed wrong
or simulator broken), or the per-position MI estimator is biased enough
to be confused by composition alone (ipp_bulk too large), or the
discriminator is too weak (ratio too small). Stdout pinpoints which.

## Reproducibility

np.random.seed(42) and an explicitly seeded default_rng(42) are set
at the start of run_experiment(). Re-running the script reproduces
the table exactly.

## Verdict (this run)

Test A.2 PASSED. Across the 24 non-uniform cells (pi_AT_skew, pi_GC_skew,
all eps, all L >= 100), separation ratios were 450-1500x and marginal-match
deviations were below 0.005. The diagnostic apparatus correctly distinguishes
Mode 1 templating from bulk-property biasing participants.
