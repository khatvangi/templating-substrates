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
