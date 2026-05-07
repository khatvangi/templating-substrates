# Dispatch H5: M2 sweet spot characterization — variation-rate principle (Test H5)

## Working directory
`/storage/kiran-stuff/templating_framework/`

## Goal
Test H3 found an unanticipated sweet spot in M2 (lineage fixation with
re-draw rate r): plateau height at r=0.10 is 0.53, exceeding both
r=0.00 (pure fixation, ~0.47) and r=0.50 (~0.42). The framework's
intuition that re-draws should hurt M1 monotonically (P_H3_3) was
wrong. H5 characterizes this finding with finer resolution and tests
whether the sweet spot location is universal across selection regimes
or system-specific.

The conceptual stakes: if the sweet spot location depends on the
mutation rate of the comparison M4 (μ = 0.01 in H3), the framework's
Draft B revision (variation-rate based, mechanism-agnostic) gets
empirical teeth — lineage-level re-draws and individual-level mutation
are interchangeable variation sources at matched rates. If the sweet
spot is at μ ≈ 0.01 regardless of selection regime, that's a
unification result. If it shifts with selection sharpness β or target
length L, the principle is more nuanced.

This is the single most underexplored finding in the project. The
dispatch summary explicitly flagged it as worth a follow-up.

## Background

Read these first:
1. `code/test_h3_inheritance_landscape.py` — simulation harness; H5
   reuses it.
2. `results/test_h3_isolated_dynamics_v1.csv` — M2 plateau heights at
   r ∈ {0.00, 0.10, 0.50, 1.00} that motivated this dispatch.
3. `results/test_h3_v3_inheritance_revision.md` §1 — the plateau
   table showing the r=0.10 sweet spot.
4. Eigen's quasispecies model and Kimura's optimal mutation rate work
   for theoretical context. The sweet spot is plausibly an
   exploration-exploitation optimum analogous to optimal mutation
   rate in classical population genetics. If this analogy holds, the
   sweet spot should be at r ≈ μ_M4.

## Pre-registered predictions

Document in script header.

P_H5_1 (sweet spot fine structure): The plateau-height-vs-r curve at
        β=10, L=32 has a single global maximum in the range r ∈
        [0.01, 0.30]. The location of this maximum is r* ∈
        [0.05, 0.15] — the H3 r=0.10 finding refined.

P_H5_2 (rate-matching to M4 mutation): The sweet spot location r* is
        comparable to M4's per-genome mutation rate (μ × L = 0.01 ×
        32 = 0.32 expected mutations per genome per generation,
        though r is per-reproduction not per-position). If the
        analogy holds, r* should be at the value where M2's expected
        per-reproduction template-replacement rate matches M4's
        per-reproduction mutation impact.

P_H5_3 (β invariance): The sweet spot r* is approximately invariant
        across β ∈ {2, 5, 10, 20}. If r* shifts substantially (e.g.,
        moves from 0.05 at β=2 to 0.20 at β=20), the variation-rate
        principle is contingent on selection regime.

P_H5_4 (L invariance): The sweet spot r* is approximately invariant
        across L_TARGET ∈ {16, 32, 64, 128}. If r* tracks L (e.g.,
        scales as 1/L), the principle is more naturally stated in
        terms of per-position rates.

P_H5_5 (M2 plateau height at r*): The plateau height at r* across all
        regimes lies between M1 (no variation) and M4 (full
        individual-level mutation). M2 cannot exceed M4's plateau —
        if it does, the framework needs a deeper revision because
        lineage-level re-draws would then be a *better* variation
        source than individual mutation.

The strongest test is P_H5_3 + P_H5_5 jointly. If r* is regime-
invariant and M2's plateau at r* sits cleanly between M1 and M4, the
variation-rate principle is robust and v3's Draft B framing is
empirically supported. If r* shifts with regime, Draft A's mechanism-
specific framing is the safer position.

## Concrete steps

### Step 1: Fine-grained r sweep at β=10, L=32

Sweep r ∈ {0.00, 0.01, 0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25,
0.30, 0.50, 0.75, 1.00}. 30 replicates per cell, K=400, N_GEN=1000,
EPS_NOISE_M2 = 0.05 (matching H3).

For each r, record per-generation mean fitness, max fitness, lineage
count. Compute plateau height (mean fitness at gen 1000 averaged over
30 reps) and 95% CI.

Output: `results/test_h5_r_sweep_v1.csv`.

Identify r* (the r value with maximum plateau height).

### Step 2: Selection regime sensitivity (P_H5_3)

For β ∈ {2, 5, 10, 20}, run the r sweep at the four β values
restricted to r ∈ {0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20} (the
neighborhood of the H3 sweet spot). 30 reps per cell.

Output: `results/test_h5_beta_sensitivity_v1.csv`.

Identify r*(β) and check whether r* shifts with β. Plot r*(β) vs β.

### Step 3: Target length sensitivity (P_H5_4)

For L_TARGET ∈ {16, 32, 64, 128}, run the r sweep at fixed β=10,
restricted to r ∈ {0.02, 0.05, 0.10, 0.15, 0.20}. 30 reps per cell.

Output: `results/test_h5_L_sensitivity_v1.csv`.

Identify r*(L) and check whether r* shifts with L. Plot r*(L) vs L.

### Step 4: M4 mutation rate sweep for comparison (the rate-matching test)

For M4 alone (no M2 competition), sweep μ ∈ {0.001, 0.005, 0.01, 0.02,
0.05, 0.10}. K=400, N_GEN=1000, β=10, L_TARGET=32. 30 reps per cell.

Compute M4 plateau height vs μ. Identify μ* (M4 mutation rate at peak
plateau height). Compare μ* to r*.

If μ* ≈ r* (after appropriate per-genome vs per-reproduction
normalization), the rate-matching analogy is supported and lineage
re-draw and individual mutation are equivalent variation sources at
matched rates. If μ* and r* differ significantly, the analogy
breaks and the two mechanisms are not interchangeable.

Output: `results/test_h5_M4_mutation_sweep_v1.csv`.

### Step 5: M2 vs M4 head-to-head at matched variation rates

Run M2 (at r=r*) vs M4 (at μ=μ*) competition with N_M2 ∈ {20, 200, 380}
out of K=400. 30 reps per cell. β=10, L=32, N_GEN=1000.

If M2 at r* matches M4 at μ* in pairwise competition (~50% win rate at
N_M2=200), the variation rates are functionally equivalent and the
mechanism distinction is empirically meaningless at matched rates. If
M4 wins ≥80% even at matched rates, individual-level mutation has an
intrinsic advantage beyond rate.

Output: `results/test_h5_M2_vs_M4_matched_v1.csv`.

This is the cleanest test of Draft B vs Draft A. If matched rates give
near-tie, Draft B (rate-based) is supported. If M4 still dominates,
Draft A (mechanism-specific) holds.

### Step 6: v3 statement

Generate `results/test_h5_v3_statement.md` containing:
1. The empirical sweet spot location r* with regime sensitivity.
2. The rate-matching result (M2 at r* vs M4 at μ*).
3. Interpretation: which of Draft A vs Draft B does H5 support, and
   how strongly.
4. A drop-in paragraph for v3 Discussion connecting M2_r=0.10 to
   the variation-rate principle, with the H5 numbers filled in.

If H5 supports Draft B (rate-matching wins), the v3 hybrid framing
(Draft A in theorem, Draft B in Discussion) becomes well-supported
because Draft B is now the empirical generalization of Draft A.

If H5 supports Draft A (M4 wins even at matched rates), the hybrid
framing weakens: Draft B is conjectural rather than supported, and
v3 might lean harder on Draft A.

## Output artifacts

- `code/test_h5_sweet_spot.py`
- `results/test_h5_r_sweep_v1.csv`
- `results/test_h5_beta_sensitivity_v1.csv`
- `results/test_h5_L_sensitivity_v1.csv`
- `results/test_h5_M4_mutation_sweep_v1.csv`
- `results/test_h5_M2_vs_M4_matched_v1.csv`
- `results/test_h5_v3_statement.md`
- `figures/test_h5_r_sweep_curve.png` (plateau vs r at β=10, L=32)
- `figures/test_h5_regime_sensitivity.png` (r*(β) and r*(L))
- `figures/test_h5_rate_matching.png` (M2 plateau vs r overlaid with
  M4 plateau vs μ on shared axis)
- `figures/test_h5_M2_vs_M4_matched.png` (head-to-head at matched
  rates)

Update `results/CANONICAL_RESULTS.md` to add Test H5.

## Constraints

- Naming: `test_h5` prefix, `_v1` suffix.
- Match H3 parameters everywhere except what's swept.
- Wall-time budget: Step 1: 14 r × 30 reps = 420 sims. Step 2:
  4 β × 7 r × 30 reps = 840 sims. Step 3: 4 L × 5 r × 30 reps =
  600 sims. Step 4: 6 μ × 30 reps = 180 sims. Step 5: 3 N_M2 ×
  30 reps = 90 sims. Total ~2100 sims at ~15-30 sec each = 1-2 hours.
  DISPATCH NOHUP WITH PROGRESS POLLING.

  Step 3's L=128 cells will be slower than L=16 cells; use 30 sec/sim
  conservative estimate.

## Reporting

Final reply should include:
1. r* at canonical β=10, L=32 with 95% CI.
2. r*(β) and r*(L) sensitivity tables.
3. The rate-matching outcome: does M2 at r* tie M4 at μ*?
4. Recommended Draft A vs Draft B verdict from H5.
5. Path to v3_statement.md.
