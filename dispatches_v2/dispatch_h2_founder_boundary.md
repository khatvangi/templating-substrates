# Dispatch H2: Founder-loss boundary + graded copyability (Test H2)

## Working directory
`/storage/kiran-stuff/templating_framework/`

## Goal
Test H v1 found P_H3 (the strict-abundance-independence claim) failed in
1/10 reps of scenario B, with Mode 1 going stochastically extinct from
40 founders. Test H2 maps the founder-loss boundary across initial
abundances and tests whether graded-copyability implementations of
Mode 6 change the outcome. The aim is a quantitatively characterized
result for v3, not a single point estimate from 10 reps.

This is Test H2. It supersedes the strict P_H3 statement of H v1 with a
boundary-mapped claim suitable for the v3 paper's Discussion.

## Background

Read these first:
1. `code/test_h_competition.py` — Test H v1 implementation. H2 reuses
   the simulation harness with parameter sweeps.
2. `results/test_h_predictions_v1.md` — full pre-registered prediction
   set with empirical numbers, including the founder-loss interpretation
   for the rep that failed P_H3.
3. `paper/drafts/draft_v2.md` Discussion line 127: the published claim
   is "an origin-of-life simulation with multiple modes equally available
   should not develop Mode 6-based inheritance." Test H2 maps when this
   claim holds and when founder loss makes it locally false.

## Pre-registered predictions (BEFORE simulation runs)

Document in script header. Document the parameterization rule too.

P_H2_1 (founder-loss threshold): For Mode 1 vs Mode 6 competition with
        Mode 1 initial counts N_1 ∈ {20, 40, 80, 120, 160} out of K=400
        and 30 reps per cell, the probability of Mode 1 reaching >95%
        by gen 1000 increases monotonically with N_1, transitioning
        from <50% at N_1 = 20 to >95% at N_1 ≥ 80.

P_H2_2 (Implementation A robustness): The founder-loss boundary
        identified in P_H2_1 holds for Implementation A (binary
        copyability, fixed template per replicate) — the canonical
        case in v1.

P_H2_3 (Implementation B robustness): For Implementation B (per-lineage
        fixed template, no inheritance of mutations), Mode 6 should
        behave identically to Implementation A because neither implements
        copy. The boundary should match A's at all N_1 values, within
        replicate noise.

P_H2_4 (Implementation C graded copyability): For Implementation C
        (heritability h ∈ {0.0, 0.25, 0.5, 0.75, 1.0}), Mode 6 should
        increasingly resemble Mode 1 as h → 1. At h = 0.5, Mode 6
        should compete with Mode 1 but still lose because Mode 1's
        full heritability gives it a fitness-climb advantage.
        At h = 1.0, Mode 6 IS Mode 1 (the framework's claim is binary,
        so this should be a degenerate case).

P_H2_5 (graded crossover): The crossover N_1 (initial Mode 1 count where
        P(Mode 1 wins) crosses 0.5) should shift right as Mode 6's
        heritability h increases. At h = 0, the crossover N_1 is the
        founder-loss threshold from P_H2_1. At h = 1, the crossover
        is at N_1 = 200 (the 50/50 starting point) by symmetry.

The strongest test is P_H2_5: if h = 0.5 doesn't shift the crossover
significantly relative to h = 0, then framework's "binary copyability
matters" claim is challenged by the data — the system would behave
as if only ε-difference mattered, not the binary copyability.

## Concrete steps

### Step 1: Founder-loss boundary mapping (Implementation A)

Sweep Mode 1 initial count N_1 ∈ {20, 40, 60, 80, 100, 120, 160, 200,
240, 280, 320, 360, 380} out of K=400. Mode 6 takes K - N_1 initial
slots. 30 replicates per cell, 1000 generations each, target length 32,
matching v1's parameters otherwise.

For each cell, record:
- P(Mode 1 reaches >95% by gen 1000)
- Median crossover generation (across reps where Mode 1 won)
- P(Mode 1 founder-extinction event, defined as Mode 1 → 0 within
  first 50 gens)
- Mean final fitness conditional on Mode 1 winning vs losing

Output: `results/test_h2_founder_boundary_v1.csv`.

Plot: `figures/test_h2_founder_boundary.png` — P(Mode 1 wins) vs N_1,
with 95% CI shading from beta-binomial on the 30 reps. Annotate the
N_1 at which P crosses 0.5 and 0.95.

### Step 2: Implementation B (per-lineage fixed template, no mutation inheritance)

Mode 6 logic: each Mode 6 lineage gets a randomly initialized template
at founding; that template is fixed for the lineage's entire history.
Mutations to the agent's "genotype" do not propagate to descendants
(different from Implementation A where the template is fixed per
replicate). New lineages get new random templates.

Run the founder-loss sweep (N_1 ∈ {20, 40, 80, 120, 160, 200} for
brevity here) × 30 reps for Implementation B.

Output: `results/test_h2_implementation_B_v1.csv`.

The framework's prediction (P_H2_3): Implementation B's boundary should
match A's because neither has copy mechanism. Confirming this strengthens
the claim that copy mechanism, not template specifics, is what matters.

### Step 3: Implementation C (graded copyability)

For h ∈ {0.0, 0.25, 0.5, 0.75, 1.0}: Mode 6 reproduction copies the
parent's genotype to offspring with probability h, and replaces with a
fresh draw from the fixed template + noise distribution with probability
(1 - h). At h = 0, this is Implementation A; at h = 1, Mode 6 has full
copy and behaves like Mode 1 with the same mutation rate.

For each h, run the founder-loss sweep:
- N_1 ∈ {20, 40, 80, 120, 160, 200, 240, 280, 320, 360, 380}
- 30 reps per cell
- 1000 generations, K=400, target length 32

Output: `results/test_h2_implementation_C_v1.csv` (5 h values × 11 N_1
values × 30 reps = 1650 simulations).

Plot: `figures/test_h2_graded_copyability.png` — P(Mode 1 wins) vs N_1
for each h, on one axis.

The crossover N_1 vs h gives the curve P_H2_5 predicts. If h = 0.5
shifts the crossover only marginally relative to h = 0, the binary-
copyability claim is supported. If the shift scales smoothly with h,
the claim is graded — useful for the v3 Discussion's nuanced position.

### Step 4: Pre-registered prediction evaluation

For each of P_H2_1 through P_H2_5, compute the empirical observation
and report whether the prediction is confirmed, falsified, or refined
by the data. Use 95% CIs where appropriate.

Output: `results/test_h2_predictions_v1.md` with the table and per-
prediction interpretation.

### Step 5: Synthesis and v3 Discussion paragraph

Generate `results/test_h2_v3_discussion_paragraph.md` containing:
1. The empirical founder-loss threshold (a number, with CI).
2. The claim the v3 Discussion should make: "Mode 1 wins when initial
   abundance is above the founder threshold of N_1 ≈ ?? out of 400;
   below this threshold, demographic stochasticity can drive Mode 1 to
   extinction in some replicates regardless of copyability advantage."
3. The implication for the framework: the inheritance theorem holds in
   the deterministic limit; Wright-Fisher founder loss is a population-
   genetics constraint orthogonal to the framework's copyability claim.
4. The graded-copyability finding (binary vs graded distinction).

This paragraph goes into v3 Discussion to replace the v2 line 127
claim with a properly characterized version.

## Output artifacts

- `code/test_h2_competition_sweep.py`
- `results/test_h2_founder_boundary_v1.csv`
- `results/test_h2_implementation_B_v1.csv`
- `results/test_h2_implementation_C_v1.csv`
- `results/test_h2_predictions_v1.md`
- `results/test_h2_v3_discussion_paragraph.md`
- `figures/test_h2_founder_boundary.png`
- `figures/test_h2_implementation_comparison.png` (A vs B side-by-side)
- `figures/test_h2_graded_copyability.png`

Update `results/CANONICAL_RESULTS.md` to add Test H2; retain Test H v1.

## Constraints

- Naming: `test_h2` prefix, `_v1` suffix.
- Do not modify Test H v1 outputs.
- Match Test H v1 parameters exactly except for what's swept: μ_Mode1 =
  0.01, target length = 32, K = 400, N_GEN = 1000, β = 10.0,
  N_MODULES_5 = 8, EPS_NOISE_5 = 0.05, EPS_NOISE_6 = 0.05.
- 30 replicates per cell as standard for H2; if the boundary is sharp
  (Step 1's CI is tight), 30 may be enough. If diffuse (CI > 0.2),
  flag in the output and note that 100 reps would be needed for
  publication confidence.
- Wall-time budget: Step 1: 13 N_1 × 30 reps = 390 sims. Step 2: 6 ×
  30 = 180 sims. Step 3: 5 × 11 × 30 = 1650 sims. Total ~2200
  simulations. Each sim is ~15 sec at K=400, N_GEN=1000. Total
  ~9 hours wall-time. THIS NEEDS NOHUP DISPATCH WITH PROGRESS POLLING.

  Use the pattern: write the script, dispatch with `nohup ... > log.txt
  2>&1 &`, write progress.txt every N completed sims, poll status, read
  final results.

## Reporting

Final reply should include:
1. Founder-loss threshold N_1 with 95% CI.
2. Whether Implementations A and B match (P_H2_3).
3. Whether the crossover shifts smoothly with h (P_H2_5).
4. The v3 Discussion paragraph draft.
5. Path to predictions_v1.md.
