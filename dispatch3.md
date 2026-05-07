# Dispatch 3: Origin-of-life Mode 1 vs Mode 6 competition (Test H)

## Working directory
`/storage/kiran-stuff/templating_framework/`

## Goal
Test the framework's strongest causal claim: in a shared-environment
simulation where Mode 1 and Mode 6 populations compete under selection,
Mode 6 goes extinct because it lacks a complementarity-based copy
mechanism. The framework predicts that Mode 6's information-capacity
ceiling and inability to inherit variation cause its extinction within
$N$ generations regardless of initial abundance. Test H pre-registers
this prediction and runs the competition.

This is the inheritance theorem made dynamical. Test D showed each mode's
fitness ceiling under isolation; Test H shows the modes compete under
selection in shared substrate, and Mode 1 wins.

## Background

Read these first:
1. `code/test_d_v2.py` and `results/test_d_v2_*.csv` — population-dynamics
   sim with five mode populations evolving in isolation. Test H reuses the
   reproduction/mutation logic but couples the populations.
2. `paper/drafts/draft_v2.md` Results section "Only Modes 1 and 2 carry
   generational information" — Test H sharpens that claim from "in
   isolation" to "in competition."
3. `paper/drafts/draft_v2.md` Discussion (line 127): the framework's claim
   that Mode 6's exclusion from biology is causal rather than contingent
   is stated AS A FALSIFIABLE PREDICTION in the published draft —
   "an origin-of-life simulation with multiple modes equally available
   should not develop Mode 6-based inheritance." Test H is the literal
   pre-registered execution of that prediction. Failure here is a
   substantive claim against the paper, not a footnote.
4. `framework_formal_v1.md` inheritance theorem — copyability, capacity
   scaling, variation-preserves-copyability, genotype-phenotype linkage.

## Pre-registered predictions

Define these BEFORE running the simulation. Document in script header.

P_H1: Starting with equal-population Mode 1 and Mode 6 (200 agents each),
       Mode 6 frequency declines monotonically and reaches <5% of total
       population within 200 generations under selection toward a fixed
       random target sequence of length 32.

P_H2: The crossover generation (Mode 1 frequency exceeds 95%) is bounded
       above by 500 generations, and the crossover dynamics are dominated
       by Mode 1's adaptation rate, not by Mode 6's slower extinction.

P_H3: When Mode 6 starts with 90% of initial population (Mode 1 at 10%),
       Mode 1 still climbs to >95% by generation 1000. The framework's
       claim is that initial abundance does not matter — copyability
       does.

P_H4: When Mode 1 is replaced with Mode 5 (length-bounded inheritance),
       the competition is closer: Mode 5 dominates Mode 6 but plateaus
       at the Mode 5 ceiling (0.41 fitness), and may not reach >95%
       depending on Mode 6's initial abundance. Mode 5 vs Mode 1 still
       favors Mode 1 because Mode 1's ceiling is higher.

The strongest version of the framework's claim is P_H1 + P_H3. If P_H3
fails — if Mode 6 can persist at high abundance when initially abundant —
the framework's claim that copyability alone determines outcomes is wrong,
and initial conditions or selection dynamics matter more than the framework
admits.

## Concrete steps

### Step 1: Set up coupled population dynamics

Reuse Test D v2's reproduction and mutation logic. Couple populations by:
- Shared substrate pool: total population capped at K = 400.
- Selection: each generation, agents reproduce proportional to fitness
  toward a fixed random target of length 32. Both modes face the same
  target.
- Replacement: cap at K, with selection determining which agents
  reproduce next-gen.

Mode 1 logic: agents carry a length-32 genotype. Reproduction copies
genotype with mutation rate μ = 0.01 per position. Phenotype = genotype.

Mode 6 logic: agents carry a length-32 sequence representing the templated
2D surface pattern at one timepoint, but no copy mechanism. Reproduction
draws phenotype freshly from a 2D-surface generative process — the
framework's claim is that Mode 6 produces the same pattern each generation
because the surface is fixed, but inherited variation does not propagate.
Implement this as: each Mode 6 agent draws phenotype from a fixed
generative process (e.g., a fixed transition matrix over the surface);
mutations to the agent's "genotype" do not propagate to descendants
because there is no genotype that copies. Operationally: Mode 6 agents'
phenotypes are i.i.d. draws from a fixed distribution centered on a
randomly initialized template; the population mean phenotype does not
shift under selection because variation is not inheritable.

This implementation is the framework's claim about Mode 6: not that
Mode 6 produces random patterns, but that Mode 6 produces consistent
patterns the population cannot evolve away from. If this implementation
is wrong, the framework's claim is wrong.

### Step 2: Run the four competition scenarios

Each for 1000 generations, 10 replicates per scenario:

Scenario A: Mode 1 (200) vs Mode 6 (200), random initial sequences.
Scenario B: Mode 1 (40) vs Mode 6 (360) — Mode 1 starts at 10%.
Scenario C: Mode 5 (200) vs Mode 6 (200) — Mode 5 with N_modules = 8.
Scenario D: Mode 1 (200) vs Mode 5 (200) — head-to-head among
            inheritance-capable modes.

Track per generation:
- Population frequency by mode
- Mean fitness by mode
- Total population (should hit K)

### Step 3: Test predictions

For each scenario, compute:
- Generation of Mode 6 (or Mode 5) extinction (frequency < 5%)
- Generation of dominant mode crossover (>95%)
- Final fitness of dominant mode

Compare against predictions P_H1–P_H4. Report each as confirmed or
falsified with the empirical numbers.

### Step 4: Sensitivity analysis

If Mode 6 implementation matters (which it does), run two additional
implementations to test robustness:

Implementation B: Mode 6 phenotype = randomly initialized template fixed
across the lineage; mutations to genotype produce drift but no inheritable
adaptation. (Slightly weaker form of "no inheritance.")

Implementation C: Mode 6 phenotype evolves but at a heritability < 1;
mutations propagate to descendants with probability 0.5 (vs 1.0 for
Mode 1). This tests whether the framework's claim is binary
(copyable / non-copyable) or graded (copyability strength matters).

Document which implementation is being tested. The framework's
prediction is strongest under Implementation A (binary copyability);
Implementation C explores the boundary.

## Output artifacts

- `code/test_h_competition.py`
- `results/test_h_scenario_A_v1.csv`, `..._B_v1.csv`, `..._C_v1.csv`,
  `..._D_v1.csv`
- `results/test_h_summary_v1.csv` (one row per scenario × replicate × mode)
- `results/test_h_predictions_v1.md` (the prediction table with
  empirical numbers and confirmed/falsified column)
- `figures/test_h_scenario_A_freq_vs_gen.png` (population frequencies
  over time, mean ± SD across 10 reps)
- `figures/test_h_all_scenarios.png` (small multiples)

Update `results/CANONICAL_RESULTS.md`.

## Constraints

- Naming: `test_h` prefix, `_v1` suffix.
- Do not modify Test D v2 outputs.
- Mutation rate μ = 0.01, target length 32, K = 400 — match Test D v2.
- 10 replicates is minimum; if budget allows, run 30.
- Wall-time budget: 4 scenarios × 10 reps × 1000 generations × ~400
  agents. Should fit within 15-min wall-time. If it projects longer,
  dispatch nohup with progress.txt.

## Reporting

Final reply should include:
1. Pre-registered predictions table with empirical numbers and pass/fail.
2. Crossover generation for each scenario.
3. The strongest claim that survives or fails: does P_H3 (initial
   abundance doesn't matter) hold?
4. Path to predictions_v1.md.
