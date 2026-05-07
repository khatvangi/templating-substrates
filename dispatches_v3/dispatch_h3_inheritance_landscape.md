# Dispatch H3: Inheritance carrier landscape — beyond binary copyability (Test H3)

## Working directory
`/storage/kiran-stuff/templating_framework/`

## Goal
Test H2 falsified P_H2_3, the framework's prediction that Implementations
A and B should give identical results because both lack copy mechanism.
The empirical result: at N_1=20, Implementation A gives Mode 1 wins 93%,
Implementation B gives Mode 1 wins 13%. Mode 6 wins 87% in B because
per-lineage fixed templates create lineage-level inheritance via
selection on chance-favorable founder draws, even though no individual
agent inherits anything from its parent.

This means the framework's binary copyability claim conflated two things:
individual-level copy (what Mode 1 has) and lineage-level fixation
(what Implementation B accidentally implemented). Test H3 enumerates the
inheritance-mechanism landscape with finer resolution and tests each
mechanism's evolutionary capacity across multiple selection regimes.

This is Test H3. Its purpose is to give v3 the conceptual machinery to
correctly distinguish what individual-level copy adds beyond what
lineage-level fixation already provides.

## Background

Read these first:
1. `code/test_h2_competition_sweep.py` and `results/test_h2_*_v1.csv`
   — the simulation that revealed the issue. Implementation A vs B
   divergence (max |P_A-P_B|=0.80) is the empirical basis for H3.
2. `results/test_h2_predictions_v1.md` — pre-registered predictions
   table; P_H2_3 falsified is the entry-point.
3. `paper/drafts/draft_v2.md` Discussion line 127 (the published claim
   that Mode 6 cannot develop in shared-environment competition) and
   the inheritance theorem in `framework_formal_v1.md` condition (i).
4. `framework_formal_v1.md` — inheritance theorem statement. Note
   condition (i) says "heritable copying" without specifying whether
   the heritability is individual-level or lineage-level. H3 makes
   this distinction operational.

## The inheritance mechanism landscape

Define five inheritance mechanisms ranging from no-inheritance to
full-Mode-1 copying. Each is implemented as a population-dynamics
agent with a specific reproduction rule. Within H3, "Mode 6" is no
longer a single thing — it's a label that admits at least three
sub-variants.

### Mechanism M0: pure stateless (true no-inheritance)
- Each agent's phenotype is freshly drawn from a fixed global
  distribution every generation, regardless of parent.
- No founder fixation, no mutation, no lineage state.
- This is what the framework should have meant by "no copy mechanism."

### Mechanism M1: lineage-level fixation, no mutation (Test H2 Implementation B)
- At lineage founding, draw a random template from the global
  distribution. Lineage members all share that template.
- Reproduction copies the *parent's lineage label* but the offspring's
  phenotype is the lineage's fixed template plus per-individual noise.
- Lineages compete; favorable templates get selected for.
- This is what Implementation B turned out to be.

### Mechanism M2: lineage-level fixation with branch re-draw
- Same as M1, but every reproduction event triggers a fresh template
  draw with probability r ∈ {0, 0.1, 0.5, 1}.
- At r=0, M2 reduces to M1. At r=1, M2 reduces to M0.
- Tests whether the inheritance benefit of M1 depends on how durable
  the lineage-level fixation is.

### Mechanism M3: individual-level copy, no mutation
- Reproduction copies the parent's phenotype exactly.
- No mutation: variation can only enter via the initial population.
- Tests whether individual copy without mutation adds anything beyond M1.

### Mechanism M4: individual-level copy with mutation (Mode 1 / Test H2 Implementation A)
- Reproduction copies the parent's genotype with per-position mutation
  rate μ.
- Mutations propagate to descendants.
- This is the canonical Mode 1.

The framework's claim — refined for H3 — should be: only M4 supports
unbounded fitness climbing under selection. M1 and M3 plateau because
they cannot generate new variation past initial diversity; M0 cannot
respond to selection at all. M2 with r > 0 inherits the M1 ceiling
modulated by re-draw rate.

## Pre-registered predictions

Document in script header BEFORE running.

### Per-mechanism fitness ceiling

P_H3_1: M0 mean fitness plateaus at chance level (mean fitness
        ≈ initial-distribution mean), unaffected by generation count
        or selection strength β.

P_H3_2: M1 mean fitness rises above M0 due to lineage selection but
        plateaus at the maximum fitness present in the *initial
        population's* lineage templates (no new variation enters).
        Plateau height should depend on initial population size K
        and target sharpness.

P_H3_3: M2 mean fitness interpolates between M0 and M1: at r=1,
        ≈M0; at r=0, ≈M1. Plateau height monotonically decreasing in r.

P_H3_4: M3 mean fitness should equal M1's plateau when starting from
        identical initial distributions. Both lack new variation;
        both plateau at initial-distribution best.

P_H3_5: M4 mean fitness rises monotonically toward target without
        plateau (within sim length), exceeding all other mechanisms'
        ceilings by gen 1000 in all selection regimes.

### Pairwise competition outcomes

P_H3_6: M4 vs M1 competition: M4 wins ≥95% by gen 500 across all
        N_1 ∈ {20, 80, 200} initial conditions. (This is the
        framework's preserved claim — individual copy with mutation
        beats lineage-level fixation eventually.)

P_H3_7: M4 vs M0 competition: M4 wins ≥99% at all N_1, fastest
        timescale of all pairings.

P_H3_8: M1 vs M0 competition: M1 wins by lineage selection;
        crossover N_1 (where P_M1_wins = 0.5) at moderate values.
        This is the result Test H2 stumbled on.

P_H3_9: M3 vs M1 competition: outcome ambiguous because both lack new
        variation; predicted near-tie at all N_1.

### The key conceptual claim

P_H3_10: Across all selection regimes (target sharpness β ∈ {2, 5, 10, 20}),
         only M4 reaches mean fitness > 0.95 of theoretical maximum
         within 1000 generations. M1, M2, M3 plateau below this
         threshold. M0 stays at chance level. This is the framework's
         strongest defensible claim: *the ability to generate new
         variation through copy-with-mutation is what makes individual-
         level inheritance distinct from lineage-level fixation, and
         this distinction is what allows unbounded fitness climbing.*

If P_H3_10 fails for any selection regime — if any non-M4 mechanism
reaches > 0.95 of max — the framework's individual-vs-lineage
distinction itself fails, and a deeper revision is needed.

## Concrete steps

### Step 1: Per-mechanism isolated dynamics

For each mechanism M0–M4, run isolated population dynamics:
- K = 400 agents, all of the same mechanism
- Target length L = 32, mutation rate μ = 0.01 for M4 (only)
- N_GEN = 1000 generations
- Selection regimes: β ∈ {2, 5, 10, 20} (target sharpness)
- For M2: r ∈ {0, 0.1, 0.5, 1.0}
- 30 replicates per cell

For each cell, record per-generation:
- Mean fitness
- Maximum fitness
- Variance in fitness
- Number of distinct lineages (for M1, M2)

Output: `results/test_h3_isolated_dynamics_v1.csv`.

Compute per-mechanism plateau height (mean fitness at gen 1000 averaged
over reps). This tests P_H3_1 through P_H3_5.

### Step 2: Pairwise competition matrix

Run head-to-head competitions. For each pair (M_i, M_j) with i < j:
- N_i ∈ {20, 80, 200} out of K=400; N_j = K - N_i
- 30 replicates per N_i
- N_GEN = 1000, β = 10 (canonical Test H setting)

Pairs to run: (M0,M1), (M0,M2 at r=0), (M0,M3), (M0,M4), (M1,M2 at r=0.5),
(M1,M3), (M1,M4), (M2,M4 at r=0), (M3,M4). Total: 9 pairs × 3 N_i × 30
reps = 810 sims.

Record per cell:
- P(M_i wins by 95% threshold)
- Median crossover generation
- Mean final fitness for each mechanism

Output: `results/test_h3_pairwise_competition_v1.csv`.

This tests P_H3_6 through P_H3_9.

### Step 3: Selection regime sensitivity

Run M0–M4 in isolation across β ∈ {2, 5, 10, 20} (covered in Step 1)
plus L_TARGET ∈ {16, 32, 64, 128} at fixed β=10. The latter tests
whether target length affects mechanism rankings.

Output: `results/test_h3_selection_regime_v1.csv`.

This is the discriminating test for P_H3_10.

### Step 4: Lineage tracking and effective heritability

For M1 and M2 specifically, track the effective heritability h_eff:
the correlation between parent-phenotype and offspring-phenotype
across generations. This tests whether lineage-level fixation
produces effective heritability indistinguishable from M3
(individual-level copy without mutation).

Compute h_eff at generations 50, 200, 500 for M1, M2 (r=0, 0.1, 0.5),
and M3.

Output: `results/test_h3_effective_heritability_v1.csv`.

If M1's h_eff matches M3's, the framework's revised statement should
be "lineage-level and individual-level inheritance are equivalent
when both lack new variation"; if they differ, the framework should
say something more specific.

### Step 5: v3 conceptual reformulation

Generate `results/test_h3_v3_inheritance_revision.md` containing:

1. The empirical inheritance-mechanism landscape (M0–M4 plateau heights,
   per-mechanism fitness trajectories).

2. The reformulated inheritance condition that replaces v2's
   condition (i). Drafts to test:

   Draft A: "A substrate is a primary Darwinian inheritance carrier iff
   it admits *individual-level copying with mutation* (M4-class) and
   the additional capacity-scaling, variation-preservation, and
   genotype-phenotype linkage conditions (ii)–(iv)." This rules out
   M1, M2, M3 by mechanism specification.

   Draft B: "A substrate is a primary Darwinian inheritance carrier iff
   it can generate new heritable variation under selection at a rate
   bounded away from zero." This is mechanism-agnostic but requires
   defining "new heritable variation" precisely.

   Draft C: "A substrate is a primary Darwinian inheritance carrier iff
   the population's accessible phenotype space scales unboundedly with
   simulation length." Operational, testable, but implicit about
   mechanism.

   The data should suggest which draft is most defensible.

3. The implication for Mode 6: the framework's claim should be that
   real Mode 6 (surface-position templating) implements M0-class
   inheritance, not M1-class. Whether this is biologically accurate
   for actual Mode 6 examples is a literature question, flagged for
   later. If real Mode 6 implements M1-class (lineages have stable
   templates), the framework should classify it as a weaker but
   genuine inheritance carrier, not as a non-inheritor.

4. The implication for the published v2 line 127 claim: it survives
   in spirit (M4 dominates competitive selection) but the specific
   "Mode 6 should not develop inheritance" claim becomes "Mode 6
   should not develop M4-class inheritance," which is weaker.

This document is the v3 Discussion's central revision. Draft each
of A, B, C as candidate replacement statements, with the data each
draft is best supported by.

## Output artifacts

- `code/test_h3_inheritance_landscape.py`
- `results/test_h3_isolated_dynamics_v1.csv`
- `results/test_h3_pairwise_competition_v1.csv`
- `results/test_h3_selection_regime_v1.csv`
- `results/test_h3_effective_heritability_v1.csv`
- `results/test_h3_v3_inheritance_revision.md`
- `figures/test_h3_isolated_plateau_heights.png`
- `figures/test_h3_pairwise_outcomes_heatmap.png`
- `figures/test_h3_selection_regime.png`
- `figures/test_h3_effective_heritability.png`

Update `results/CANONICAL_RESULTS.md` to add Test H3.

## Constraints

- Naming: `test_h3` prefix, `_v1` suffix.
- Do not modify Test H or H2 outputs.
- Match Test H2 parameters where overlapping: K=400, μ=0.01 (M4 only),
  L_TARGET=32, β=10 canonical.
- Wall-time budget: 5 mechanisms × {4 β + 4 L_TARGET} × 30 reps for
  Step 1+3 = 1200 sims. Step 2: 810 sims. Step 4 reuses Step 1 data.
  Total ~2000 sims. Each sim ~15 sec. Total ~8 hours wall-time.
  DISPATCH NOHUP WITH PROGRESS POLLING.

  Use the H2 nohup pattern. Write progress.txt every 50 completed sims.

## Reporting

Final reply should include:
1. The 5×5 plateau-height table (mean fitness at gen 1000 per mechanism
   per β regime).
2. Whether P_H3_10 holds across all β regimes.
3. Effective heritability comparison: M1 (lineage fixation) vs M3
   (individual copy without mutation).
4. Recommended draft (A, B, or C) for the revised inheritance condition,
   with the data justifying the choice.
5. Path to v3_inheritance_revision.md.
