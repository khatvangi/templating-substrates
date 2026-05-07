# Dispatch H4: Long-horizon M4 vs M1 convergence (Test H4)

## Working directory
`/storage/kiran-stuff/templating_framework/`

## Goal
Test H3 found P_H3_6 falsified: M4 vs M1 win rate is only 83% at equal
initial population (N_M1=200) at gen=1000. The question is whether this
is finite-time artifact (M4 wins eventually but the 1000-gen window is
too short) or genuine M1 stability (M1 reaches a plateau that prevents
M4 from displacing it within any horizon). Test H4 extends the horizon
to 5000 generations and tracks M4's convergence rate to settle this.

This is the empirical hook for v3's defense of the framework's
mechanism-specific theorem statement (Draft A). If M4's win rate
converges to ≥99% by gen=5000, P_H3_6 was a finite-time artifact and
v3 can state "M4 dominates eventually" without qualification. If M4's
win rate plateaus below 99% at long horizons, M1 has unexpected
durability and v3 needs to acknowledge the bound.

## Background

Read these first:
1. `code/test_h3_inheritance_landscape.py` — the simulation harness.
   H4 reuses it with extended horizon.
2. `results/test_h3_pairwise_competition_v1.csv` — the M4 vs M1 data
   from H3 at N_M1 ∈ {20, 80, 200} × 30 reps × 1000 gens.
3. `results/test_h3_v3_inheritance_revision.md` — the inheritance
   reformulation H4 is sharpening.

## Pre-registered predictions

Document in script header.

P_H4_1 (convergence): At N_M1=200 (equal start), M4 win rate at
        gen=5000 is ≥99% (vs gen=1000's 83%). The deficit is finite-time
        artifact, not stable-coexistence.

P_H4_2 (M1 plateau): M1's mean fitness within 30 reps plateaus by
        gen=500 at ≤0.55 (matching H3 plateau height). M1 does not
        keep climbing past gen=500. (This isolates the artifact: if
        M1 plateaus while M4 keeps climbing, the eventual outcome is
        determined.)

P_H4_3 (M4 climbing rate): M4's mean fitness at gen=5000 reaches ≥0.98
        across all initial conditions (matching M4 isolated-dynamics
        plateau in H3 Step 1). Cross-check: M4's per-generation
        improvement rate slows as it approaches the plateau, which
        should match the isolated-dynamics curve.

P_H4_4 (crossover generation): The median crossover generation
        (where M4 frequency exceeds M1 frequency) at N_M1=200 is
        identifiable, with 95% CI bounded by gen=3000.

P_H4_5 (small N_M1 boundary): At N_M1=20, P_H3_6 was 100% (M4
        wins). H4 should preserve this at gen=5000. If it doesn't —
        if M1 ever wins at small initial population given long enough —
        that's a deeper issue worth flagging.

The strongest test is P_H4_1. If M4 wins ≥99% at gen=5000, the v3
framework is intact and the 83% number was finite-time noise. If M4
wins <99%, M1 has unexpected long-horizon stability and v3 needs to
acknowledge it.

## Concrete steps

### Step 1: Extended-horizon M4 vs M1 sweep

Sweep:
- N_M1 ∈ {20, 80, 200} out of K=400 (matching H3 cells)
- 30 replicates per cell
- N_GEN = 5000 (5x H3's 1000)
- All other parameters match H3 exactly: μ_M4 = 0.01, β = 10,
  L_TARGET = 32, EPS_NOISE_M1 = 0.05

For each cell, record per-generation:
- Population frequency by mechanism
- Mean fitness by mechanism
- Maximum fitness by mechanism

Output: `results/test_h4_long_horizon_v1.csv`.

### Step 2: Convergence rate analysis

For each replicate at N_M1=200, identify:
- Crossover generation (M4 frequency exceeds M1 frequency for the
  first time and stays above)
- Generation at which M4 reaches ≥95% of population
- Mean fitness at gen 1000, 2000, 3000, 5000 for each mechanism

Output: `results/test_h4_convergence_v1.csv`.

Plot: `figures/test_h4_convergence_curves.png` showing M4 frequency
trajectories with 30-rep mean and 95% CI shading, plus the H3
gen=1000 horizon marked as a vertical line.

### Step 3: M1 plateau check

Compute per-replicate M1 fitness trajectory at N_M1=200. Confirm M1
plateaus by gen=500 (H3 prediction). Check whether M1's plateau height
varies across replicates (driven by initial-population variance) and
whether high-plateau M1 replicates are the ones where M4's win is
delayed.

Output: `results/test_h4_m1_plateau_v1.csv`.

This isolates the mechanism: M1 wins (or delays M4's win) when its
initial population happens to draw a high-fitness lineage by chance.
M4 wins eventually because mutation lets it climb past M1's
chance-distributed ceiling.

### Step 4: Verdict and pre-registered prediction evaluation

For each P_H4_1 through P_H4_5, evaluate empirically and report
confirmed/refined/falsified with the data.

Output: `results/test_h4_predictions_v1.md`.

### Step 5: v3-ready statement

Generate `results/test_h4_v3_statement.md` containing two short
paragraphs for v3:

1. **Discussion paragraph** (next to the M2_r=0.10 sweet spot
   discussion): "M4 displaces M1 at long horizons but the displacement
   is gradual when M1's initial population happens to draw high-fitness
   lineages. At gen=1000 (Test H3), M4 wins 83% at equal start; at
   gen=5000, M4 wins XX% — confirming finite-time artifact rather than
   stable coexistence."

2. **Methods/Supplementary note**: "The H3 finite-time win-rate (83%)
   reflects M4's gradual climbing from a flat-fitness initial population
   versus M1's chance-favorable lineage selection. At extended horizon
   (5000 generations, Test H4), M4's win rate converges to XX% across
   all tested N_M1, supporting the asymptotic dominance claim of
   condition (i)."

Fill in the actual XX values from Step 2.

## Output artifacts

- `code/test_h4_long_horizon.py`
- `results/test_h4_long_horizon_v1.csv`
- `results/test_h4_convergence_v1.csv`
- `results/test_h4_m1_plateau_v1.csv`
- `results/test_h4_predictions_v1.md`
- `results/test_h4_v3_statement.md`
- `figures/test_h4_convergence_curves.png`
- `figures/test_h4_m1_plateau_distribution.png`

Update `results/CANONICAL_RESULTS.md` to add Test H4.

## Constraints

- Naming: `test_h4` prefix, `_v1` suffix.
- Do not modify any earlier outputs.
- Match H3 parameters exactly except N_GEN (1000 → 5000).
- Wall-time budget: 3 cells × 30 reps × 5000 gens × K=400 = 90 sims at
  ~75 sec/sim = ~115 minutes. DISPATCH NOHUP WITH PROGRESS POLLING.

  Use the H2/H3 nohup pattern. Write progress.txt every 10 completed
  sims (more frequent than H3's 50 because total is smaller).

## Reporting

Final reply should include:
1. M4 win rate at gen=5000 for each N_M1 (with 95% CIs).
2. Whether P_H4_1 confirms (M4 wins ≥99% at N_M1=200).
3. Median crossover generation at N_M1=200.
4. The v3 Discussion drop-in paragraph filled with empirical numbers.
5. Path to v3_statement.md.
