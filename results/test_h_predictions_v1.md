# Test H -- Pre-registered prediction results (v1)

Pre-registered in `code/test_h_competition.py` header BEFORE simulation runs.
Mode 6 implementation: A (binary copyability, fixed-template noise draw) per dispatch3.md §Step 1.

## Predictions table

| ID    | Verdict     | Empirical observation |
|-------|-------------|------------------------|
| P_H1  | **CONFIRMED** | 10/10 reps reach <5% by gen 200; median extinction gen across reps that reached <5% = 3 (range 2-8); reps reaching <5% by end of run = 10/10 |
| P_H2  | **CONFIRMED** | 10/10 reps reach Mode1>95% by gen 500; median crossover gen = 3 (range 2-8) |
| P_H3  | **FALSIFIED** | 9/10 reps reach Mode1>95% by gen 1000; final Mode1 frequency mean = 0.900, median crossover gen (where reached) = 3 |
| P_H4  | **CONFIRMED** | scenario C final freqs: Mode5=1.000, Mode6=0.000 (Mode5 dominates Mode6: True); scenario D final freqs: Mode1=1.000, Mode5=0.000 (Mode1 dominates Mode5: True) |

## Prediction statements

### P_H1 -- CONFIRMED

**Statement:** Scenario A: Mode 6 frequency declines to <5% within 200 generations in all 10 replicates.

**Empirical:** 10/10 reps reach <5% by gen 200; median extinction gen across reps that reached <5% = 3 (range 2-8); reps reaching <5% by end of run = 10/10

**Notes:** selection toward fixed random target len 32, K=400, 10 reps

### P_H2 -- CONFIRMED

**Statement:** Scenario A: Mode 1 reaches >95% frequency by generation 500 in all reps.

**Empirical:** 10/10 reps reach Mode1>95% by gen 500; median crossover gen = 3 (range 2-8)

**Notes:** strongest version: bounded by 500 gens; weaker forms (by 1000) reported in summary csv

### P_H3 -- FALSIFIED

**Statement:** Scenario B (Mode 1 starts at 10%): Mode 1 reaches >95% by generation 1000 in all reps. THIS IS THE STRONGEST VERSION OF THE FRAMEWORK'S CLAIM (initial abundance does not matter).

**Empirical:** 9/10 reps reach Mode1>95% by gen 1000; final Mode1 frequency mean = 0.900, median crossover gen (where reached) = 3

**Notes:** if FALSIFIED, framework's copyability-only claim is wrong

### P_H4 -- CONFIRMED

**Statement:** Scenario C (Mode 5 vs Mode 6): Mode 5 dominates Mode 6; Scenario D (Mode 1 vs Mode 5): Mode 1 dominates Mode 5.

**Empirical:** scenario C final freqs: Mode5=1.000, Mode6=0.000 (Mode5 dominates Mode6: True); scenario D final freqs: Mode1=1.000, Mode5=0.000 (Mode1 dominates Mode5: True)

**Notes:** Mode 5 may not reach 95% because of its 0.41 fitness ceiling

## Crossover generation per scenario

| Scenario | Dominant mode | Reps reaching >95% | Median crossover gen | Range |
|----------|---------------|---------------------|-----------------------|--------|
| A | Mode1 | 10/10 | 3 | 2--8 |
| B | Mode1 | 9/10 | 3 | 3--12 |
| C | Mode5_N8 | 10/10 | 4 | 2--36 |
| D | Mode1 | 10/10 | 5 | 4--6 |

## Per-scenario final frequency (mean across 10 replicates)

| Scenario | Mode | Mean final freq | SD final freq | Mean final fitness |
|----------|------|------------------|----------------|---------------------|
| A | Mode1 | 1.0000 | 0.0000 | 0.9718 |
| A | Mode6_implA | 0.0000 | 0.0000 | nan |
| B | Mode1 | 0.9000 | 0.3000 | 0.9724 |
| B | Mode6_implA | 0.1000 | 0.3000 | 0.3384 |
| C | Mode5_N8 | 1.0000 | 0.0000 | 0.4162 |
| C | Mode6_implA | 0.0000 | 0.0000 | nan |
| D | Mode1 | 1.0000 | 0.0000 | 0.9715 |
| D | Mode5_N8 | 0.0000 | 0.0000 | nan |

## Reproducibility

- Module seed: `np.random.seed(42)` and `np.random.default_rng(42)`
- Per-replicate seeds: `42 + r` for `r in {0..9}`
- Mode 6 fixed template per replicate: derived from the per-replicate seed
- K = 400, L_TARGET = 32, N_GEN = 1000, BETA = 10.0
- MU_MODE1 = 0.01, MU_MODE5 = 0.02, N_MODULES_5 = 8
- EPS_NOISE_5 = 0.05, EPS_NOISE_6 = 0.05
- 10 replicates per scenario, 4 scenarios, 1000 generations each

## Mode 6 implementation note

Implementation A (binary copyability) is the primary test. Mode 6 has a single fixed template per replicate (the 2D-surface pattern), and every Mode 6 agent draws its phenotype each generation as `template + i.i.d. noise`. There is no copy mechanism, so selection within Mode 6 cannot move the population mean phenotype away from the fixed template. Implementations B (per-lineage fixed template) and C (heritability < 1) from dispatch3.md §Step 4 are deferred to a sensitivity follow-up.

## Interpretation -- where the framework's claim succeeds and where it fails

Three of the four pre-registered predictions held cleanly. The framework's qualitative claim (Mode 6 cannot accumulate inherited adaptation; modes that can will out-compete it) is supported by every single replicate of scenarios A, C, and D, and by 9 of 10 in scenario B. Once Mode 1 (or Mode 5) gets any foothold, it climbs and Mode 6 is excluded. The crossover happens fast (median gen 3-5) because Mode 1's mean fitness with a fresh random initial pool is comparable to Mode 6's mean fitness around its fixed template, and Mode 1 starts climbing immediately.

P_H3 was falsified at the strict "all 10 reps" level by a single replicate (scenario B, rep 6, seed 48) where Mode 1 went stochastically extinct from its starting count of 40 before its mutation/selection feedback could produce any high-fitness mutants. In that replicate Mode 6's randomly-initialized fixed template happened to score mean fitness 0.34 against the target, while Mode 1's initial random pool of 40 agents scored mean 0.26 with max ~0.40. The selection pressure in favor of Mode 6's bulk fitness pulled Mode 1 from 40 → 16 → 10 → 5 → 3 → 1 → 0 over six generations, after which there is no Mode 1 left to evolve. In the remaining 9 reps Mode 1 reached >95% by gen 12 at the latest.

What this means for the framework's discussion claim:

- The framework's *qualitative* claim — that copyable modes win when present — survives in every single replicate where a copyable mode persists past the founder bottleneck. No replicate in any scenario showed Mode 6 *evolving* high fitness; whenever Mode 6 won, it won by stochastic drift, not by inheritance.

- The framework's *strict* claim that initial abundance does not matter is falsified by the small-founder regime: when the copyable population is small enough (40 / 400 = 10% here), demographic stochasticity can drive it to extinction before its adaptive advantage manifests, and the non-evolving mode then occupies the niche by default.

- This is a meaningful refinement: copyability is necessary but not sufficient at low initial abundance. The framework's published claim ("an origin-of-life simulation with multiple modes equally available should not develop Mode 6-based inheritance") is true in the sense that Mode 6 never developed inheritance in any replicate — it simply persisted at its initial fitness level when Mode 1 went extinct. The published wording "should not develop Mode 6-based inheritance" is technically defended by every replicate. The stronger pre-registered version (P_H3, "Mode 1 reaches >95% even from 10% initial abundance") is what failed, and the failure mode is mundane population genetics (Wright-Fisher founder loss), not a flaw in the framework's structural argument.

## Suggested follow-ups (not run in v1)

- Implementations B and C of Mode 6 (per dispatch3.md §Step 4) to test whether the framework's claim is binary or graded with respect to copyability strength.
- Scenarios B' / B'' with Mode 1 initial counts of 80, 120 to map the founder-loss boundary.
- Re-run scenario B with N_REPLICATES = 30 to estimate the founder-loss probability with tighter confidence bounds.
