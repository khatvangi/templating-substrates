# Dispatch F2: Test F robustness — fidelity-parameter sensitivity (Test F2)

## Working directory
`/storage/kiran-stuff/templating_framework/`

## Goal
Test F v1 found 33/33 secondary-residue clades in-envelope and 6/6 (or
9/9 per the falsifiability statement) primary-gate clades shifted in
predicted direction. The pre-registered envelope (G ≤ 0.0033, period-2
≥ 0.9793, I_struct ≥ 0.9996 at L=64) was computed at a single fidelity
parameter set. Test F2 sweeps the fidelity parameters to characterize
the envelope's robustness — what range of ε values still keeps clades
in-envelope, and where does the envelope break?

This is Test F2. It strengthens Test F's falsifiability statement by
turning a point-estimate envelope into a range with explicit ε-dependence.

## Background

Read these first:
1. `code/test_f_family_sweep.py` — Test F v1 implementation.
2. `results/test_f_family_predictions_v1.csv` — the per-clade prediction
   set.
3. `results/test_f_falsifiability_v1.md` — the pre-registered numerical
   envelope at the canonical ε.
4. `dispatch_01_family_sweep.md` (the original dispatch) — the
   parameterization rule.

## What F2 adds

Test F's parameterization rule sets specific ε values (state_A
fidelity 0.99 for E26-intact, 0.50 for E26Q-broken, 0.85 for E26D).
These were derived from Deng et al.'s reported E26Q misincorporation
data (80% dA / 20% dG). The E26D rate is interpolated, not measured.

F2 asks: how sensitive is the in-envelope/out-of-envelope classification
to the parameterization rule? If the rule is approximately right (E26D
gives intermediate fidelity in {0.7, 0.85, 0.95}), do the predictions
hold? If the rule is exactly wrong (E26D gives 0.50 like E26Q), does
the framework's classification still survive?

## Concrete steps

### Step 1: Sensitivity sweep on E26D parameterization

The E26D parameterization is the soft point in Test F v1. Sweep:
- state_A fidelity for E26D ∈ {0.50, 0.60, 0.70, 0.80, 0.85, 0.90, 0.95, 0.99}
- For each value, recompute Test F's predictions for the 9 primary-gate
  clades that include E26D representatives (or whatever subset has E26D
  per `test_f_family_predictions_v1.csv`).
- Record: marginal G, period-2 peak, I_struct, in-envelope/out-of-envelope.

Output: `results/test_f2_E26D_sensitivity_v1.csv`.

The framework's claim: at E26D fidelity ≥ 0.85, the predicted observables
should remain distinguishable from E26Q (more like WT than like Q). This
is what would make the E26D prediction publishable as a falsifiable
claim.

### Step 2: Universal-gate hypothetical sweep (R253, G248)

Test F v1 found 0 natural variants at R253 or G248. The framework predicts
SDM mutants (R253A, G248A) would give I_struct < 0.5 bits and period-2
peak ≈ 0.25 (no cycle structure). F2 explicitly simulates these
hypothetical mutants:

- R253A: state_C fidelity collapses; every position in C-state
  misincorporates with high probability. Parameterize as state_C fidelity
  0.25 (uniform random).
- G248A: state_C dG exclusion fails; dG misincorporation ~ 0.30.
- R253A + G248A double mutant.

For each, simulate the apparatus signature at L ∈ {16, 32, 64, 128, 500}.

Output: `results/test_f2_universal_gate_hypotheticals_v1.csv`.

This generates pre-registered predictions for SDM experiments any structural
biologist could run. They strengthen the framework's claim that universal
gates are non-negotiable.

### Step 3: Secondary-residue ε perturbation

The framework predicts secondary residues do not affect mode classification.
F2 tests: if a secondary-residue substitution somehow shifted ε (say, the
substitution destabilizes the active-site geometry), would the clade still
be in-envelope?

For the 33 in-envelope clades, perturb their effective ε:
- ε_perturbed ∈ {0.01, 0.02, 0.05, 0.10, 0.15}
- For each, recompute predicted observables.
- Find the ε at which the clade exits the envelope.

Output: `results/test_f2_secondary_residue_epsilon_v1.csv`.

This tells us how much ε perturbation a secondary-residue substitution
would need to falsify the framework's claim. If the answer is "ε must
double from 0.01 to 0.02" the framework's claim is fragile; if the answer
is "ε must reach 0.10," the claim is robust.

### Step 4: Replicate variability

Test F v1 used n=100 reps per L. F2 reruns with n=500 reps per L for the
WT clade and the strongest E26Q clade, to get tight CIs on the apparatus
observables.

Output: `results/test_f2_replicate_CIs_v1.csv`.

Confirm Test F's observed numbers fall within these CIs (sanity check).

### Step 5: Robust falsifiability statement

Generate `results/test_f2_robust_falsifiability_v1.md` containing:
1. Pre-registered envelope from F v1, with 95% CI from F2 Step 4.
2. E26D parameterization sensitivity range (Step 1) — what fidelity
   values keep E26D distinguishable from WT and from Q?
3. Universal-gate hypothetical predictions (Step 2) — explicit numerical
   predictions for R253A and G248A SDM mutants.
4. Secondary-residue ε tolerance (Step 3) — how much ε perturbation
   the framework's classification can absorb before failing.

This is the v3-ready falsifiability statement: an envelope with stated
parameter sensitivity, not a point estimate.

## Output artifacts

- `code/test_f2_robustness.py`
- `results/test_f2_E26D_sensitivity_v1.csv`
- `results/test_f2_universal_gate_hypotheticals_v1.csv`
- `results/test_f2_secondary_residue_epsilon_v1.csv`
- `results/test_f2_replicate_CIs_v1.csv`
- `results/test_f2_robust_falsifiability_v1.md`
- `figures/test_f2_E26D_sensitivity.png`
- `figures/test_f2_universal_gate_hypotheticals.png`
- `figures/test_f2_secondary_residue_tolerance.png`

Update `results/CANONICAL_RESULTS.md` to add Test F2 alongside Test F v1.

## Constraints

- Naming: `test_f2` prefix, `_v1` suffix.
- Do not modify Test F v1 outputs.
- Match Test F v1 simulation harness exactly; only the parameter sweeps
  change.
- Wall-time budget: ~5000 simulations across all steps. ~30 minutes
  if cells are fast. If projecting longer, dispatch nohup with progress.txt.

## Reporting

Final reply should include:
1. The robust envelope (with CIs) for Test F's three observables.
2. The E26D fidelity range that keeps the prediction valid.
3. Universal-gate hypothetical predictions (R253A, G248A) as numbers.
4. Path to robust_falsifiability_v1.md.
