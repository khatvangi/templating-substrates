# Dispatch F3: Universal-gate architecture vs fidelity characterization (Test F3)

## Working directory
`/storage/kiran-stuff/templating_framework/`

## Goal
Test F2 found that R253A and G248A predictions diverge in an unexpected
way. R253A: I_struct=0.0090 (cycle architecture collapses, Mode 3 fails).
G248A: I_struct=0.9966 (cycle architecture preserved, but marginal G
shifts to 0.1516 and period-2 peak drops to 0.7469). The dispatch had
predicted both would collapse I_struct similarly. The data show R253
is *architectural* (its substitution destroys the cycle) while G248 is
*selectivity-only* (its substitution degrades fidelity within an intact
cycle).

This is a substantive refinement of the framework's universal-gate claim,
not just an unexpected number. v2 treated R253 and G248 as
indistinguishable — both invariant, both universal gates. F3
characterizes the architectural-vs-selectivity distinction with finer
resolution and generates separate testable predictions for each gate
type.

This is Test F3. Its purpose is to give v3 the language to distinguish
architectural gates (whose substitution disrupts mode classification)
from selectivity gates (whose substitution degrades fidelity within
the mode).

## Background

Read these first:
1. `code/test_f2_robustness.py` and `results/test_f2_universal_gate_hypotheticals_v1.csv`
   — the simulation that revealed R253 vs G248 divergence.
2. `results/test_f2_robust_falsifiability_v1.md` §3 — the SDM
   prediction table.
3. `paper/drafts/draft_v2.md` Results section on Drt3b conservation
   pattern (R253 100%, G248 99.9%) — the universal-gate claim being
   refined.
4. `framework_formal_v1.md` Mode 3 mechanism description — the
   architectural claim about cyclic protein conformations encoding
   sequence specificity.

## The architectural-vs-selectivity distinction

Reframe the universal gates:

**Architectural gates**: residues whose substitution disrupts the
cycle's state structure itself. The N-state cycle requires R253 (per
the Deng et al. cation-π interaction with the templating dC residue
of the ncRNA). Without R253, the cycle does not transition cleanly
between A-state and C-state; the system loses its alternation
entirely. Substitution → loss of mode classification. I_struct
collapses because there is no cycle to transfer information.

**Selectivity gates**: residues whose substitution degrades fidelity
within an otherwise intact cycle. G248 enforces dG exclusion at
the C-state by steric occlusion. Without G248, the cycle still
alternates between states but the C-state's dG selectivity is broken.
The system remains Mode 3 N=2 but with degraded ε. I_struct stays
high (cycle present); marginal G shifts (selectivity broken);
period-2 peak modestly degrades (cycle intact but noisier).

This distinction matters because it generates two qualitatively
different SDM predictions. Any structural biologist running R253A
should observe loss of cycle architecture (I_struct < 0.1 bits, no
period-2 peak). Any structural biologist running G248A should observe
preserved cycle but degraded selectivity (I_struct > 0.95 bits,
period-2 peak ~ 0.75, marginal G shifted).

## Pre-registered predictions

Document in script header.

P_F3_1 (R253A architecture loss): I_struct < 0.1 bits at L=64,
        period-2 peak < 0.30, marginal G ≈ 0.125 (chance-level).
        Cycle architecture lost; system reduces to Mode-1-like
        random output without template — i.e., AbiK-like.

P_F3_2 (G248A selectivity degradation): I_struct > 0.95 bits at L=64,
        period-2 peak in [0.70, 0.80], marginal G in [0.13, 0.17].
        Cycle architecture preserved; only the C-state's dG selectivity
        is degraded.

P_F3_3 (other universal-position substitutions): For the other primary
        gates conserved at >90% (E26, Y289, Y650), substitution should
        produce predictions intermediate between architecture-loss
        (R253-like) and selectivity-degradation (G248-like). Specifically:
        - E26 substitutions affect state_A's dA specificity → marginal
          G shift but cycle preserved → "selectivity" type
        - Y289 substitutions affect state_C's pyrimidine recognition
          → fidelity degradation but cycle preserved → "selectivity"
          type
        - Y650 substitutions affect priming, not catalysis → predicted
          minimal change in cycle observables, possible failure to
          initiate at all

P_F3_4 (R253A vs G248A double mutant): The double mutant should look
        R253A-like, not G248A-like, because architecture loss precedes
        selectivity degradation. Predicted: I_struct < 0.1 bits.

P_F3_5 (intermediate R253 substitutions): R253K (conservative,
        positive charge preserved) should partially preserve cycle
        architecture: I_struct intermediate between WT (0.99) and R253A
        (0.01), perhaps 0.5–0.8. R253A is the strongest architectural
        disruption; K is a chemical near-equivalent.

If P_F3_2 fails — if G248A actually breaks the cycle (I_struct < 0.5
bits) — the architectural/selectivity distinction is not stable and
both gates should be classified together. If P_F3_1 fails — if R253A
preserves the cycle (I_struct > 0.9 bits) — the framework's claim
about R253's architectural role is wrong and the cycle's structural
basis needs reidentification.

## Concrete steps

### Step 1: Refined universal-gate sweep

Re-run the F2 universal-gate hypotheticals with finer parameter
resolution. For each gate, sweep substitution severity:

For R253:
- WT (R253-R): state_C fidelity 0.99 (canonical)
- R253K: state_C fidelity 0.85 (conservative substitution; cycle
  partially preserved)
- R253H: state_C fidelity 0.60 (less conservative, intermediate)
- R253A: state_C fidelity 0.25 (uniform random; full architectural
  disruption)

For G248:
- WT (G248-G): C-state dG misincorp 0.005 (canonical exclusion)
- G248A: C-state dG misincorp 0.30 (steric exclusion lost; cycle
  intact)
- G248V: C-state dG misincorp 0.15 (intermediate; partial steric
  preservation)
- G248D: C-state dG misincorp 0.20 (charged residue; geometry shifted
  but cycle present)

For each substitution, simulate apparatus signature at L ∈ {16, 32,
64, 128, 500} with n=1000 reps per cell.

Output: `results/test_f3_gate_substitution_sweep_v1.csv`.

### Step 2: Architecture-vs-selectivity classifier

For each substitution, classify as:
- "architectural disruption": I_struct < 0.5 bits (cycle lost)
- "selectivity degradation": I_struct > 0.95 bits AND marginal G > 0.05
  AND period-2 peak < 0.95 (cycle intact, fidelity degraded)
- "near-WT": I_struct > 0.95 bits AND marginal G < 0.01 AND period-2
  peak > 0.95 (cycle and fidelity both intact, near-WT)
- "ambiguous": none of above

Output: `results/test_f3_classification_v1.csv` with substitution,
classification, three observable values.

This generates the architecture-vs-selectivity decision boundary in
parameter space.

### Step 3: Other primary-gate predictions (E26, Y289, Y650)

For each, simulate the dispatched substitution series:

For E26:
- E26-E (WT): state_A fidelity 0.99
- E26-D: state_A fidelity 0.85 (the F2 prediction)
- E26-Q: state_A fidelity 0.50 (the Test E v2 anchor)
- E26-A: state_A fidelity 0.25

For Y289:
- Y289-Y (WT): state_C pyrimidine recognition 0.99
- Y289-F: state_C pyrimidine recognition 0.90 (conservative; aromatic
  preserved)
- Y289-A: state_C pyrimidine recognition 0.50

For Y650:
- Y650-Y (WT): priming intact, no cycle effect
- Y650-F: priming intact, possible minor effect
- Y650-A: priming lost; predict no product (zero observation count
  rather than degraded apparatus)

Apply the same architecture/selectivity classifier from Step 2.

Output: `results/test_f3_other_primary_gates_v1.csv`.

### Step 4: Double mutants and synergy

For four pairs, simulate the double mutant and compare to single-
mutant predictions:
- R253A + G248A
- R253A + E26Q
- G248A + E26Q
- E26D + G248A

Predicted (from F3 logic): R253A's architectural disruption dominates
all double mutants involving it. Pairs not involving R253 should be
additive in fidelity space.

Output: `results/test_f3_double_mutants_v1.csv`.

### Step 5: SDM-ready prediction table for v3

Generate `results/test_f3_sdm_predictions_v1.md` containing a single
table organized by gate type:

Table 1: Architectural gates (single substitution disrupts cycle)
- Gates: R253 (and any other identified by Step 2 classification)
- For each: WT signature, predicted SDM signature, expected
  classification, falsifiability claim

Table 2: Selectivity gates (single substitution degrades fidelity)
- Gates: G248, E26, Y289, possibly others
- For each: WT signature, predicted SDM signature, expected
  classification, falsifiability claim

Table 3: Priming gates (substitution prevents activity entirely)
- Gates: Y650
- For each: WT signature, predicted "no product" signature

Each row is one SDM experiment a structural biologist could run.

Output: `results/test_f3_sdm_predictions_v1.md`.

### Step 6: v3 Results paragraph draft

Generate `results/test_f3_v3_results_paragraph.md` containing a
drop-in paragraph for v3 Results that:
1. States the architecture/selectivity distinction explicitly
2. Cites the specific predictions for each gate type
3. Notes that the distinction generates orthogonal SDM predictions:
   architectural gates predict cycle loss; selectivity gates predict
   intact cycle with shifted observables
4. Frames this as a strengthening of the framework's universal-gate
   claim, not a retraction

This paragraph augments F2's findings into v3's Results.

## Output artifacts

- `code/test_f3_gate_substitution_sweep.py`
- `results/test_f3_gate_substitution_sweep_v1.csv`
- `results/test_f3_classification_v1.csv`
- `results/test_f3_other_primary_gates_v1.csv`
- `results/test_f3_double_mutants_v1.csv`
- `results/test_f3_sdm_predictions_v1.md`
- `results/test_f3_v3_results_paragraph.md`
- `figures/test_f3_gate_substitution_signatures.png` (joint scatter
  of I_struct vs marginal G for all substitutions, color-coded by
  classification)
- `figures/test_f3_classification_decision_boundary.png`

Update `results/CANONICAL_RESULTS.md` to add Test F3.

## Constraints

- Naming: `test_f3` prefix, `_v1` suffix.
- Do not modify Test F or F2 outputs.
- Match F2's simulation harness; only the substitution parameters
  change.
- Wall-time budget: ~700 simulations across all steps. ~15 minutes if
  cells are fast. Standard (non-nohup) dispatch.

## Reporting

Final reply should include:
1. The architecture-vs-selectivity classification table for all
   simulated substitutions.
2. Whether the distinction holds (P_F3_1 and P_F3_2 confirmed?).
3. The full SDM-ready prediction table from Step 5.
4. Path to v3_results_paragraph.md.
