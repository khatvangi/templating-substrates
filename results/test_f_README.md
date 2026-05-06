# Test F -- family-level Mode 3 prediction sweep (v1)

## What this tests

Generalizes Test E v2 (one biological anchor: Drt3b WT + E26Q + AbiK) to ~K
natural clade variants present in the 1,232-sequence Deng et al. 2026 alignment
of Drt3b homologs. Each clade is defined by its 10-residue signature at the
catalytic positions (E26, R168, Y170, G248, R253, Y289, T335, T338, R408, Y650).

## Parameterization rule (the prediction)

The framework's prediction is encoded entirely in the parameterization rule
mapping clade signature -> per-state channel. See module docstring of
`code/test_f_family_sweep.py` for the full rule. Headlines:

  - E26 = E         -> state_A intact (P(G)=0.0033, eps_A=0.01)
  - E26 = D         -> state_A half-broken (P(G)=0.10, eps_A=0.15)
  - E26 = Q / other -> state_A broken    (P(G)=0.20, eps_A=0.20)
  - G248 != G       -> state_C broken (dG misincorp)
  - Y289 != Y       -> state_C broken (pyrimidine recognition lost)
  - R253 != R       -> cycle disrupted; predict I_struct < 0.5 bits
  - secondary residues (R168, R408, Y170, T335, T338) -> NO parameter change
    (the framework predicts these do not shift observables)

## Outputs

  - test_f_family_predictions_v1.csv      one row per clade, with the
                                          framework's predicted observables
  - test_f_per_clade_simulations_v1.csv   raw L-sweep simulation outputs
  - test_f_falsifiability_v1.md           pre-registered numerical thresholds
  - figures/test_f_clade_predictions.png  scatter of predictions in observable space

## How to falsify

For any clade row in test_f_family_predictions_v1.csv, express the named
representative and run cDIP-seq. If a secondary-residue-only clade shifts
its observables outside the in-envelope thresholds (G > 0.05, period-2 < 0.93,
I_struct < 0.95), the framework's residue partitioning is wrong.

## Run summary

  - 1232 sequences in alignment
  - 42 clades found with n >= 5 members
  - sim time: 503.5s
  - L sweep: 4, 8, 16, 32, 64, 128, 256, 500
  - 30 reps x 5000 samples per (clade, L)

## Pass criterion

PASS = secondary-only clades all stay within the WT-like envelope at L=64
AND primary-gate-degraded clades all shift in the predicted direction.

note: Y650 is partitioned as a primary residue but the parameterization rule
documents that Y650 controls initiation, not the per-state elongation channel
we measure. so Y650-only-substituted clades (Y650F) are predicted to stay in
the elongation envelope; the pass criterion treats them accordingly.

See test_f_falsifiability_v1.md for explicit numerical thresholds.

## Run verdict (2026-05-05)

  OVERALL: PASS

  Secondary-only clades:                       33 total, 0 envelope violations
  Per-state primary clades (E26/Y289/G248):     6 total, 6 shifted in predicted direction
                                                (all six are E26->D substitutions:
                                                 marginal G ~0.052, period-2 peak ~0.857,
                                                 I_struct ~0.996 -- all consistent with
                                                 the parameterization rule's E26D entry)
  Initiation-only primary clades (Y650F):       3 total, 0 violations
                                                (all stay in elongation envelope as
                                                 framework predicts)
  Universal-gate clades:                        0 (R253 and G248 are 100%/99.9%
                                                conserved across the family, so the
                                                strongest framework prediction --
                                                I_struct collapse on universal-gate
                                                substitution -- is not testable in
                                                natural variants and would require
                                                site-directed mutagenesis)

