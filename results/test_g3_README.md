# Test G3 — apparatus channel-ensemble robustness

## what
G2 defined I_struct^chan against a 5-channel ensemble {Drt3a-WC,
Drt3b-N=2, Drt3b-N=3, Drt3b-N=5, AbiK-uniform}. The numerical value of
I_struct^chan depends on this choice. G3 sweeps the channel-ensemble
choice over a 10-channel pool and tests whether the mode-classification
claims are robust.

## why
Reviewers will ask whether the apparatus produces different verdicts
under different ensembles. P_G3_3 (mode-classification stability) is
the strongest test: if classifications flip across ensembles, v3 needs
a canonical ensemble or to report ranges.

## predictions
- P_G3_1: pairwise KL is ensemble-invariant (mathematically guaranteed;
          verified empirically as the 10×10 matrix).
- P_G3_2: focal-channel I_chan DOES depend on ensemble (mixture changes).
- P_G3_3: mode classifications stable across ensembles. STRONGEST.
- P_G3_4: G2's ΔI_chan(Drt3a vs Drt3b-N=2) ≈ 1.0 bits is preserved
          across ensembles ranging from 2-channel to 10-channel.

## PASS criteria
- P_G3_1: 10×10 KL matrix has zero diagonal and finite off-diagonal —
  PASS by construction (any non-degeneracy of empirical sample). Step 1
  populates `test_g3_cross_channel_kl_matrix_v1.csv`.
- P_G3_3: PASS if all/most channels in the pool retain the same Mode 3
  binary verdict across all ensembles in which they appear.
  Result: 10/10 channels stable.
- P_G3_4: PASS if ΔI_chan(Drt3a vs Drt3b-N=2) is sign-consistent across
  all 7 ensembles tested.
  Result: sign-consistent = True; range
  0.976 to 1.025 bits.

## artifacts
- `test_g3_cross_channel_kl_matrix_v1.csv` (Step 1; 10×10 = 100 rows)
- `test_g3_ensemble_growth_v1.csv` (Step 2; ~108 rows)
- `test_g3_classification_stability_v1.csv` (Step 3; 40 rows)
- `test_g3_g2_result_robustness_v1.csv` (Step 4; 7 rows)
- `test_g3_v3_methods_note.md` (Step 5)
- `figures/test_g3_kl_matrix_heatmap.png`
- `figures/test_g3_ensemble_growth.png`
- `figures/test_g3_classification_stability.png`

## deviations from dispatch
- L: dispatch states "L=64 (matching G2)" but G2 actually used L=6 for
  the row-distribution comparison (see
  `results/test_g2_mode_separation_matrix_v1.csv`). At L=64 with
  n_samples=5000 the empirical row distributions for Mode 1 random,
  Mode 5 NRPS, and AbiK become trivially singleton-supported in
  4^64 ≈ 3e38 outcomes. We use L=6 to match G2's actual harness.
- Mode 5 NRPS labels: dispatch says "Mode 5 NRPS-like L=8" and "L=16"
  but the comparison space is L=6. We re-purpose those labels to
  carry distinct N_modules values (L8 → N_modules=6, L16 → N_modules=3)
  so the two Mode-5 channels have distinct module-count signatures
  while sharing the L=6 output length.
- Drt3a-WC: matches G2's alt_2phase parameterization (2-element
  template population), not single-fixed-template, so I_chan is non-
  trivial and the comparison to G2 is direct.
