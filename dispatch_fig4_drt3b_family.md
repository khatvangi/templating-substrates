# Templating Substrates Fig 4: Drt3b family-level apparatus application

## Purpose

Fig 4 shows that the Mode 3 classification is robust across the Drt3b family
(1,232 homologs from Deng et al. 2026 supplementary alignment), with two
specific predictions: (A) 33 secondary-residue clades fall in-envelope at the
WT-Drt3b apparatus signature, confirming that secondary residue diversity does
not break the cyclic-conformational mode; (B) 6 E26→D primary-gate clades
shift in framework-predicted direction (f_monomer^A degraded to 0.71–0.84
while f_phase preserved at 0.94–0.98), consistent with selectivity-gate
disruption rather than architectural collapse.

Reader takeaway: the apparatus is family-portable; predictions hold across
1,232 sequences, not just at the reference.

## Inputs

- /storage/kiran-stuff/templating_framework/results/test_f_family_predictions_v1.csv
- /storage/kiran-stuff/templating_framework/results/test_f_per_clade_simulations_v1.csv
- /storage/kiran-stuff/templating_framework/results/test_f2_secondary_residue_epsilon_v1.csv
- /storage/kiran-stuff/templating_framework/results/test_f2_E26D_sensitivity_v1.csv

## Existing diagnostic plots

- figures/test_f_clade_predictions.png
- figures/test_f2_secondary_residue_tolerance.png
- figures/test_f2_E26D_sensitivity.png

## Panel specification

### Panel A: 33 secondary-residue clades in-envelope

- Panel type: scatter plot in 2D apparatus space
- X-axis: f_phase
- Y-axis: f_monomer (averaged across A and C states, or pick one)
- Each point: one secondary-residue clade
- Reference shaded region: WT-Drt3b apparatus envelope (f_phase > 0.95,
  f_monomer > 0.95) drawn as a green-tinted rectangle in upper-right
- All 33 clade points should fall inside the envelope
- Annotation: "33/33 clades in-envelope"
- Palette: single neutral color for clade points (dark gray); green tint
  for envelope; emphasize WT reference point in accent color
- Library: seaborn
- Dimensions: 85 mm × 60 mm

### Panel B: 6 E26→D primary-gate clades, framework-direction shift

- Panel type: paired-points or arrow plot
- X-axis: f_phase
- Y-axis: f_monomer^(A-state) (the one that should degrade with E26→D)
- WT-Drt3b reference at top-right (f_phase ≈ 0.99, f_monomer^A ≈ 0.99) in
  accent color
- 6 E26→D clade points showing degraded f_monomer^A (0.71–0.84) but
  preserved f_phase (0.94–0.98), drawn in vermilion
- Optional: arrows from WT reference to each E26→D point to show direction
- Annotation: "6/6 clades shift in framework-predicted direction"
- Palette: Okabe-Ito blue (WT reference), vermilion (E26→D clades)
- Library: seaborn
- Dimensions: 85 mm × 60 mm

## Style requirements

Same as Fig 1.

## Script and output locations

- Scripts: /storage/kiran-stuff/templating_framework/code/figures_v7/
  - fig4_panel_a_secondary_clades_in_envelope.py
  - fig4_panel_b_e26d_predicted_shift.py
- Outputs: /storage/kiran-stuff/templating_framework/paper/figures/v7/fig4/

## Validation

Same as Fig 1.

## Assembly notes

2 panels side-by-side, 180 mm total width. Panel letters A, B top-left.
Caption: clarify "in-envelope" criterion explicitly (f_phase > 0.95,
f_monomer > 0.95).
