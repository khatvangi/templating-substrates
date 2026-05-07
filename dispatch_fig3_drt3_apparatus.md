# Templating Substrates Fig 3: Apparatus signature on the Drt3 system

## Purpose

Fig 3 is the empirical anchor of the manuscript. It shows that the diagnostic
apparatus correctly classifies four real biological systems with known but
distinct mechanisms: WT Drt3b (Mode 3 cyclic-conformational), E26Q Drt3b
(gate-broken, predicted dG misincorporation), AbiK (same fold, non-templating),
and Drt3a (Mode 1 sequence-template). Three panels: (A) periodicity peaks
distinguishing WT from E26Q; (B) the apparatus triple separating all four
systems in joint-observable space; (C) channel-MI separating Drt3a Mode 1 from
Drt3b Mode 3.

Reader takeaway: the apparatus distinguishes mechanisms that share output
statistics; the joint signature, not any scalar summary, is what classifies.

## Inputs

- /storage/kiran-stuff/templating_framework/results/test_e_v2_results.csv
- /storage/kiran-stuff/templating_framework/results/test_g2_dual_observable_v1.csv
- /storage/kiran-stuff/templating_framework/results/test_g2_drt3a_three_ensembles_v1.csv
- /storage/kiran-stuff/templating_framework/results/test_g2_mode_separation_matrix_v1.csv

Note: the analytical E26Q prediction is marginal G fraction = 0.10; the
observed value in Deng et al. 2026 Fig. 4J is 0.1016. State this match in
the figure caption (not in the panel itself).

## Existing diagnostic plots

- figures/test_e_v2_periodicity.png
- figures/test_e_v2_three_systems.png
- figures/test_g2_dual_observable_signature.png
- figures/test_g2_mode_separation_heatmap.png

## Panel specification

### Panel A: Drt3b WT vs E26Q periodicity peaks

- Panel type: paired bar plot or two-line autocorrelation overlay
- X-axis: lag (0 to 6 positions)
- Y-axis: autocorrelation coefficient
- Two lines: WT Drt3b (peak ~0.98 at lag 2) and E26Q (peak ~0.83 at lag 2)
- Annotate: "WT alternation: 0.99"; "E26Q alternation: 0.83"
- Palette: Okabe-Ito two-color (sky-blue for WT, vermilion for E26Q)
- Dimensions: 60 mm × 60 mm

### Panel B: Apparatus triple separating four systems

- Panel type: scatter or bar group across three observables
- X-axis: three observables: I_struct^pop, periodicity peak, marginal G fraction
- Y-axis: value of observable (different per observable; use grouped bars or
  faceted panels)
- Four colored markers/bars per observable: Drt3b WT, E26Q, AbiK, Drt3a
- Best layout: three side-by-side mini-panels, one per observable, with the
  same four-system color key
- Palette: Okabe-Ito 4-color qualitative
- Library: seaborn or plotly-static
- Dimensions: 60 mm × 60 mm

### Panel C: Channel-MI separation Drt3a Mode 1 vs Drt3b Mode 3

- Panel type: bar chart with annotation
- X-axis: 5 channels in the comparison ensemble (Drt3a Mode 1, Drt3b Mode 3 WT,
  Drt3b E26Q, Drt3a homopolymer null, AbiK random)
- Y-axis: I_struct^chan against the 5-channel ensemble (bits)
- Drt3a and Drt3b WT both at 1.0 bits (saturate against the ensemble)
- Annotation: "Δ I_struct^chan(Drt3a, Drt3b) = 1.0 bits"
- Palette: 5-color qualitative, but emphasize Drt3a and Drt3b WT in saturated
  color and the others in muted
- Dimensions: 60 mm × 60 mm

## Style requirements

Same as Fig 1.

## Script and output locations

- Scripts: /storage/kiran-stuff/templating_framework/code/figures_v7/
  - fig3_panel_a_wt_vs_e26q_periodicity.py
  - fig3_panel_b_apparatus_triple.py
  - fig3_panel_c_channel_mi_separation.py
- Outputs: /storage/kiran-stuff/templating_framework/paper/figures/v7/fig3/

## Validation

Same as Fig 1.

## Assembly notes

3 panels side-by-side, 180 mm total. Panel letters A, B, C top-left.
Caption should include the E26Q analytical-vs-observed match
(0.10 vs 0.1016).
