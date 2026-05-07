# Templating Substrates Fig 5: Population-dynamic plateaus from substrate condition failure

## Purpose

Fig 5 supports Corollary 2 (population-dynamic plateaus). Four panels test
the prediction that mechanisms failing R3 or R4 at the individual level are
structurally capped under finite-population dynamics: (A) plateau heights
across M0–M4 at four selection sharpness β values (test H3); (B) M2 sweet-spot
at r ≈ 0.10, regime-invariant across β and L (test H5); (C) M4 vs M2 head-
to-head at matched per-offspring drift rate (test H5 matched comparison); (D)
long-horizon (5000 gen) M4 vs M1 stability of plateau ranking (test H4).

Reader takeaway: the structural cap predicted by Corollary 2 holds across
parameter regimes; it is not an artifact of any specific selection sharpness,
sequence length, or simulation horizon.

## Inputs

- /storage/kiran-stuff/templating_framework/results/test_h3_isolated_dynamics_v1.csv
- /storage/kiran-stuff/templating_framework/results/test_h3_selection_regime_v1.csv
- /storage/kiran-stuff/templating_framework/results/test_h5_r_sweep_v1.csv
- /storage/kiran-stuff/templating_framework/results/test_h5_beta_sensitivity_v1.csv
- /storage/kiran-stuff/templating_framework/results/test_h5_L_sensitivity_v1.csv
- /storage/kiran-stuff/templating_framework/results/test_h5_M2_vs_M4_matched_v1.csv
- /storage/kiran-stuff/templating_framework/results/test_h4_convergence_v1.csv
- /storage/kiran-stuff/templating_framework/results/test_h4_long_horizon_v1.csv

## Existing diagnostic plots

- figures/test_h3_isolated_plateau_heights.png
- figures/test_h3_selection_regime.png
- figures/test_h5_r_sweep_curve.png
- figures/test_h5_regime_sensitivity.png
- figures/test_h5_M2_vs_M4_matched.png
- figures/test_h4_convergence_curves.png

## Panel specification

### Panel A: M0–M4 plateau heights across β

- Panel type: grouped bar chart or line plot
- X-axis: selection sharpness β ∈ {2, 5, 10, 20}
- Y-axis: plateau fitness at gen 1000 (averaged over 30 replicates)
- 5 lines/groups: M0, M1, M2 (at r=0.10), M3, M4
- Reference: dashed line at chance (0.250) and at theoretical max (1.0)
- Match the literal numbers from main paper Table tab:plateau_heights:
  M0=0.249, M1=0.472, M2(r=0.10)=0.533, M3=0.498, M4=0.972 at β=10.
- Palette: Okabe-Ito qualitative 5-color
- Library: seaborn
- Dimensions: 85 mm × 60 mm

### Panel B: M2 sweet-spot regime invariance

- Panel type: line plot, multiple curves
- X-axis: re-draw rate r, log scale, 0.001 to 1.0
- Y-axis: M2 plateau fitness at gen 1000
- One curve per (β, L) combination spanning β ∈ {2, 5, 10, 20} and
  L ∈ {16, 32, 64, 128} — pick a representative subset (e.g., 4–6 curves)
- All curves should peak around r ≈ 0.10 with peak position invariant within
  ±0.05 across regimes
- Annotation: vertical dashed line at r* = 0.10
- Palette: sequential viridis (one shade per regime)
- Library: seaborn
- Dimensions: 85 mm × 60 mm

### Panel C: M4 vs M2 head-to-head at matched drift rate

- Panel type: paired bar chart or scatter
- Comparison: M2 at r* = 0.10 (introduces 2.40 new positions/offspring) vs
  M4 at μ* = 0.001 (introduces 0.02 new positions/offspring)
- Two metrics shown side by side:
  1. Plateau height: M2 = 0.548 vs M4 = 0.997
  2. Direct competition outcome: M4 wins 67%, M2 wins 33%
- Annotation: state "M2 drift rate is 100× M4's, but M4 still wins"
- Palette: M4 = Okabe-Ito blue, M2 = Okabe-Ito orange
- Library: seaborn
- Dimensions: 85 mm × 60 mm

### Panel D: Long-horizon M4 vs M1 stability

- Panel type: line plot
- X-axis: generation, 0 to 5000
- Y-axis: mean fitness across 30 replicates
- Two curves: M4 (rises monotonically to ~0.99) and M1 (plateaus at ~0.47)
- Shaded bands for ±1 SEM across replicates
- Annotation: "M4 wins outright in 28/30 replicates"
- Palette: M4 = Okabe-Ito blue, M1 = Okabe-Ito orange
- Library: seaborn
- Dimensions: 85 mm × 60 mm

## Style requirements

Same as Fig 1.

## Script and output locations

- Scripts: /storage/kiran-stuff/templating_framework/code/figures_v7/
  - fig5_panel_a_m0_m4_plateaus.py
  - fig5_panel_b_m2_sweet_spot.py
  - fig5_panel_c_m4_vs_m2_matched.py
  - fig5_panel_d_long_horizon.py
- Outputs: /storage/kiran-stuff/templating_framework/paper/figures/v7/fig5/

## Validation

Same as Fig 1.

## Assembly notes

2×2 grid in Inkscape, 180 mm × 130 mm total. Panel letters A, B, C, D.
