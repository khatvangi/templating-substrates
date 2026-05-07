# Templating Substrates Fig 2: Mode-specific information-capacity scaling

## Purpose

Fig 2 supports Corollary 1 (bounded heredity under condition failure) by
showing the framework-predicted scaling laws for each substrate mode. Four
panels: (A) Mode 3 saturates at log_2(N) regardless of L (test B); (B) Mode
1 scales linearly in L while Mode 3 saturates (test B mode contrast); (C)
Mode 3 output autocorrelation peaks at lag = N (test B periodicity); (D) Mode
5 has a length cap at the module count N (test C length limit).

Reader takeaway: each mode's empirical capacity matches the analytical bound
given by which (R1)–(R4) condition that mode fails.

## Inputs

- /storage/kiran-stuff/templating_framework/results/test_b_results.csv
- /storage/kiran-stuff/templating_framework/results/test_c_results.csv
- /storage/kiran-stuff/templating_framework/results/test_c_length_limit.csv

## Existing diagnostic plots to inspect

- figures/test_b_saturation.png
- figures/test_b_capacity.png
- figures/test_b_l_scaling.png
- figures/test_b_mode_contrast.png
- figures/test_b_periodicity.png
- figures/test_c_N_scaling.png
- figures/test_c_length_limit.png

## Panel specification

### Panel A: Mode 3 saturation at log_2 N

- Panel type: line plot
- X-axis: N (cycle states), {2, 3, 4, 5, 6, 8, 10}, log scale
- Y-axis: I_struct^pop (bits)
- Lines: one per L ∈ {N, 5N, 10N, 50N}
- Reference: dashed black line at log_2(N)
- Show that all curves converge to the saturation regardless of L
- Palette class: sequential (different L = different shade); viridis
- Library: seaborn
- Dimensions: 85 mm × 60 mm

### Panel B: Mode 1 linear vs Mode 3 saturating

- Panel type: line plot, two lines on shared axes
- X-axis: L (template length), 16 to 128
- Y-axis: I_struct^pop (bits)
- Line 1: Mode 1 (analytical L · log_2(4) = 2L, plus empirical points from
  test_a1_results.csv)
- Line 2: Mode 3 with N=2 (saturates at log_2 2 = 1 bit, empirical points
  from test_b_results.csv)
- Annotate the gap between the two as "open-ended (R2 satisfied) vs
  bounded (R2 fails)"
- Palette: Okabe-Ito blue (Mode 1), Okabe-Ito orange (Mode 3)
- Library: seaborn
- Dimensions: 85 mm × 60 mm

### Panel C: Mode 3 autocorrelation peak

- Panel type: line plot
- X-axis: lag (positions), 0 to ~10
- Y-axis: autocorrelation coefficient
- Lines: one per N ∈ {2, 3, 4, 6}
- Show peak at lag = N for each curve
- Palette: sequential viridis by N
- Library: seaborn
- Dimensions: 85 mm × 60 mm

### Panel D: Mode 5 module-bounded length

- Panel type: line plot
- X-axis: requested output length L
- Y-axis: actual output length (positions)
- Lines: one per module count N ∈ {3, 5, 10}
- Show the cap at output_length = N regardless of requested L
- Reference: dashed y = x line for the unbounded case
- Palette: sequential cividis by N (works in grayscale)
- Library: seaborn
- Dimensions: 85 mm × 60 mm

## Style requirements

Same as Fig 1: 8 pt axis labels, 7 pt ticks, axis linewidth ≥ 1.0 pt, vector
PDF + 600 dpi PNG per panel, separate files.

## Script and output locations

- Scripts: /storage/kiran-stuff/templating_framework/code/figures_v7/
  - fig2_panel_a_mode3_saturation.py
  - fig2_panel_b_mode1_vs_mode3.py
  - fig2_panel_c_mode3_autocorr.py
  - fig2_panel_d_mode5_length_cap.py
- Outputs: /storage/kiran-stuff/templating_framework/paper/figures/v7/fig2/
  - fig2_panel_a.{pdf,png}
  - fig2_panel_b.{pdf,png}
  - fig2_panel_c.{pdf,png}
  - fig2_panel_d.{pdf,png}

## Validation

Same checklist as Fig 1.

## Assembly notes

2×2 grid in Inkscape, total 180 mm × 130 mm. Panel letters A, B, C, D top-
left of each panel.
