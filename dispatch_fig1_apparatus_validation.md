# Templating Substrates Fig 1: Apparatus validation and bulk-matched control

## Purpose

Fig 1 supports the manuscript's claim that the diagnostic apparatus
(I_struct^pop, I_struct^chan, f_phase, f_monomer) is a sound, calibrated
instrument. It does this in three panels: (A) shows that the empirical
I_struct^pop tracks its analytical prediction across an error-rate × length
sweep (test A1); (B) shows that the apparatus separates real templating from
bulk-matched controls by orders of magnitude (test A2); (C) shows that
classification verdicts are robust across alternative comparison ensembles
(test G3).

Reader takeaway: the apparatus is calibrated against analytics, controlled
against compositional confounds, and stable across reasonable analyst
choices.

## Inputs

- /storage/kiran-stuff/templating_framework/results/test_a1_results.csv
- /storage/kiran-stuff/templating_framework/results/test_a2_results.csv
- /storage/kiran-stuff/templating_framework/results/test_g3_classification_stability_v1.csv
- /storage/kiran-stuff/templating_framework/results/test_g3_ensemble_growth_v1.csv

## Existing diagnostic plots to inspect (do not embed verbatim — re-render)

- figures/test_a1_scaling.png
- figures/test_a2_comparison.png
- figures/test_g3_classification_stability.png
- figures/test_g3_ensemble_growth.png

## Panel specification

### Panel A: Empirical vs analytical I_struct^pop

- Panel type: line plot with error bars or scatter with diagonal y=x reference
- X-axis: analytical prediction (bits)
- Y-axis: empirical estimate (bits)
- Group/color by: per-position misincorporation rate ε ∈ {0.001, 0.01, 0.1}
- Marker shape by: L ∈ {16, 32, 64, 128}
- Reference line: y = x (1:1)
- Annotation: report Pearson r and slope from a linear fit on the points
- Palette class: sequential (one ε value = one color shade); use viridis or
  cividis
- Library: seaborn or plotly-static
- Dimensions: 60 mm × 60 mm (panel A is one of three side-by-side at
  180 mm total figure width)

### Panel B: Bulk-matched control discrimination

- Panel type: bar chart or paired strip plot
- X-axis: condition (3 template biases: π_uniform, π_AT-skew, π_GC-skew)
- Y-axis: I_struct^pop (bits, log scale)
- Two bars/strips per x position: "genuine templating" and "bulk-matched
  control"
- Annotation: separation ratio (genuine / bulk) above each pair, in text
- Pass criteria from test A2 README: bulk < 0.05 bits; separation ratio > 20×
- Palette class: qualitative two-color (one for genuine, one for control);
  use Okabe-Ito blue and orange or similar
- Library: seaborn
- Dimensions: 60 mm × 60 mm

### Panel C: Comparison-ensemble robustness

- Panel type: violin or box plot
- X-axis: ensemble identifier (1 = canonical, 2-11 = alternative draws)
- Y-axis: I_struct^chan (bits)
- One distribution per ensemble across the test_g3 classification points
- Reference line: I_struct^chan from the canonical 5-channel ensemble
- Annotation: state in caption that the qualitative classification (Drt3a as
  Mode 1, Drt3b as Mode 3) is preserved across all 11 ensembles
- Palette: single neutral color for the alternative ensembles, accent color
  for canonical
- Library: seaborn
- Dimensions: 60 mm × 60 mm

## Style requirements

- Use a project style module if available; otherwise set seaborn defaults to
  context="paper", style="ticks", palette as specified per panel.
- Font: sans-serif, size 8 pt for axis labels, 7 pt for tick labels, 9 pt
  for panel letters (A, B, C — added in Inkscape, not in Python).
- Axis linewidth ≥ 1.0 pt.
- Output: each panel as a SEPARATE PDF (vector) AND PNG (600 dpi) file. No
  Python-side composite.
- All text must be selectable in the PDF (no flattened raster text).

## Script and output locations

- Scripts: /storage/kiran-stuff/templating_framework/code/figures_v7/
  - fig1_panel_a_apparatus_calibration.py
  - fig1_panel_b_bulk_matched.py
  - fig1_panel_c_ensemble_robustness.py
- Outputs: /storage/kiran-stuff/templating_framework/paper/figures/v7/fig1/
  - fig1_panel_a.pdf, fig1_panel_a.png
  - fig1_panel_b.pdf, fig1_panel_b.png
  - fig1_panel_c.pdf, fig1_panel_c.png

## Validation per panel

Each script must verify before writing output:

1. Output file path exists after save
2. PDF opens (pdfinfo or PyPDF2)
3. PNG opens (Pillow)
4. All axis labels are non-empty
5. No text-overlap warnings
6. Minimum font size ≥ 7 pt
7. File size < 5 MB

If any check fails, report the failure in the script's stderr log and exit
with non-zero status.

## Assembly notes

User assembles in Inkscape:
- 3 panels side-by-side, 180 mm total width
- Add panel letters A, B, C in 10 pt bold sans-serif at top-left of each panel
- Add a single shared figure caption below
- Final output: fig1.pdf (composite) at the manuscript's figures directory
