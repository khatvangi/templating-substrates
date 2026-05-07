# Templating Substrates Fig 6: Mode 6 reframing and ciliate cortex bounded-heredity retrodiction

## Purpose

Fig 6 supports the manuscript's claim that the eukaryotic ciliate cortex is
the framework's empirical instance of Mode 6: a 2D pattern that satisfies R1
and R2 at the surface level but fails R3 and R4, and is therefore parasitic
on Mode 1 (the genome) for inheritance. Three panels: (A) schematic of
Beisson-Sonneborn 1965 cortical inversion; (B) the M2 mechanism analog (per-
generation re-draw of cortical units from genome-encoded parts; existing
pattern biases assembly); (C) predicted M2 plateau height matches the
literature-reported cortical inheritance plateau.

Reader takeaway: Mode 6 is real biology, structurally bounded, and the
framework retrodicts the bound from Corollary 2 applied to a substrate
with R3 and R4 failures.

NOTE: Panel A (the Beisson-Sonneborn schematic) is hand-drawn in Inkscape
by the user, not generated on Boron. Boron only handles panels B and C.

## Inputs (for panels B and C only)

- /storage/kiran-stuff/templating_framework/results/test_h5_r_sweep_v1.csv
  (for the M2 r-sweep)
- /storage/kiran-stuff/templating_framework/results/test_h3_isolated_dynamics_v1.csv
  (for plateau heights)

Literature reference for panel C target value:
- Beisson & Sonneborn 1965 PNAS 53:275 — cortical inheritance propagated
  across hundreds of generations
- Frankel 1989 — cortical inheritance bounded; no novel structures accumulate

## Panel specification

### Panel A (USER, NOT BORON): Cortical inheritance schematic

- Panel type: hand-drawn schematic in Inkscape
- Content (3 sub-illustrations stacked or side-by-side):
  1. Wild-type Paramecium with rows of cilia oriented in one direction
  2. Same cell after grafting an inverted patch — cilia in the patch face
     the opposite direction
  3. Daughter cells after several divisions — the inverted patch propagates
     across generations without genomic change
- Style: simple line art, no shading; cilia as small directional triangles
- Reference: Beisson & Sonneborn 1965 PNAS 53:275; Frankel 1989 figure
- Boron does not generate this; user creates it directly in Inkscape

### Panel B: M2 mechanism analog of cortical inheritance

- Panel type: schematic + small data inset
- Schematic (top): show one cell division as
  "existing cortical pattern → new cortex assembled from genome-encoded
  parts, with existing pattern biasing assembly → daughter inherits pattern"
- Data inset (bottom): bar showing where on the M2 r-sweep the cortex sits
  (r ≈ 0.10 sweet spot, with annotation "cortex falls in M2 regime: low
  re-draw rate, lineage-fixation-dominant")
- Library: seaborn for the data inset; user adds the schematic in Inkscape
  on top
- Dimensions: 85 mm × 80 mm (taller than other panels to accommodate both
  schematic and inset)

### Panel C: Predicted vs observed plateau

- Panel type: bar comparison plot
- X-axis: two categories: "Framework prediction (M2, r=0.10)" and "Literature
  observation (cortical inheritance)"
- Y-axis: plateau-like quantity (fitness in simulation; relative
  generational stability in literature)
- Annotation: explicit caveat in caption that the literature value is
  qualitative (cortical patterns persist for hundreds of generations
  without accumulating novel structures), not a quantitative fitness;
  the comparison is therefore directional, not numerical match
- Palette: Okabe-Ito two-color
- Library: seaborn
- Dimensions: 85 mm × 60 mm

## Style requirements

Same as Fig 1.

## Script and output locations

- Scripts: /storage/kiran-stuff/templating_framework/code/figures_v7/
  - fig6_panel_b_m2_inset.py
  - fig6_panel_c_predicted_vs_observed.py
- Outputs: /storage/kiran-stuff/templating_framework/paper/figures/v7/fig6/

## Validation

For panels B and C: same as Fig 1. Panel A is user-drawn; no Boron-side
validation.

## Assembly notes

3 panels in a 1×3 row. Panel A is the largest (cortex schematic with three
sub-cells); panels B and C are equal width. Total ~180 mm. User assembles
all three in Inkscape, with Panel A being directly drawn by the user.

Caption should explicitly state:
1. Panel A is a redrawn schematic after Beisson-Sonneborn 1965, not new data
2. Panel B's "M2 regime" claim is a structural mapping argument, not a
   quantitative fit
3. Panel C's comparison is directional (the framework predicts bounded
   heredity; the literature reports bounded heredity), not a numerical match
