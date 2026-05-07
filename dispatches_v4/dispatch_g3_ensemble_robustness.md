# Dispatch G3: Apparatus channel-ensemble robustness (Test G3)

## Working directory
`/storage/kiran-stuff/templating_framework/`

## Goal
Test G2 defined I_struct^chan against a 5-channel ensemble {Drt3a-WC,
Drt3b-N=2, Drt3b-N=3, Drt3b-N=5, AbiK-uniform}. The numerical value of
I_struct^chan depends on this choice. Reviewers will ask: are the
mode-classification claims robust to ensemble choice, or does the
apparatus produce different verdicts under different ensembles?

G3 sweeps the channel-ensemble choice and tests robustness. This is
defensive but necessary for v3 publication. Cheap.

## Background

Read these first:
1. `code/test_g2_dual_observable.py` — G2 simulation harness; G3
   reuses it.
2. `results/test_g2_v3_methods_paragraph.md` — the apparatus
   definition that requires this robustness check.
3. `results/test_g2_mode_separation_matrix_v1.csv` — the 5×5 cross-
   mode separation matrix from G2 that depends on the ensemble.

## Pre-registered predictions

Document in script header.

P_G3_1 (mode-pair invariance): The cross-mode separation D_KL between
        any two channels is *invariant* to the broader ensemble
        choice. The KL divergence between Drt3a's and Drt3b-N=2's
        product distributions is a property of those two channels
        alone; adding or removing other channels from the ensemble
        does not change it. This is mathematically guaranteed by the
        definition.

P_G3_2 (focal-channel I_struct^chan dependence): The focal-channel
        contribution D_KL(P_focal || P_mixture) DOES depend on the
        ensemble because the mixture distribution depends on which
        channels are included. As ensemble size grows, the mixture
        becomes more "average" and a focal channel that's distinct
        from average gets a higher D_KL.

P_G3_3 (mode-classification stability): Despite ensemble-dependence
        of absolute I_struct^chan values, the *relative* ranking of
        mode classifications (which channel is most distinct, which
        pair is most separable) is invariant across reasonable
        ensemble choices. If a channel is classified as Mode 3 N=2
        under one ensemble, it should be classified as Mode 3 N=2
        under any ensemble that includes the relevant alternatives.

P_G3_4 (Drt3a vs Drt3b sensitivity): The G2 finding that Drt3a-E1
        (single fixed template) separates from Drt3b-N=2 by ΔI_chan
        ≈ 1.0 bits is preserved under ensembles ranging from
        2-channel (just Drt3a + Drt3b-N=2) to 10-channel.

The strongest test is P_G3_3. If mode classification is stable
across ensembles, the apparatus is robust. If it isn't, v3 needs to
specify a canonical ensemble or report ranges.

## Concrete steps

### Step 1: Cross-channel KL matrix is ensemble-invariant

Verify P_G3_1 directly. Compute D_KL(P_i || P_j) for every pair of
channels in the full 10-channel pool:

Pool: {Drt3a-WC, Drt3b-N=2, Drt3b-N=3, Drt3b-N=5, Drt3b-N=8,
Mode 1 random L=32, Mode 1 random L=64, Mode 5 NRPS-like L=8,
Mode 5 NRPS-like L=16, AbiK-uniform}

For each pair, compute D_KL on the empirical product distributions
at L=64 (matching G2). 10 channels × 9 partners = 90 pairwise KLs.

Output: `results/test_g3_cross_channel_kl_matrix_v1.csv`.

This is the bedrock — a 10×10 matrix of pairwise channel separations
that's independent of any specific ensemble.

### Step 2: Ensemble-size sensitivity

Form ensembles of growing size: starting from 2 channels (Drt3a + 1
alternative), grow to 10 channels by adding alternatives one at a
time in a fixed order. For each ensemble:

- Compute I_struct^chan for each focal channel (the focal channel's
  D_KL against the mixture).
- Record how each focal channel's I_chan value changes as the
  ensemble grows.

Specifically, run two ensemble-growth orders:
- Order A: closely-related channels first (Drt3a, Drt3b-N=2,
  Drt3b-N=3, ...)
- Order B: maximally-distant channels first (Drt3a, AbiK, Drt3b-N=8,
  Drt3b-N=2, ...)

Output: `results/test_g3_ensemble_growth_v1.csv` with columns
ensemble_size, ensemble_order, focal_channel, I_chan, mixture_entropy.

### Step 3: Mode-classification stability under ensemble swaps

Define a binary classification rule: a channel is "in mode M3" if its
I_chan signature against the rest of the ensemble exceeds 0.5 bits AND
its periodicity peak is in [0.95, 1.0]. Apply this rule under 10
different ensemble choices (random subsets of size 4 from the 10-
channel pool). For each ensemble, classify Drt3b-N=2 and check
whether the classification flips.

Output: `results/test_g3_classification_stability_v1.csv`.

Predicted: Drt3b-N=2 classifies consistently as Mode 3 across all
ensembles in which it appears. AbiK consistently classifies as
non-Mode-3. Drt3a's classification depends on whether the ensemble
contains a comparison Mode 3 channel.

### Step 4: G2 result reproducibility

Reproduce the G2 mode separation matrix under the canonical 5-channel
ensemble + 4 alternative ensembles (3-channel, 4-channel, 6-channel,
7-channel). Show the Drt3a-vs-Drt3b separation across ensembles.

Output: `results/test_g3_g2_result_robustness_v1.csv`.

Confirm or refine P_G3_4.

### Step 5: v3 Methods note

Generate `results/test_g3_v3_methods_note.md` containing a short
paragraph (3-4 sentences) for v3 Methods specifying:
1. The canonical ensemble used for headline numbers.
2. The reported sensitivity range across plausible alternatives.
3. The mode-classification invariance result (P_G3_3 outcome).
4. Citation to the underlying KL matrix in supplementary data.

This paragraph appears in v3 Methods alongside the G2 apparatus
definition. It addresses reviewer concerns proactively.

## Output artifacts

- `code/test_g3_ensemble_robustness.py`
- `results/test_g3_cross_channel_kl_matrix_v1.csv`
- `results/test_g3_ensemble_growth_v1.csv`
- `results/test_g3_classification_stability_v1.csv`
- `results/test_g3_g2_result_robustness_v1.csv`
- `results/test_g3_v3_methods_note.md`
- `figures/test_g3_kl_matrix_heatmap.png`
- `figures/test_g3_ensemble_growth.png` (I_chan vs ensemble size for
  each focal channel, two growth orders)
- `figures/test_g3_classification_stability.png`

Update `results/CANONICAL_RESULTS.md` to add Test G3.

## Constraints

- Naming: `test_g3` prefix, `_v1` suffix.
- Match G2's simulation harness exactly. n_samples = 5000 per channel
  per ensemble cell.
- Wall-time budget: Step 1: 90 KL computations on existing simulated
  data, near-instant. Step 2: 2 orders × 9 ensemble sizes × ~8 focal
  channels ≈ 144 cells, each computed from existing simulated data,
  fast. Step 3: 10 random ensembles × 4 channels per ensemble = 40
  classifications, fast. Step 4: 5 ensembles × ~5 channels = 25
  cells. Total ~5-10 minutes. Standard dispatch.

## Reporting

Final reply should include:
1. The 10×10 cross-channel KL matrix (or its top-eigenvalue summary).
2. Whether mode classification is stable across ensembles (P_G3_3).
3. The G2 ΔI_chan = 1.0 result robustness across ensembles.
4. The v3 Methods note paragraph.
