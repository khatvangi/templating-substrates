# v3 Methods note — apparatus channel-ensemble robustness (drop-in)

(Insert this paragraph in v3 Methods immediately after the
$I_\text{struct}^\text{chan}$ definition.)

The numerical value of $I_\text{struct}^\text{chan}$ depends on
the channel ensemble $\mathcal{C}$ used to compute the mixture
distribution $\bar P(\phi_Y(Y))$. Headline numbers in this paper use
the canonical 5-channel ensemble {Drt3a-WC, Drt3b N=2, Drt3b N=3,
Drt3b N=5, AbiK-uniform}. To verify that mode-classification claims
are not artefacts of this choice, we repeated the analysis under
2-channel through 10-channel ensembles drawn from a 10-channel pool
that adds Drt3b N=8, two Mode 1 random-template variants, and two
Mode 5 NRPS-like channels (Test G3,
`results/test_g3_*_v1.csv`). The pairwise KL matrix $D_\text{KL}(P_i \| P_j)$
between any two channels is, by definition, independent of the
broader ensemble; we verified this empirically (Step 1, 10×10 KL
matrix in `test_g3_cross_channel_kl_matrix_v1.csv`). The focal-channel
contribution $D_\text{KL}(P_\text{focal} \| \bar P)$ does
shift with ensemble size (Step 2), but the 10/10 channels
in the pool retain a stable Mode 3 binary classification (criterion:
$I_\text{struct}^\text{chan} > 0.5$ bits AND periodicity peak in
$[0.95, 1.0]$) across 10 randomly drawn 4-channel ensembles
(Step 3). The headline G2 separation between the
fixed-ACACAC Drt3a probe and Drt3b-N=2 (the
$\Delta I_\text{struct}^\text{chan} \approx 1.0$ bits result of
G2 §Step 3) ranges from 0.976 to 1.025 bits across
the 7 alternative ensemble choices (2-channel through 10-channel; Step 4)
and is sign-consistent —
i.e., the fixed-template Drt3a always classifies as more distinct
from the mixture than Drt3b-N=2, regardless of the broader ensemble. We therefore report
$I_\text{struct}^\text{chan}$ values against the canonical
5-channel ensemble while noting that the qualitative mode-classification
verdicts (which channel is Mode 1, which is Mode 3, which is non-
templating) are robust to the ensemble choice within this pool.
