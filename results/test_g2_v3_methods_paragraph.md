# v3 Methods — apparatus paragraph (drop-in)

(this paragraph replaces the single-observable apparatus definition in v2
Methods/Results lines 35–46 of `templating_substrates_draft_v2.md`.)

## Apparatus definition

We define a templating event as the tuple
$(X, O, S, \Delta G; \phi_X, \phi_Y) \to Y$, with notation as in v2. The framework's
apparatus uses two complementary mutual-information observables, evaluated
descriptor-relatively against the bulk-matched and structure-scrambled controls of v2.

**Observable 1: population mutual information.** With $X$ taken as the *realization*
of the template (a specific 1D sequence, conformer, or module sequence drawn from a
template population $\mathcal{X}$ with distribution $\Pi(X)$), the population MI is

$$ I_\text{struct}^\text{pop}(\mathcal{X}) \;=\; H(\phi_Y(Y)) \;-\; \mathbb{E}_{X \sim \Pi}\left[H(\phi_Y(Y) \mid X)\right]. $$

$I_\text{struct}^\text{pop}$ is bounded above by the data-processing inequality at
$H(\Pi)$, the entropy of the template population. It is identically zero for a fixed
template ($|\mathcal{X}| = 1$, e.g., a single ncRNA sequence in a cell) regardless of
the per-sample mechanism. $I_\text{struct}^\text{pop}$ measures cross-realization
transferable information; it is the *informational capacity of the system as an
encoder*.

**Observable 2: channel-as-X mutual information.** With $X$ taken as the *identity
of the channel itself* (the parameterized templating apparatus, including substrate-
alphabet, mechanism class, and channel parameters), the channel-as-X MI against a
comparison ensemble $\mathcal{C} = \{c_1, \ldots, c_K\}$ of channels is

$$ I_\text{struct}^\text{chan}(\mathcal{C}) \;=\; \sum_{c \in \mathcal{C}} P(C=c)\, D_\text{KL}\!\left(P(\phi_Y(Y) \mid C=c) \,\|\, \bar P(\phi_Y(Y))\right), $$

where $\bar P(\phi_Y(Y)) = \sum_c P(C=c) P(\phi_Y(Y) \mid C=c)$ is the mixture row
distribution. The focal-channel summand,
$D_\text{KL}(P(\phi_Y(Y) \mid C=c_\text{focal}) \,\|\, \bar P(\phi_Y(Y)))$, is
the focal channel's contribution and is the natural per-channel signature.
$I_\text{struct}^\text{chan}$ is bounded above by $\log_2 K$. It is non-zero
whenever the focal channel's product distribution differs from the ensemble
mixture, regardless of whether the focal channel's template population is degenerate.
$I_\text{struct}^\text{chan}$ measures *how distinguishable the channel is from
its alternatives at the product-distribution level*.

**Conditions of validity.** $I_\text{struct}^\text{pop}$ requires a template
population with $|\mathcal{X}| > 1$ to be informative; for a fixed template it is
zero by construction. $I_\text{struct}^\text{chan}$ requires a stated comparison
ensemble $\mathcal{C}$ that includes the focal channel and at least one alternative;
its numerical value depends on the ensemble choice and that choice must be reported
alongside the value.

**Joint signature.** The refined apparatus output is the triple
$(I_\text{struct}^\text{pop},\; I_\text{struct}^\text{chan},\; \bar f)$, where
$\bar f = \mathbb{E}[\mathbf{1}\{Y_i = \text{WC}(X_i)\}]$ is the per-base mechanism
fidelity averaged over positions and realizations (or the appropriate channel-
specific per-position correctness probability for non-WC channels). $\bar f$ is not
a mutual information but is a non-negotiable mechanism observable: it distinguishes
Drt3a's per-base WC pairing (≈ 0.99) from AbiK's chance-level output (= 0.25)
regardless of template degeneracy or channel-ensemble choice.

**Mode classification with the joint signature.**
- Mode 1 (Drt3a): $\bar f \approx 1 - \varepsilon$ regardless of template ensemble;
  $I_\text{struct}^\text{pop}$ scales linearly with $L$ when $\mathcal{X}$ is
  non-degenerate (E3 case) and is zero/bounded when $\mathcal{X}$ is degenerate
  (E1/E2 cases); $I_\text{struct}^\text{chan}$ depends on the ensemble.
- Mode 3 (Drt3b N=2): $\bar f \approx 1 - \varepsilon$ as a *channel*-fidelity;
  $I_\text{struct}^\text{pop}$ saturates at $\log_2 N$ regardless of $L$;
  $I_\text{struct}^\text{chan}$ separates from Mode 1 when the channel ensemble
  contains both.
- Random / non-templating (AbiK): $\bar f = 1/|\mathcal{A}|$;
  $I_\text{struct}^\text{pop} = 0$; $I_\text{struct}^\text{chan}$ is the
  divergence of the uniform distribution from the ensemble mixture.

**Controls.** The bulk-matched and structure-scrambled controls of v2 apply unchanged
to $I_\text{struct}^\text{pop}$. For $I_\text{struct}^\text{chan}$ the analogous
control is $D_\text{KL}(P(\phi_Y(Y_\text{scrambled})) \,\|\, \bar P(\phi_Y(Y)))$,
which probes whether the channel's signature against the ensemble survives positional
scrambling of the focal channel's outputs.
