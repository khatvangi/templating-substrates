# Test G — Drt3a Mode 1 / Mode 3 boundary resolution (v1)

## Statement of competing predictions

Drt3a (Sharma et al. 2026, *Science*) synthesizes poly(GT) by Watson-Crick
templating off a 6-nt alternating ACACAC region of an ncRNA. Naively the
framework places this in Mode 1 with $L = 6$, alphabet 4, and predicts
$I_\mathrm{struct} \le L \log_2 4 = 12$ bits.

But the template is alternating: only 1 bit of phase information
(ACAC… vs CACA…) is non-redundant. Two predictions:

- **Prediction A (mechanism-determined):** the apparatus measures per-position
  transfer at the joint distribution level. Predicted $I_\mathrm{struct}
  \approx L \log_2 4 = 12$ bits at $L = 6$, regardless of template degeneracy.
  At $\varepsilon = 0.01$ the realistic ceiling is
  $L \cdot (2 - H(Y|X)) = 6 \cdot 1.9034 = 11.4201$ bits.
- **Prediction B (template-information-limited):** the apparatus is bottlenecked
  by the template's intrinsic information content. Predicted
  $I_\mathrm{struct} \le 1$ bit total, because that is all the alternating
  template's population carries (1 bit of phase).

## Empirical result (sweep 1: alternating vs random templates)

At $L_T = 6$, $\varepsilon = 0.01$, $n = 1000$ samples,
per-position plug-in MI estimator (Test A.1 convention):

| template type                          | $I_\mathrm{struct}$ (bits) | bulk-matched | structure-scrambled-Y | verdict |
|----------------------------------------|:--------------------------:|:------------:|:---------------------:|:--------|
| alternating, single fixed ACAC...      | 0.0000 | 0.0000 | 0.0000 | B (template-information-limited, ~1 bit) |
| alternating, 2-phase population        | 5.8438 | 0.0155 | 0.0144 | B (template-information-limited, ~1 bit) |
| random uniform template (A.1 baseline) | 11.4640 | 0.0416 | 0.4051 | A (mechanism-determined, ~12 bits) |

Full sweep (`results/test_g_alternating_vs_random_template_v1.csv`):

| $L_T$ | template          | $I$ (bits) | bulk | scr-Y |
|------:|-------------------|----------:|-----:|------:|
|     2 | alternating_2phase |   1.9131 | 0.0093 | 0.0067 |
|     2 | alternating_fixed |   0.0000 | 0.0000 | 0.0000 |
|     2 | random            |   3.8164 | 0.0247 | 0.9064 |
|     4 | alternating_2phase |   3.8159 | 0.0065 | 0.0058 |
|     4 | alternating_fixed |   0.0000 | 0.0000 | 0.0000 |
|     4 | random            |   7.6795 | 0.0255 | 0.4480 |
|     6 | alternating_2phase |   5.8438 | 0.0155 | 0.0144 |
|     6 | alternating_fixed |   0.0000 | 0.0000 | 0.0000 |
|     6 | random            |  11.4640 | 0.0416 | 0.4051 |
|     8 | alternating_2phase |   7.7231 | 0.0304 | 0.0241 |
|     8 | alternating_fixed |   0.0000 | 0.0000 | 0.0000 |
|     8 | random            |  15.0995 | 0.0602 | 0.3112 |
|    12 | alternating_2phase |  11.4986 | 0.0299 | 0.0351 |
|    12 | alternating_fixed |   0.0000 | 0.0000 | 0.0000 |
|    12 | random            |  22.9036 | 0.0720 | 0.2346 |
|    24 | alternating_2phase |  23.0448 | 0.0485 | 0.0730 |
|    24 | alternating_fixed |   0.0000 | 0.0000 | 0.0000 |
|    24 | random            |  45.8992 | 0.1272 | 0.2804 |

## Empirical result (sweep 2: Drt3a vs Drt3b L-scaling)

Sequence-level joint MI estimator (Test B convention), $X$ = scalar
phase descriptor with 2 states for both modes:

| $L$ | Drt3a (alt-2phase) | Drt3b (Mode 3 N=2) |
|----:|-------------------:|-------------------:|
|   2 |             0.9936 |             0.9810 |
|   4 |             0.9997 |             0.9988 |
|   6 |             0.9993 |             0.9974 |
|   8 |             0.9995 |             0.9998 |
|  12 |             0.9995 |             0.9988 |
|  24 |             0.9991 |             0.9965 |

Linear-fit slopes (informational): Drt3a slope = +1.1565e-04 bits/position,
Drt3b slope = +3.1313e-04 bits/position.

Both curves saturate near $\log_2 2 = 1$ bit. Drt3a's joint MI does NOT
scale linearly with L when X is taken to be the scalar phase index — it is
bounded above by $H(X) = 1$ bit by the data-processing inequality.

## Resolution: which prediction does the apparatus confirm?

The apparatus confirms a **per-position transfer interpretation only when the
template ensemble carries the relevant entropy.** Specifically:

1. Random uniform template (A.1 baseline): the per-position MI estimator
   gives ≈ 11.46 bits at $L_T = 6$, matching Prediction A. This
   confirms the apparatus reads per-position WC transfer correctly when
   the template population is non-degenerate.

2. Alternating-fixed template (single ACAC… for all samples): the per-position
   MI estimator gives ≈ 0.00 bits at $L_T = 6$. This is
   Prediction-B-like. Mechanically, every product position is determined by
   a template position via Watson-Crick — but with X constant across the
   population, the empirical I(X_i; Y_i) is identically zero by definition of
   plug-in MI. The apparatus measures cross-sample variability, not per-sample
   determinism.

3. Alternating-2phase template (X drawn from a 2-element population): the
   per-position MI estimator gives ≈ 5.84 bits at $L_T = 6$.
   This is bounded by H(X) = 1 bit per position (one binary descriptor) and
   by L_T * 1 bit = 6 bits if positions are treated independently — but
   positions are perfectly correlated given the phase, so the joint
   sequence-level MI is also bounded by H(X) = 1 bit (sweep 2 confirms).

**Conclusion: the apparatus confirms Prediction B for the literal Drt3a setup
(alternating template, fixed or 2-phase). Prediction A holds only counterfactually,
when one substitutes the random-uniform-template population as the comparison
ensemble.**

This is a substantive boundary-case finding for the framework. Two readings:

- **Reading 1 (apparatus-as-stated is correct, framing was sloppy):** the
  framework's I_struct is a *population-level* mutual information, not a
  per-sample mechanism descriptor. If the population of templates is
  degenerate, transferable information is by definition limited by template
  entropy. Calling Drt3a 'Mode 1' on the basis of mechanism (per-base WC
  pairing happens) is fine descriptively, but the I_struct numeric must be
  reported relative to the population the descriptor varies over. For Drt3a's
  natural ncRNA-encoded template (one fixed ACACAC), I_struct measured against
  any reasonable comparison ensemble is ~0–1 bit — *not* 12 bits.

- **Reading 2 (apparatus needs a per-sample mechanism observable):** the
  framework should add a complementary observable that captures per-position
  WC transfer at the single-sample level — e.g., the per-position fidelity
  $f_i = P(Y_i = WC(X_i) | X_i)$ averaged over realizations. This *would*
  give ~12 bits worth of 'transfer determinism' at $L_T = 6$ regardless of
  template degeneracy, but it would not be a mutual information.

## Implication for the Mode 1 / Mode 3 boundary

Sweep 2 shows that **with the I_struct apparatus as currently specified,**
Drt3a (Mode 1, alternating template) and Drt3b (Mode 3, N=2) are
*indistinguishable*: both saturate at ~1 bit and neither scales linearly in L.
Drt3a slope = +1.1565e-04 bits/position, Drt3b slope = +3.1313e-04
bits/position. The framework's prose claim ('Drt3a I scales linearly with L,
Drt3b saturates at 1 bit') only holds for the random-template counterfactual,
not for Drt3a as it actually exists.

**The Mode 1 / Mode 3 boundary, as the apparatus measures it, is
output-determined for degenerate-template Mode 1 systems, not
mechanism-determined.** Distinguishing Drt3a from Drt3b therefore requires
either (a) additional observables beyond I_struct (e.g., per-base fidelity,
template separability, ncRNA dependence experiments) or (b) restating the
framework's mode-classification claim to acknowledge that I_struct cannot
distinguish them when both produce the same (AC)_n output and X-population
entropies are matched.

## Provenance

- code: `code/test_g_drt3a_boundary.py`
- sweep 1 csv: `results/test_g_alternating_vs_random_template_v1.csv`
- sweep 2 csv: `results/test_g_drt3a_vs_drt3b_comparison_v1.csv`
- figures: `figures/test_g_alternating_vs_random.png`, `figures/test_g_drt3a_vs_drt3b_scaling.png`
- estimator: per-position plug-in MI on 4×4 joints (Test A.1) for sweep 1; joint sequence-level $H(Y) - H(Y|X)$ (Test B) for sweep 2
- seed: `np.random.seed(42)`, `np.random.default_rng(42)`
- n_samples = 1000, epsilon = 0.01
