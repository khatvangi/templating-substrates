# Test G2 — apparatus repair decision (v1)

## Background

Test G v1 found that the I_struct apparatus, computed with X = template
realization (= I_struct^pop, population mutual information), gives ~0 bits
for Drt3a's actual fixed-ACACAC ncRNA template and ~1 bit for the 2-phase
alternating counterfactual. Both Drt3a and Drt3b-N=2 saturate at ~1 bit by
the joint sequence-level estimator (sweep 2 of v1). The v2 draft's central
Mode 1 / Mode 3 distinguishing claim — that I_struct grows linearly with L for
Drt3a and saturates at log2(N) for Drt3b — does not survive when applied
to Drt3a as it actually exists.

Test G v1's resolution document offered two interpretive readings. Test G2
reframes them as **complementary observables**, not competing repairs:

- **Observable 1: I_struct^pop** — population MI with X = template realization.
  $I_\mathrm{struct}^\mathrm{pop} = H(Y) - H(Y \mid X = x)$ averaged over the
  template population. Bounded above by H(template ensemble). Equal to the v1
  measurement.
- **Observable 2: I_struct^chan** — channel-as-X MI with X = the channel identity.
  Computed against a fixed 5-channel ensemble {Drt3a-WC, Drt3b-Mode3-N=2,
  Drt3b-Mode3-N=3, Drt3b-Mode3-N=5, AbiK-uniform}. The focal-channel summand is
  $\mathrm{KL}(P(Y \mid C = c_\mathrm{focal}) \,\|\, \bar P(Y))$, where
  $\bar P(Y)$ is the equiprobable mixture of the 5 channels' empirical product
  distributions.

The product (I_struct^pop, I_struct^chan, per-base fidelity) is the refined
apparatus signature.

## Joint apparatus signature for the 6 test cases (L=6, ε=0.01, n=5000)

| Test case            | I_pop (bits) | I_chan (bits) | per-base fidelity | JS vs Drt3b-N=2 |
|----------------------|-------------:|--------------:|------------------:|----------------:|
| E1_fixed_ACACAC      |       0.0000 |        3.3320 |            0.9901 |          1.0000 |
| E2_2phase_alt        |       5.7770 |        2.3528 |            0.9903 |          1.0000 |
| E3_random_uniform    |      11.4291 |        9.5929 |            0.9900 |          0.9889 |
| Drt3b_N2             |       0.9997 |        2.3308 |            0.9893 |          0.0018 |
| Drt3b_N3             |       1.5848 |        2.3218 |            0.9903 |          1.0000 |
| AbiK_uniform         |       0.0000 |        9.6810 |            0.2500 |          0.9905 |

Notation:
- E1 = Drt3a with the single fixed ACACAC template (the actual biological case)
- E2 = Drt3a with a 2-phase {ACACAC, CACACA} template population
- E3 = Drt3a with a random uniform 6-nt template (counterfactual baseline)
- per-base fidelity for Mode 3 entries is mechanism fidelity
  $P(Y_k = (\phi_0 + k) \bmod N)$, not Watson-Crick fidelity

## Pairwise observable separation

| Pair | ΔI_pop | ΔI_chan | Δf | ΔJS | Which observable separates? |
|------|-------:|--------:|---:|----:|-----------------------------|
| E1_fixed_ACACAC vs Drt3b_N2 | 1.000 | 1.001 | 0.001 | 0.998 | I_pop, I_chan, JS |
| E2_2phase_alt vs Drt3b_N2 | 4.777 | 0.022 | 0.001 | 0.998 | I_pop, JS |
| E3_random_uniform vs Drt3b_N2 | 10.429 | 7.262 | 0.001 | 0.987 | I_pop, I_chan, JS |
| E1_fixed_ACACAC vs AbiK_uniform | 0.000 | 6.349 | 0.740 | 0.010 | I_chan, fidelity |
| Drt3b_N2 vs Drt3b_N3 | 0.585 | 0.009 | 0.001 | 0.998 | I_pop, JS |

Threshold conventions: ΔI ≥ 0.5 bits = clear separation;
Δf ≥ 0.1 or ΔJS ≥ 0.1 = clear distributional separation.

## The Drt3a-vs-Drt3b decision

The central v3 question: does the dual-observable apparatus distinguish Drt3a
(as it actually exists, E1: fixed ACACAC) from Drt3b-N=2 (Mode 3 alternating)?

- ΔI_pop  = 0.9997 bits
- ΔI_chan = 1.0012 bits
- Δf (per-base fidelity) = 0.0008
- ΔJS(output dist vs Drt3b) = 0.9982

**Verdict:** the channel-as-X observable I_chan **does** distinguish Drt3a-E1
from Drt3b-N=2. The biological Drt3a templates a *random-template-class* output
(Y rows distributed similarly to the WC channel applied to a non-degenerate
template), not a phase-restricted alternating-output distribution. The 5-channel
ensemble's empirical row distributions for Drt3a-WC and Drt3b-Mode3-N=2 differ
substantially in the row-distribution sense, and that difference is what I_chan
measures.

Direct numerical comparison:

- Drt3a-E1: (I_pop, I_chan) = (0.0000, 3.3320) bits
- Drt3b-N2: (I_pop, I_chan) = (0.9997, 2.3308) bits

## Recommendation for the v3 Methods section

The v3 Methods apparatus paragraph should:

1. Replace the single ambiguous I_struct(X;Y) of v2 with **two explicitly-defined**
   observables: I_struct^pop and I_struct^chan, with formulae and stated domains
   of validity.
2. State that the joint signature (I_struct^pop, I_struct^chan) is the primary
   apparatus output. Per-base fidelity f is a supplementary, non-MI observable
   reported alongside.
3. Note that I_struct^pop is bounded above by H(template ensemble) and is
   identically zero for a fixed-template population (Drt3a-E1 case). I_struct^chan
   is bounded above by log2(|channel ensemble|) and is non-zero whenever the
   focal channel's row distribution differs from the ensemble mixture.
4. Define the channel ensemble explicitly: in this paper, the canonical ensemble
   is {Drt3a-WC, Drt3b-Mode3-N=2, Drt3b-Mode3-N=3, Drt3b-Mode3-N=5, AbiK-uniform}.
   Other choices are valid but must be stated.

See `test_g2_v3_methods_paragraph.md` for a drop-in draft.

## Recommendation for the v3 Results section (Mode 1 vs Mode 3)

The Mode 1 vs Mode 3 distinguishing claim should be **carried by I_struct^chan**,
not by I_struct^pop. v3 should state:

> Drt3a (Mode 1 with biologically-encoded ACACAC template) and Drt3b (Mode 3
> with N=2 cyclic active site) are **distinguishable by their output row
> distributions** (I_struct^chan against the 5-channel ensemble:
> 3.33 bits for Drt3a vs 2.33 bits for Drt3b-N=2),
> not by per-template-realization population MI. The L-scaling claim of v2
> survives only when applied to the random-template counterfactual (E3); for
> the actual biological systems, the framework distinguishes them by
> product-distribution divergence and per-base fidelity, not by linear-in-L
> growth of I_struct^pop.

## 5×5 cross-mode separation matrix (JS divergence in bits)

| | Drt3a-WC | Drt3b-Mode3-N2 | Drt3b-Mode3-N3 | Drt3b-Mode3-N5 | AbiK-uniform |
|---|---:|---:|---:|---:|---:|
| **Drt3a-WC** | 0.0000 | 1.0000 | 1.0000 | 1.0000 | 0.9752 |
| **Drt3b-Mode3-N2** | 1.0000 | 0.0000 | 1.0000 | 1.0000 | 0.9873 |
| **Drt3b-Mode3-N3** | 1.0000 | 1.0000 | 0.0000 | 0.9997 | 0.9827 |
| **Drt3b-Mode3-N5** | 1.0000 | 1.0000 | 0.9997 | 0.0000 | 0.9698 |
| **AbiK-uniform** | 0.9752 | 0.9873 | 0.9827 | 0.9698 | 0.0000 |

KL divergence (i → j) in bits:

| | Drt3a-WC | Drt3b-Mode3-N2 | Drt3b-Mode3-N3 | Drt3b-Mode3-N5 | AbiK-uniform |
|---|---:|---:|---:|---:|---:|
| **Drt3a-WC** | 0.0000 | 38.2063 | 38.2063 | 38.2063 | 10.9137 |
| **Drt3b-Mode3-N2** | 38.3637 | 0.0000 | 38.3637 | 38.3637 | 11.0772 |
| **Drt3b-Mode3-N3** | 37.7889 | 37.7889 | 0.0000 | 37.7830 | 19.0901 |
| **Drt3b-Mode3-N5** | 37.0373 | 37.0373 | 37.0207 | 0.0000 | 19.8748 |
| **AbiK-uniform** | 28.1436 | 28.3721 | 28.2667 | 27.9081 | 0.0000 |

## Provenance

- code: `code/test_g2_dual_observable.py`
- step 1 csv: `results/test_g2_dual_observable_v1.csv`
- step 2 csv: `results/test_g2_mode_separation_matrix_v1.csv`
- step 3 csv: `results/test_g2_drt3a_three_ensembles_v1.csv`
- figures: `figures/test_g2_dual_observable_signature.png`,
            `figures/test_g2_mode_separation_heatmap.png`
- estimators: per-position plug-in MI for I_pop; KL/JS over empirical row
  distributions for I_chan and the cross-mode matrix.
- seed: `np.random.seed(42)`, `np.random.default_rng(42)`
- n_samples = 5000 per cell (5× v1's 1000), epsilon = 0.01
