# Dispatch G2: Apparatus repair — population MI and channel-as-X side by side (Test G2)

## Working directory
`/storage/kiran-stuff/templating_framework/`

## Goal
Test G v1 found that the I_struct apparatus, computed with X = template
realization, gives ~0 bits for Drt3a's actual ACACAC template (Reading 1)
and cannot distinguish Drt3a Mode 1 from Drt3b Mode 3. The two readings
documented in `results/test_g_drt3a_boundary_v1.md` propose two
complementary observables. Test G2 implements BOTH and reports them
jointly so the data shows which channel separates which mode pair under
which template ensemble. The two readings are not competing repairs;
they are complementary observables that together constitute the
refined apparatus.

This is Test G2. It supersedes Test G v1's resolution-doc framing of
"two readings, pick one." The v3 paper will report both observables.

## Background

Read these first:
1. `results/test_g_drt3a_boundary_v1.md` — full statement of the two
   readings and the empirical sweep that produced them.
2. `code/test_g_drt3a_boundary.py` — Test G v1 implementation. Test G2
   reuses the simulation harness; only the observables change.
3. `paper/drafts/draft_v2.md` lines 35–46 (apparatus subsection),
   lines 51 (Mode 1 scaling claim), line 61 (Mode 1 vs Mode 3
   distinguishing claim). These are the claims the apparatus refinement
   has to support or correct.
4. `framework_formal_v1.md` — apparatus definition. Note: the formal
   spec uses I_struct(X;Y) without specifying whether X is template
   realization or template-generating channel. G2 makes both
   interpretations explicit.

## The two observables

### Observable 1: I_struct^pop (population mutual information, Reading 1)

X = template realization (a specific ACACAC, or a specific random
sequence). For a population of templates, I_struct^pop is the standard
mutual information across the population:

  I_struct^pop = H(Y) − H(Y | X = x), averaged over x in the template population

For a population of one fixed template, I_struct^pop ≡ 0 by construction.
For a population of K equiprobable templates, I_struct^pop ≤ log2 K.
This is what Test G v1 measured.

### Observable 2: I_struct^chan (channel-as-X mutual information, Reading 2 reformulated)

X = the Markov channel itself (transition matrix, conformational state
sequence, gate identity). For Drt3a Mode 1, X = the Watson-Crick channel
parameterized by template length L and ε. For Drt3b Mode 3, X = the
2-state cyclic channel parameterized by ε.

For a templating system whose channel is fixed (one organism's enzyme,
one ncRNA), I_struct^chan against an alternative-channel ensemble
(other channels of comparable kind) measures how much the channel
identity determines product distribution. Operationally:

  I_struct^chan = H(Y | channel-class) − H(Y | specific channel)

where "channel-class" is the population of channels under consideration
(e.g., {Drt3a-WC, Drt3b-Mode3, AbiK-uniform}) and "specific channel"
is one member.

This recovers the per-base fidelity intuition without abandoning mutual
information: each channel determines a distinct product distribution,
and the MI across the channel ensemble measures how much information
the channel identity carries about the product.

For Drt3a (ACACAC, ε=0.01): Y = poly(GT) with high specificity;
distinct from poly(AC) (Drt3b output) and from random (AbiK output).
The channel identity carries log2 3 ≈ 1.585 bits about Y given a
3-channel ensemble, and the per-position resolution is what gives
the L-scaling: at L=6 against a 3-channel ensemble, I_struct^chan
should approach log2 3 × 6 / 6 = 1.585 bits if all positions equally
identify the channel, or 6 × something less if some positions are
ambiguous. The exact predicted value is what G2 measures.

### Joint reporting

Every test case is reported with both observables:

(template_type, mode_label, L, I_struct^pop, I_struct^chan, separation_pop, separation_chan, periodicity, marginal_freqs)

The product (I_struct^pop, I_struct^chan) is the refined apparatus
signature. Mode classification uses the joint signature, not either
observable alone.

## Concrete steps

### Step 1: Reproduce Test G v1 sweep with both observables

Reuse `code/test_g_drt3a_boundary.py` simulation. For every (L_T,
template_type) cell already run in v1:
- Compute I_struct^pop as before (preserves v1 numbers).
- Compute I_struct^chan against the channel ensemble:
  {Drt3a-WC, Drt3b-Mode3-N=2, Drt3b-Mode3-N=3, Drt3b-Mode3-N=5, AbiK-uniform}.
  This 5-channel ensemble is fixed across L_T.

Output: `results/test_g2_dual_observable_v1.csv` with columns
template_type, L_T, I_pop, I_chan, sep_pop, sep_chan, period_peak,
margin_A, margin_C, margin_G, margin_T.

### Step 2: Cross-mode separation matrix

For each pair of modes in {Drt3a-WC, Drt3b-Mode3-N=2, Drt3b-Mode3-N=3,
Drt3b-Mode3-N=5, AbiK-uniform}, compute the cross-mode separation:
how distinguishable are their product distributions under each
observable?

Output: `results/test_g2_mode_separation_matrix_v1.csv` (5x5 matrix of
JS divergences, plus a parallel matrix of I_struct^chan contributions).

The matrix entries should show:
- Drt3a (alt template) vs Drt3b (Mode 3 N=2): under I_struct^pop,
  ~0 separation (both produce same alternating output). Under
  I_struct^chan, separable because the channel distributions differ
  (Mode 1 has a template position descriptor; Mode 3 has a phase index).
- Drt3a vs AbiK: under both observables, large separation.
- Drt3b N=2 vs Drt3b N=3: both observables separate them by output period.

### Step 3: Drt3a-specific apparatus signature

Compute the full apparatus signature for Drt3a under three template
ensembles:
- E1: single fixed ACACAC template (the actual biological case)
- E2: 2-phase population (ACACAC, CACACA)
- E3: random uniform 6-nt templates (counterfactual baseline)

For each ensemble, report:
- I_struct^pop: the population mutual information
- I_struct^chan: the channel-as-X mutual information (against the
  5-channel ensemble from Step 1)
- per-base fidelity f_i = P(Y_i = WC(X_i) | X_i): the per-position
  determinism (always 0.99 for Drt3a regardless of ensemble)
- product distribution divergence from Drt3b output: JS(Y_Drt3a || Y_Drt3b)

### Step 4: Decision matrix for v3 apparatus claim

Write `results/test_g2_apparatus_decision_v1.md` with:

1. The full joint-observable signature for each system (Drt3a-E1, Drt3a-E2,
   Drt3a-E3, Drt3b-N=2, Drt3b-N=3, AbiK).
2. A statement of which observable separates which pair.
3. A recommendation for the v3 Methods section: which apparatus claim
   to make, with what qualifications.
4. A recommendation for the v3 Results section on Mode 1 vs Mode 3:
   which observable carries the distinguishing claim.

The decision should be data-driven. The v1 conclusion was "the apparatus
cannot distinguish Drt3a from Drt3b." G2 tests whether adding the
channel-as-X observable repairs this without abandoning mutual
information.

### Step 5: Apparatus claim for the v3 paper

Generate `results/test_g2_v3_methods_paragraph.md` containing a draft
paragraph for v3 Methods that:
- Defines the two observables explicitly with formulae.
- States the conditions under which each is well-defined.
- States the joint signature (I_struct^pop, I_struct^chan, per-base
  fidelity) as the refined apparatus.
- Identifies which observable carries which classification claim.

This paragraph IS the v3 apparatus repair. Get it right.

## Output artifacts

- `code/test_g2_dual_observable.py`
- `results/test_g2_dual_observable_v1.csv`
- `results/test_g2_mode_separation_matrix_v1.csv`
- `results/test_g2_drt3a_three_ensembles_v1.csv`
- `results/test_g2_apparatus_decision_v1.md`
- `results/test_g2_v3_methods_paragraph.md`
- `figures/test_g2_dual_observable_signature.png` (joint scatter of
  I_pop vs I_chan for all test cases, color-coded by mode)
- `figures/test_g2_mode_separation_heatmap.png` (5x5 cross-mode
  separation matrix as heatmap)

Update `results/CANONICAL_RESULTS.md` to add Test G2; mark Test G v1 as
"superseded by G2 for v3, retained for provenance."

## Constraints

- Naming: `test_g2` prefix, `_v1` suffix.
- Do not modify Test G v1 outputs; G2 supplements rather than replaces.
- Use the same Markov-chain harness as G v1 — only the observables change.
- ε = 0.01 for all systems unless explicitly varied; n_samples = 5000
  per cell (5x v1's 1000 because the cross-channel ensemble needs more
  resolution to separate close cases).
- Wall-time budget: 5 channels × 6 L values × 5000 samples ≈ small. <5
  minutes total.

## Reporting

Final reply should include:
1. The full joint-observable signature table for the 6 test cases.
2. The cross-mode separation matrix (5x5).
3. Whether the channel-as-X observable distinguishes Drt3a from Drt3b
   (yes/no with numbers).
4. The recommended v3 apparatus formulation, in one paragraph.
5. Path to the v3 Methods paragraph draft.
