# Dispatch 2: Drt3a Mode 1 with degenerate template (Test G)

## Working directory
`/storage/kiran-stuff/templating_framework/`

## Goal
Resolve the Mode 1 / Mode 3 boundary case raised in the handoff. Drt3a
synthesizes poly(GT) using a 6-nt ACACAC region of the ncRNA as
Watson-Crick template. Naïvely, Mode 1 with L = 6 and alphabet 4 predicts
$I_\text{struct} \le 12$ bits. But the template is alternating AC: only
1 bit of phase information is non-redundant. This dispatch tests whether
the apparatus reports the naïve Mode 1 ceiling or the
information-content-limited ceiling, and what that says about how Mode 1
and Mode 3 relate when the Mode 1 template is fully degenerate.

This is Test G. Conceptual stakes are real: the framework distinguishes
modes by mechanism, not by output statistics, but the apparatus measures
output statistics. The boundary case stress-tests this gap.

## Background

Read these first:
1. `code/test_a1.py` and `code/test_a2.py` — Mode 1 sims with bulk-matched
   control. Test G reuses this skeleton.
2. `paper/drafts/draft_v2.md` Results sections on Mode 1 length-scaling
   and on Drt3 — Test G is positioned as a diagnostic, not a fifth
   biological anchor.
3. `framework_formal_v1.md` — apparatus definition: $I_\text{struct}(X;Y)
   = H(\phi_Y(Y)) - H(\phi_Y(Y) | \phi_X(X))$ with bulk-matched and
   structure-scrambled controls.

## Two competing predictions to test

Prediction A (mechanism-determined):
The apparatus measures $I_\text{struct}$ at the joint distribution level.
Even with a degenerate template, every product position is determined by a
template position via Watson-Crick. The mutual information equals
$H(\phi_Y(Y)) - H(\phi_Y(Y) | \phi_X(X))$ where the conditioning on X
collapses the product distribution at each position to a near-delta on the
WC complement. Predicted $I_\text{struct} \approx L \log_2 4 = 12$ bits at
L = 6, regardless of template degeneracy.

Prediction B (template-information-limited):
The template's intrinsic information content (counting phase as 1 bit) is
the ceiling on transferable information. Predicted $I_\text{struct} \le 1$
bit total because that is all the information the template carries.

The framework as currently written predicts A: the apparatus is a
descriptor-relative measure on joint distributions, and a low-entropy
template still mediates high-fidelity per-position transfer. The product
of templating an ACACAC template is GT GT GT — and conditional on
knowing the template, every product position is determined.

If the simulation confirms A and the apparatus gives ~12 bits at L = 6,
that resolves the boundary in favor of mechanism-based classification:
Drt3a is Mode 1 even with a degenerate template, because the apparatus
reads the per-position complementarity, not the template's information
content. Mode 3 with N = 2 saturates at 1 bit not because the template
is degenerate but because there is no L-scaling.

If the simulation gives B-like results (~1 bit), the apparatus is
inadvertently measuring template information rather than transfer
information, and the apparatus needs revision.

## Concrete steps

### Step 1: Build Drt3a templating sim

Markov-chain simulation:
- Template X = ACACAC (length 6)
- For each position i, Y[i] = WC_complement(X[i]) with fidelity 0.99,
  uniform misincorporation 0.01/3 to the other three nucleotides.
- Sweep over template length L_T ∈ {2, 4, 6, 8, 12, 24} where the template
  is the alternating AC pattern repeated to length L_T. For L_T > 6,
  treat as a longer alternating template (Drt3a's actual template is 6 nt
  via translocation, but this isolates the degeneracy effect from the
  length effect).
- Sweep over a comparison case: random-sequence template of the same
  length L_T (this is Test A.1 territory). Record both for comparison.

For each (L_T, template_type), run 1000 reps and compute:
- Empirical joint distribution P(Y | X)
- $I_\text{struct} = H(Y) - H(Y | X)$ using empirical entropies
- Per-position $I_\text{struct}$
- Bulk-matched control: scramble X positions, measure same observable
- Structure-scrambled control: use 2nd-order Markov scramble of Y from
  P(Y) marginals while preserving compositional and nearest-neighbor
  statistics

### Step 2: Compare alternating vs random templates

Plot $I_\text{struct}$ vs L_T for both template types. Key question:
Do the two curves coincide (Prediction A: apparatus measures transfer,
not template entropy) or diverge (Prediction B: alternating template
yields lower $I_\text{struct}$)?

### Step 3: Compare Drt3a (Mode 1, alternating template) vs Drt3b (Mode 3, N=2)

The product of both is poly(GT/AC). The apparatus must distinguish them
on mechanism, not output. Run a side-by-side comparison:
- Drt3a sim: Mode 1 with ACACAC template, output GTGTGT
- Drt3b sim: Mode 3 N=2, output ACACAC

Both should give nearly identical output statistics (alternating
two-letter pattern). The apparatus distinguishes them because in Drt3a,
$\phi_X$ is the template sequence (separate molecule, descriptor present);
in Drt3b, $\phi_X$ is the cyclic state index, not the template (no
separate molecule). The descriptor-relative apparatus reports
$I_\text{struct}$ relative to a different X for each case.

The framework predicts:
- Drt3a $I_\text{struct}$ scales linearly with L_T (Mode 1).
- Drt3b $I_\text{struct}$ saturates at 1 bit regardless of L (Mode 3, N=2).

Same output statistics, different scaling. This is the framework's
mechanism-versus-output distinction made operational.

### Step 4: Document the resolution

Write `results/test_g_drt3a_boundary_v1.md` with:
1. Statement of competing predictions A and B.
2. Empirical result.
3. Resolution: which prediction does the framework make, and which does
   the apparatus confirm?
4. Implication for Mode 1 / Mode 3 boundary: is the boundary mechanism-
   determined or output-determined?

If the simulation confirms Prediction A (most likely outcome under
correct apparatus behavior): the boundary is mechanism-determined, and
the apparatus reports per-position transfer correctly even with degenerate
templates. Drt3a is Mode 1 unambiguously.

If the simulation confirms Prediction B: the apparatus needs revision.
Document the failure mode and propose a fix (probably: condition on the
full template descriptor including phase, not just on the template
sequence).

## Output artifacts

- `code/test_g_drt3a_boundary.py`
- `results/test_g_alternating_vs_random_template_v1.csv`
- `results/test_g_drt3a_vs_drt3b_comparison_v1.csv`
- `results/test_g_drt3a_boundary_v1.md`
- `figures/test_g_alternating_vs_random.png` (line plot of $I_\text{struct}$
  vs L_T for both template types)
- `figures/test_g_drt3a_vs_drt3b_scaling.png` (line plot of $I_\text{struct}$
  vs L for both modes — should diverge)

Update `results/CANONICAL_RESULTS.md`.

## Constraints

- Naming: `test_g` prefix, `_v1` suffix.
- Do not modify any earlier test outputs.
- Use the same entropy-estimation conventions as Test A — finite-sample
  bias correction (e.g., Miller-Madow or NSB) for $H(Y)$ if used in
  Test A. Match it.
- Mode 1 fidelity: 0.99 per position with uniform misincorp, matching
  Test A.
- Wall-time budget: small. Should finish in <5 minutes.

## Reporting

Final reply should include:
1. The empirical curve $I_\text{struct}$ vs L_T for alternating vs random
   templates (do they coincide? — yes/no with numbers).
2. The Drt3a vs Drt3b L-scaling comparison (linear vs saturating? — yes/no
   with slopes).
3. Which prediction the apparatus confirms.
4. Path to the resolution document.
