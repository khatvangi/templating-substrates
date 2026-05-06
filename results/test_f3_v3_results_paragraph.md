# Test F3 -- v3 Results paragraph (drop-in)

The following is a drop-in paragraph for the v3 manuscript Results section,
augmenting the F2 universal-gate findings with the F3 architecture-vs-
selectivity refinement. It is intended to follow the F2 paragraph that
reports R253A and G248A predictions diverged.

---

**Universal-gate residues divide into two qualitatively distinct functional
classes.** Test F3 refines the v2 universal-gate claim by characterizing the
unexpected divergence between R253A and G248A predictions reported in F2.
The two substitutions, both at residues conserved at >99.9% in the Drt3b
family alignment, generate qualitatively different apparatus signatures:
R253A collapses the cycle architecture entirely (I_struct =
0.0109 bits at L=64, period-2 peak =
0.4318, marginal G = 0.1271
-- consistent with a Mode-1-like averaged channel without templating),
whereas G248A leaves the cycle intact and merely degrades fidelity within
it (I_struct = 0.9687 bits, period-2 peak =
0.7465, marginal G = 0.1519).
This divergence is not a quantitative anomaly but a structural distinction:
**architectural gates** (R253) are residues whose substitution destroys the
N-state cycle itself, eliminating mode classification; **selectivity gates**
(G248, E26, Y289) are residues whose substitution degrades per-state
fidelity within an otherwise intact cycle. A third class, **priming
gates** (Y650), prevents apparatus initiation entirely and yields no product.

The architecture-vs-selectivity distinction generates orthogonal SDM
predictions. Architectural-gate substitutions predict cycle loss
(I_struct < 0.5 bits, with period-2 autocorrelation collapsing to the
chance baseline determined by the average state distribution -- here
~0.43 for the WT_state_A + uniform_state_C mix); selectivity-gate
substitutions predict an intact cycle with shifted observables
(I_struct > 0.95 bits, marginal G shifted away from chance, period-2
peak modestly degraded but still well above chance).
Test F3 generates explicit per-substitution predictions for the R253
series (R253K/H/A), the G248 series (G248A/V/D), the E26 series (E26D/Q/A),
the Y289 series (Y289F/A), the Y650 priming series, and four double-mutant
combinations spanning architectural and selectivity classes. The predictions
are recorded in `results/test_f3_sdm_predictions_v1.md` and are independently
testable by site-directed mutagenesis followed by cDIP-seq.

This refinement strengthens the framework's universal-gate claim rather than
weakening it: instead of a single class of "universal gates", the framework
now identifies three structurally distinct classes whose substitutions
predict three structurally distinct cDIP-seq signatures. R253A's predicted
collapse to I_struct = 0.0109 bits is the framework's
sharpest discriminating prediction; any SDM R253A retaining I_struct > 0.5
bits would falsify the architectural-gate identification.

---

Source data: `results/test_f3_gate_substitution_sweep_v1.csv`,
`results/test_f3_classification_v1.csv`,
`results/test_f3_other_primary_gates_v1.csv`,
`results/test_f3_double_mutants_v1.csv`. Detailed per-substitution
predictions: `results/test_f3_sdm_predictions_v1.md`.
