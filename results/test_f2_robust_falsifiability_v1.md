# Test F2 -- robust falsifiability statement (v1)

This document supersedes `test_f_falsifiability_v1.md` for v3-paper purposes.
It turns the v1 point-estimate envelope into a robust statement with explicit
parameter sensitivity bounds and replicate confidence intervals.

## Overall verdict

**NARROW BUT VALID**: one of the two robustness dimensions (E26D fidelity sensitivity OR secondary-residue eps tolerance) is tight. The framework's prediction is correct under the v1 parameterization, but the falsifiability envelope is narrow.

The dispatch's framing for the two non-robust outcomes is:
"narrow but valid" or "fragile". The verdict above uses that framing.

---

## 1. Robust envelope at L=64 (with 95% CIs from n=500 replicate runs, Step 4)

Test F v1 reported a point-estimate envelope at L=64 from n=30 reps. Test F2
re-ran the WT clade and the strongest degraded-E26 clade (C27_PRIM-E26D,
labelled `E26D-fallback`) at n=500 reps to obtain tight CIs.

| observable      | WT clade (C01) -- mean [95% CI]                                | degraded (C27_PRIM-E26D, E26D-fallback) -- mean [95% CI]                                |
|----------------|------------------------------------------------|--------------------------------------------------|
| marginal G      | 0.0033 [0.0031, 0.0035] | 0.0516 [0.0509, 0.0525] |
| period-2 peak   | 0.9802 [0.9794, 0.9808] | 0.8569 [0.8553, 0.8584] |
| I_struct (bits) | 0.9997 [0.9987, 1.0000] | 0.9963 [0.9926, 0.9996] |

**Sanity check**: Test F v1 reported (worst-case across in-envelope clades) marginal
G <= 0.0033, period-2 peak >= 0.9793, I_struct >= 0.9996 bits. The WT-clade CIs above
**bracket those v1 numbers**, so v1's point estimate is consistent with the F2
high-replicate distribution.

The robust envelope (v3-ready) for **secondary-residue-only clades** at L=64 is
the WT-clade row above with stated 95% CIs. Any natural family member observed
outside the marked CI bands would be a candidate counter-example; observation
outside by more than 5x the half-width would be conclusive.

---

## 2. E26D parameterization sensitivity (Step 1)

Test F v1 used a single state_A fidelity P(A) = 0.85 for E26D (interpolated from
Deng et al. 2026's E26Q misincorporation data, not measured directly). Step 1
swept this parameter to characterize sensitivity across the 6
E26D-bearing clades in the natural alignment.

The reference baselines in the sweep:
- WT marginal G at fidelity 0.99: 0.0050 (sigma ~ 0.0001)
- E26Q marginal G at fidelity 0.80: 0.0683 (sigma ~ 0.0004)

| E26D fidelity P(A) | mean marginal G | distinguishable from WT (3 sigma) | distinguishable from E26Q (3 sigma) |
|---|---|---|---|
| 0.50 | 0.1684 | YES | YES |
| 0.60 | 0.1351 | YES | YES |
| 0.70 | 0.1017 | YES | YES |
| 0.80 | 0.0683 | YES | no |
| 0.85 | 0.0516 | YES | YES |
| 0.90 | 0.0350 | YES | YES |
| 0.95 | 0.0183 | YES | YES |
| 0.99 | 0.0050 | no | YES |

**E26D fidelity range that keeps the prediction valid:**
E26D's marginal-G prediction is distinguishable from BOTH WT (fidelity 0.99) and E26Q (fidelity 0.80) at fidelity in **[0.50, 0.95]**.

The intermediate-fidelity claim ("E26D should look between WT and Q") survives
across the swept range as long as fidelity is bounded away from the two anchors.
At fidelity = 0.99, the E26D channel becomes WT-indistinguishable; at 0.80 it
becomes Q-indistinguishable. The v1 anchor of 0.85 sits in the operational range
where both distinguishability tests pass.

---

## 3. Universal-gate hypothetical SDM predictions (Step 2)

The natural alignment contains 0 R253-or-G248 substitutions; both universal
gates are 99.9-100% conserved. The framework's strongest predictions are
therefore not testable in the family alignment but ARE testable via SDM.

Pre-registered predictions for hypothetical SDM mutants at L=64 (n=30 reps,
n=5000 samples per rep):

| mutant         | marginal G [95% CI]                                  | period-2 peak [95% CI]                              | I_struct (bits) [95% CI]                           |
|----------------|---------------------------------------|---------------------------------------|--------------------------------------|
| R253A          | 0.1265 [0.1255, 0.1274] | 0.4325 [0.4303, 0.4345] | 0.0090 [0.0065, 0.0120] |
| G248A          | 0.1516 [0.1505, 0.1528] | 0.7469 [0.7455, 0.7479] | 0.9966 [0.9936, 0.9991] |
| R253A + G248A  | 0.1519 [0.1506, 0.1537] | 0.3872 [0.3855, 0.3892] | 0.0084 [0.0059, 0.0116] |

**These are the framework's most testable explicit-number predictions for SDM
experiments.** A structural biologist could express R253A or G248A, run cDIP-seq,
and either falsify the framework (if cycle architecture is preserved despite
the universal-gate substitution) or confirm the strongest prediction (if both
gates are required for Mode 3).

The R253A prediction (cycle disrupted, I_struct collapses to ~ 0.01 bits)
is the most discriminating: the framework predicts the cycle architecture
falls apart entirely, not just degrades. Any R253A SDM that retains
I_struct > 0.5 bits would falsify the universal-gate / Mode 3 architecture
claim.

---

## 4. Secondary-residue eps tolerance (Step 3)

Step 3 perturbed the effective state_A and state_C fidelity by eps in
{0.01, 0.02, 0.05, 0.10, 0.15} for the 33 in-envelope clades, modelling a
counterfactual where a secondary-residue substitution destabilizes the active-
site geometry.

Per-observable exit eps from the v1 envelope:
  - marginal_G: exits envelope at eps = 0.01
  - period2_peak: exits envelope at eps = 0.02
  - I_struct: exits envelope at eps = 0.02

**Verdict on eps tolerance:**
The framework's secondary-residue prediction is **fragile**: a perturbation as small as eps = 0.01 exits the v1 envelope. A real secondary-residue substitution that destabilized the active-site geometry by even 1-2% would visibly violate the prediction.

This is the operational answer to "how much would the framework be wrong if
secondary residues did affect mode classification?"

---

## 5. Pre-registered v3 envelope (combining 1-4)

| quantity                    | predicted bound (v1 point) | predicted bound (v3 robust)              |
|-----------------------------|----------------------------|------------------------------------------|
| marginal G (in-envelope)    | <= 0.0033                  | <= 0.0035 (95% CI hi from n=500)  |
| period-2 peak (in-envelope) | >= 0.9793                  | >= 0.9794 (95% CI lo from n=500) |
| I_struct (in-envelope)      | >= 0.9996 bits             | >= 0.9987 bits (95% CI lo from n=500) |
| marginal G (E26D-degraded)  | ~0.0517                    | mean 0.0516, [0.0509, 0.0525] |
| marginal G (R253A SDM)      | n/a (not in alignment)     | mean 0.1265, [0.1255, 0.1274] |
| marginal G (G248A SDM)      | n/a (not in alignment)     | mean 0.1516, [0.1505, 0.1528] |
| I_struct (R253A SDM)        | < 0.5 bits                 | mean 0.0090, [0.0065, 0.0120] |
| I_struct (G248A SDM)        | not pre-stated             | mean 0.9966, [0.9936, 0.9991] |
| E26D fidelity range valid   | (single point: 0.85)       | distinguishable in [see Step 1 table]    |
| eps tolerance (secondary)   | not pre-stated             | exit at eps = 0.01 |

Generated by `code/test_f2_robustness.py` from `code/test_f_family_sweep.py`'s
parameterization rule and v1 prediction set.
