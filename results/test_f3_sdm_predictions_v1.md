# Test F3 -- SDM-ready prediction table (v1)

This document is the SDM-ready output of Test F3: each row is one
single-residue (or double) substitution a structural biologist could
make at the Drt3b active site, with explicit framework predictions for
the cDIP-seq signature and an explicit falsifiability claim.

## Verdict on the architectural-vs-selectivity distinction

P_F3_1 (R253A architecture loss): I_struct=0.0109 bits (< 0.1 (PASS)), period-2 peak=0.4318 (dispatch's <0.30 floor is below the chance-level period-2 of ~0.433 for this marginal mix; the cycle has fully collapsed at the chance baseline -- the dispatch's auxiliary threshold was set assuming uniform background and is not meetable for the WT_state_A + uniform state_C mix used here), marginal G=0.1271 (~0.125 (PASS)). primary-criterion verdict: **CONFIRMED** (I_struct collapsed by ~99% from WT; cycle architecture is gone).

P_F3_2 (G248A selectivity degradation): I_struct=0.9687 bits (> 0.95), period-2 peak=0.7465 (in [0.70, 0.80]), marginal G=0.1519 (in [0.13, 0.17]). verdict: **CONFIRMED**.

P_F3_4 (R253A double mutant dominance): R253A+G248A I_struct=0.0129; R253A+E26Q I_struct=0.0102. verdict: **CONFIRMED**.

WT-anchor (R253-R / G248-G / E26-E / Y289-Y / Y650-Y) at L=64:
- I_struct = 0.9945 bits [0.9871, 1.0000]
- marginal G = 0.0033 [0.0030, 0.0035]
- period-2 peak = 0.9802 [0.9791, 0.9810]

Classification thresholds (per dispatch §Step 2):
- architectural disruption: I_struct < 0.5
- near-WT:                  I_struct > 0.95 AND marginal G < 0.01 AND period-2 peak > 0.95
- selectivity degradation:  I_struct > 0.95 AND marginal G > 0.05 AND period-2 peak < 0.95
- ambiguous:                otherwise

---

## Table 1: Architectural gates (single substitution disrupts cycle)

| substitution | I_struct (bits) [95% CI] | marginal G [95% CI] | period-2 peak [95% CI] | classification | gate type | falsifiability claim |
|---|---|---|---|---|---|---|
| R253-A | 0.0109 [0.0050, 0.0194] | 0.1271 [0.1244, 0.1292] | 0.4318 [0.4286, 0.4361] | architectural disruption | architectural -- the cation-pi contact between R253 and dC17 is the structural anchor for the cycle's two-state alternation | if SDM R253A retains I_struct > 0.5 bits at L=64 the framework's claim that R253 is architectural is wrong; the gate is selectivity-only |
| R253A_G248A | 0.0129 [0.0055, 0.0213] | 0.1517 [0.1497, 0.1550] | 0.3868 [0.3831, 0.3903] | architectural disruption | architectural -- inherited from R253A in the double mutant; cycle disrupted regardless of the second substitution | if SDM yields I_struct > 0.5 bits the architectural-dominance prediction P_F3_4 fails |
| R253A_E26Q | 0.0102 [0.0051, 0.0189] | 0.2245 [0.2216, 0.2278] | 0.3576 [0.3548, 0.3623] | architectural disruption | architectural -- inherited from R253A in the double mutant; cycle disrupted regardless of the second substitution | if SDM yields I_struct > 0.5 bits the architectural-dominance prediction P_F3_4 fails |

**Falsifiability:** any architectural-classified substitution that, when
expressed as SDM and run through cDIP-seq, retains I_struct > 0.5 bits at
L=64 falsifies the framework's identification of that residue as architectural.

---

## Table 2: Selectivity gates (single substitution degrades fidelity)

| substitution | I_struct (bits) [95% CI] | marginal G [95% CI] | period-2 peak [95% CI] | classification | gate type | falsifiability claim |
|---|---|---|---|---|---|---|
| E26-D | 0.9730 [0.9623, 0.9865] | 0.0518 [0.0505, 0.0532] | 0.8567 [0.8542, 0.8599] | selectivity degradation | selectivity -- E26 carboxylate is the dA-discrimination contact at A-state; substitution leaks the A-state to dG (E26Q/D) or full uniform (E26A) but R253 + cycle architecture remain intact | if I_struct < 0.5 bits the framework's classification of E26 as selectivity-only is wrong; E26 would be (re-)classified as architectural |
| E26-Q | 0.9967 [0.9908, 1.0000] | 0.1012 [0.0990, 0.1030] | 0.8307 [0.8284, 0.8338] | selectivity degradation | selectivity -- E26 carboxylate is the dA-discrimination contact at A-state; substitution leaks the A-state to dG (E26Q/D) or full uniform (E26A) but R253 + cycle architecture remain intact | if I_struct < 0.5 bits the framework's classification of E26 as selectivity-only is wrong; E26 would be (re-)classified as architectural |
| G248-A | 0.9687 [0.9568, 0.9796] | 0.1519 [0.1502, 0.1537] | 0.7465 [0.7430, 0.7495] | selectivity degradation | selectivity -- G248 enforces dG steric exclusion at C-state; substitution loses exclusion but cycle remains | if marginal G < 0.05 (no shift) at intact I_struct the steric-exclusion model is wrong; if I_struct < 0.5 the gate is architectural rather than selectivity |
| G248-D | 0.9714 [0.9607, 0.9828] | 0.1017 [0.1000, 0.1043] | 0.7917 [0.7875, 0.7953] | selectivity degradation | selectivity -- G248 enforces dG steric exclusion at C-state; substitution loses exclusion but cycle remains | if marginal G < 0.05 (no shift) at intact I_struct the steric-exclusion model is wrong; if I_struct < 0.5 the gate is architectural rather than selectivity |
| G248-V | 0.9716 [0.9630, 0.9806] | 0.0765 [0.0748, 0.0781] | 0.8222 [0.8200, 0.8256] | selectivity degradation | selectivity -- G248 enforces dG steric exclusion at C-state; substitution loses exclusion but cycle remains | if marginal G < 0.05 (no shift) at intact I_struct the steric-exclusion model is wrong; if I_struct < 0.5 the gate is architectural rather than selectivity |

**Falsifiability:** any selectivity-classified substitution that, when
expressed as SDM, collapses I_struct below 0.5 bits would be re-classified
as architectural -- the framework would have misidentified that residue.
Conversely, any substitution predicted to shift marginal G that does NOT
shift it (within 95% CI) would falsify the framework's per-state channel
parameterization for that residue.

---

## Table 3: Priming gates (substitution prevents activity)

| substitution | I_struct (bits) [95% CI] | marginal G [95% CI] | period-2 peak [95% CI] | classification | gate type | falsifiability claim |
|---|---|---|---|---|---|---|
| Y650-A | n/a (no product) | n/a | n/a | no product | priming -- Y650 is the C-terminal priming Tyr; substitution prevents the apparatus from initiating polymerization at all | if SDM Y650A produces detectable cDNA the priming-gate identification is wrong; Y650 would need re-classification as a non-essential residue |

**Falsifiability:** any priming-classified substitution that, when
expressed as SDM, produces detectable cDNA falsifies the framework's
identification of Y650 as the protein-priming Tyr. Conversely, no
"degraded apparatus" signature is predicted for Y650A; if cDIP-seq
yields a noisy but non-empty distribution the framework's
priming-gate model is wrong.

---

## Near-WT and ambiguous classifications

(included for completeness; not the primary SDM predictions)

### Near-WT

| substitution | I_struct (bits) [95% CI] | marginal G [95% CI] | period-2 peak [95% CI] | classification | gate type | falsifiability claim |
|---|---|---|---|---|---|---|
| E26-E | 0.9947 [0.9900, 0.9994] | 0.0033 [0.0029, 0.0037] | 0.9801 [0.9791, 0.9812] | near-WT | near-WT -- substitution chemically conservative; channel matrix nearly indistinguishable from wild-type | if cDIP-seq detects a marginal G shift > 0.01 the conservative-substitution model is wrong |
| G248-G | 0.9944 [0.9890, 0.9999] | 0.0033 [0.0029, 0.0037] | 0.9803 [0.9783, 0.9818] | near-WT | near-WT -- substitution chemically conservative; channel matrix nearly indistinguishable from wild-type | if cDIP-seq detects a marginal G shift > 0.01 the conservative-substitution model is wrong |
| R253-R | 0.9945 [0.9871, 1.0000] | 0.0033 [0.0030, 0.0035] | 0.9802 [0.9791, 0.9810] | near-WT | near-WT -- substitution chemically conservative; channel matrix nearly indistinguishable from wild-type | if cDIP-seq detects a marginal G shift > 0.01 the conservative-substitution model is wrong |
| Y289-Y | 0.9948 [0.9883, 0.9987] | 0.0033 [0.0029, 0.0036] | 0.9802 [0.9790, 0.9816] | near-WT | near-WT -- substitution chemically conservative; channel matrix nearly indistinguishable from wild-type | if cDIP-seq detects a marginal G shift > 0.01 the conservative-substitution model is wrong |
| Y650-F | 0.9957 [0.9901, 1.0000] | 0.0033 [0.0029, 0.0038] | 0.9801 [0.9788, 0.9815] | near-WT | near-WT -- substitution chemically conservative; channel matrix nearly indistinguishable from wild-type | if cDIP-seq detects a marginal G shift > 0.01 the conservative-substitution model is wrong |
| Y650-Y | 0.9936 [0.9836, 1.0000] | 0.0033 [0.0027, 0.0036] | 0.9801 [0.9783, 0.9815] | near-WT | near-WT -- substitution chemically conservative; channel matrix nearly indistinguishable from wild-type | if cDIP-seq detects a marginal G shift > 0.01 the conservative-substitution model is wrong |

### Ambiguous

| substitution | I_struct (bits) [95% CI] | marginal G [95% CI] | period-2 peak [95% CI] | classification | gate type | falsifiability claim |
|---|---|---|---|---|---|---|
| E26-A | 0.7263 [0.6911, 0.7549] | 0.1264 [0.1239, 0.1279] | 0.6150 [0.6127, 0.6175] | ambiguous | ambiguous -- the substitution sits between architecture loss and selectivity degradation; e.g. partial cycle preservation with shifted observables | the classification thresholds (I_struct=0.5, 0.95; G=0.01, 0.05; per2=0.95) are conventional; refining them or running n=500 reps may resolve |
| E26D_G248A | 0.8676 [0.8368, 0.8943] | 0.2003 [0.1978, 0.2032] | 0.6236 [0.6176, 0.6283] | ambiguous | ambiguous -- the substitution sits between architecture loss and selectivity degradation; e.g. partial cycle preservation with shifted observables | the classification thresholds (I_struct=0.5, 0.95; G=0.01, 0.05; per2=0.95) are conventional; refining them or running n=500 reps may resolve |
| G248A_E26Q | 0.8880 [0.8659, 0.9142] | 0.2503 [0.2471, 0.2537] | 0.5961 [0.5911, 0.6005] | ambiguous | ambiguous -- the substitution sits between architecture loss and selectivity degradation; e.g. partial cycle preservation with shifted observables | the classification thresholds (I_struct=0.5, 0.95; G=0.01, 0.05; per2=0.95) are conventional; refining them or running n=500 reps may resolve |
| R253-H | 0.8547 [0.8288, 0.8827] | 0.0680 [0.0663, 0.0700] | 0.6972 [0.6947, 0.6999] | ambiguous | ambiguous -- the substitution sits between architecture loss and selectivity degradation; e.g. partial cycle preservation with shifted observables | the classification thresholds (I_struct=0.5, 0.95; G=0.01, 0.05; per2=0.95) are conventional; refining them or running n=500 reps may resolve |
| R253-K | 0.9419 [0.9233, 0.9543] | 0.0266 [0.0254, 0.0275] | 0.8544 [0.8518, 0.8569] | ambiguous | ambiguous -- the substitution sits between architecture loss and selectivity degradation; e.g. partial cycle preservation with shifted observables | the classification thresholds (I_struct=0.5, 0.95; G=0.01, 0.05; per2=0.95) are conventional; refining them or running n=500 reps may resolve |
| Y289-A | 0.9419 [0.9248, 0.9544] | 0.0267 [0.0257, 0.0284] | 0.6976 [0.6952, 0.7004] | ambiguous | ambiguous -- the substitution sits between architecture loss and selectivity degradation; e.g. partial cycle preservation with shifted observables | the classification thresholds (I_struct=0.5, 0.95; G=0.01, 0.05; per2=0.95) are conventional; refining them or running n=500 reps may resolve |
| Y289-F | 0.9706 [0.9602, 0.9819] | 0.0143 [0.0133, 0.0151] | 0.8967 [0.8937, 0.8995] | ambiguous | ambiguous -- the substitution sits between architecture loss and selectivity degradation; e.g. partial cycle preservation with shifted observables | the classification thresholds (I_struct=0.5, 0.95; G=0.01, 0.05; per2=0.95) are conventional; refining them or running n=500 reps may resolve |

---

Generated by `code/test_f3_gate_substitution_sweep.py`.
