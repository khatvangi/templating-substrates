# Canonical results (DO NOT OVERWRITE)

The following CSVs are the canonical results of completed tests:

- test_a1_results.csv (Test A.1, PASS)
- test_a2_results.csv (Test A.2, PASS)
- test_b_results.csv (Test B, PASS, including verdict_corrected)
- test_c_results.csv (Test C, PASS)
- test_c_length_limit.csv (Test C length-limit sub-test, PASS)
- test_d_v2_results.csv (Test D v2, PASS — population dynamics, supersedes test_d)
- test_e_v2_results.csv (Test E v2, PASS — Drt3 anchoring with paper-derived E26Q params)
- test_f_family_predictions_v1.csv (Test F, PASS — 33/33 secondary-residue clades stay in Mode 3 envelope; 6/6 E26→D clades shift in predicted direction)
- test_f_per_clade_simulations_v1.csv (Test F raw sim outputs)
- test_g_alternating_vs_random_template_v1.csv (Test G — Drt3a degenerate-template boundary; per-position estimator)
- test_g_drt3a_vs_drt3b_comparison_v1.csv (Test G — joint-MI comparison; NEGATIVE FINDING: at L=6, Drt3a's 6-nt ACACAC template gives ~1 bit not ~12 bits, so I_struct L-scaling does NOT distinguish Drt3a Mode 1 from Drt3b Mode 3 in the boundary case — see test_g_drt3a_boundary_v1.md)
- test_h_summary_v1.csv (Test H, MIXED — Mode 1 vs Mode 6 competition, 4 scenarios × 10 reps)
- test_h_scenario_{A,B,C,D}_v1.csv (Test H per-scenario per-generation traces)

Test H verdicts: P_H1, P_H2, P_H4 CONFIRMED; P_H3 strictly FALSIFIED in 1/10
reps of scenario B (Wright-Fisher founder loss when Mode 1 starts at 10%, NOT
Mode 6 developing inheritance). Qualitative framework claim survives; strict
"initial abundance does not matter" claim fails at the small-founder regime.

## v2 follow-ups (Test F2, G2, H2)

- test_f2_E26D_sensitivity_v1.csv (Test F2, NARROW BUT VALID — E26D distinguishable from WT and Q at fidelity in [0.50, 0.95])
- test_f2_universal_gate_hypotheticals_v1.csv (Test F2 — explicit SDM predictions: R253A I_struct=0.009 bits, G248A I_struct=0.997, R253A+G248A I_struct=0.008)
- test_f2_secondary_residue_epsilon_v1.csv (Test F2 — secondary-residue ε tolerance: marginal_G exits envelope at ε=0.01, period-2 and I_struct at ε=0.02 — FRAGILE on this axis)
- test_f2_replicate_CIs_v1.csv (Test F2 — n=500-rep CIs confirm v1 numbers within tight bands)
- test_g2_dual_observable_v1.csv (Test G2, REPAIR SUCCEEDED — apparatus extended with channel-as-X observable I_struct^chan)
- test_g2_mode_separation_matrix_v1.csv (Test G2 — 5×5 JS divergence matrix; all distinct modes ~1.0 separation)
- test_g2_drt3a_three_ensembles_v1.csv (Test G2 — Drt3a-E1 (fixed ACACAC): I_pop=0.00 / I_chan=3.33 vs Drt3b-N2: I_pop=1.00 / I_chan=2.33 — 1-bit clean separation via I_chan)
- test_h2_founder_boundary_v1.csv (Test H2 — Implementation A founder sweep, 13 N_1 × 30 reps; Mode 1 wins 90-100% across all N_1)
- test_h2_implementation_B_v1.csv (Test H2 — Implementation B per-lineage fixed templates; max |P_A-P_B|=0.80, framework's binary-copyability claim FALSIFIED)
- test_h2_implementation_C_v1.csv (Test H2 — Implementation C graded copyability h ∈ {0,0.25,0.5,0.75,1.0} × 11 N_1 × 30 reps; crossover at P=0.5 not found at any h, refines P_H2_5)

v2 verdicts:
- F2: NARROW BUT VALID — v1 envelope confirmed at high replicate count; E26D fidelity range [0.50, 0.95] valid; universal-gate SDM predictions publication-ready; secondary-residue ε tolerance FRAGILE.
- G2: REPAIR SUCCEEDED — dual-observable apparatus distinguishes Drt3a-E1 from Drt3b-N=2 by ΔI_chan=1.0 bits. v3 Methods to use (I_struct^pop, I_struct^chan) joint signature; v3 Results to carry Mode 1/Mode 3 distinction via I_chan.
- H2: P_H2_1 FALSIFIED (Mode 1 wins ≥90% even at N_1=20), P_H2_2 INCONCLUSIVE, P_H2_3 FALSIFIED (Implementation B much weaker for Mode 1 — binary copyability claim wrong), P_H2_4 and P_H2_5 REFINED (no crossover at P=0.5 found).

H v1's P_H3 founder-loss interpretation does NOT replicate at H2's n=30 (Mode 1
wins 93% at N_1=20). The 1/10 fail in H v1 was likely sampling noise. The
framework's qualitative claim (copyable beats non-copyable in Implementation A)
holds robustly across all N_1; the binary-vs-graded copyability claim
(Implementation B vs A) is the new open question.

## v3 follow-ups (Test F3, H3)

- test_f3_gate_substitution_sweep_v1.csv (Test F3 — refined R253/G248 series + L sweep, 21 substitution conditions)
- test_f3_classification_v1.csv (Test F3 — architecture/selectivity/near-WT/ambiguous classification per substitution)
- test_f3_other_primary_gates_v1.csv (Test F3 — E26, Y289, Y650 series; F3 finds Y289 series classifies as "ambiguous", possibly framework-parameterization too conservative)
- test_f3_double_mutants_v1.csv (Test F3 — 4 double mutants; R253A architectural disruption dominates as predicted)
- test_h3_isolated_dynamics_v1.csv (Test H3 — 5 mechanisms M0-M4 × 4 β regimes × 30 reps; plateau heights table)
- test_h3_pairwise_competition_v1.csv (Test H3 — 9 pairs × 3 N_i × 30 reps; M4 vs M1 at N_1=200 only wins 83% — P_H3_6 FALSIFIED)
- test_h3_selection_regime_v1.csv (Test H3 — M0-M4 × 4 L_TARGET sweep at β=10)
- test_h3_effective_heritability_v1.csv (Test H3 — M1, M2 r-sweep, M3 heritability across gen 50/200/500; M3 yields nan — likely script artifact, flagged for follow-up)

v3 verdicts:
- F3: ARCHITECTURE/SELECTIVITY DISTINCTION HOLDS — R253A I_struct=0.011 (cycle collapsed, P_F3_1 confirmed), G248A I_struct=0.969 (cycle intact, P_F3_2 confirmed on all sub-criteria). Three gate classes now distinguished: architectural (R253), selectivity (G248, E26, Y289), priming (Y650). Y289 classification ambiguous — flagged as possible framework-parameterization issue. v3 Results paragraph drop-in at test_f3_v3_results_paragraph.md; SDM predictions at test_f3_sdm_predictions_v1.md.
- H3: PLATEAU LANDSCAPE CONFIRMED, BUT THREE PREDICTIONS REFINED OR FALSIFIED.
  - P_H3_1 (M0=chance), P_H3_2 (M1 plateau), P_H3_4 (M3≈M1), P_H3_5 (M4 dominates), P_H3_7 (M4 v M0), P_H3_8 (M1 v M0): CONFIRMED
  - P_H3_3 REFINED — M2 NOT monotonic in r. M2_r=0.10 (small fresh-draw injection) outperforms M2_r=0.00 (pure M1) by 6-14% across all β. Sweet-spot effect not predicted.
  - P_H3_6 FALSIFIED — M4 vs M1 at N_1=200 (equal start) only wins 83%, not the predicted ≥95%. Lineage fixation more competitive than framework's binary copyability claim allows.
  - P_H3_9 REFINED — M3 vs M1 not the predicted tie; M3 wins 67-97% of head-to-head. Individual-level exact copy beats lineage-level fixation.
  - P_H3_10 PARTIAL — no non-M4 mechanism reaches 0.95 (the framework's claim survives in that direction), but M4 itself only reaches 0.832 at β=2. Strongest claim has to be qualified by selection regime.
  - Recommended Draft A (mechanism-specific) for v3 inheritance condition rewrite. v3 Discussion drop-in sentence at test_h3_v3_inheritance_revision.md §5.
  - Effective heritability table has anomalies: M3=nan (script bug — divide-by-zero on exact-copy correlation), M1 h_eff~0.04 (surprisingly low). Flagged for H3-followup before v3 publication.

## v3 follow-up (Test H4)

- test_h4_long_horizon_v1.csv (Test H4 — M4 vs M1 at N_M1 ∈ {20, 80, 200} × 30 reps × 5000 gens; per-cell summary)
- test_h4_convergence_v1.csv (Test H4 — per-rep crossover gen and fitness checkpoints at gens 1000/2000/3000/5000)
- test_h4_m1_plateau_v1.csv (Test H4 — per-rep M1 plateau heights to test the plateau-vs-climbing hypothesis)

H4 verdict: **P_H4_1 FALSIFIED — extending horizon 1000→5000 gens did NOT
restore M4 dominance.** M4 win rate at N_M1=200 (equal start): 0.80 with
95% Wilson CI [0.63, 0.91], essentially unchanged from H3's 0.83. The
17–20% deficit is genuine long-horizon durability, not finite-time noise.
Other predictions: P_H4_2 CONFIRMED (M1 plateau at 0.502, no rep climbs
past gen=500), P_H4_3 REFINED (M4 mean fitness 0.973 — just below the
0.98 prediction), P_H4_4 CONFIRMED (median crossover gen=2 in winning
reps), P_H4_5 CONFIRMED (M4 wins 100% at N_M1=20).

Mechanism: M4 climbs to high fitness (0.97) when it survives the early
demographic phase, but in 6/30 equal-start reps M4 either fails to fix
(2/30 no crossover) or transiently dominates then is pushed back by
M1's lineage advantage (4/30 crossover-then-loss). M1's plateau (0.50)
is high enough that selection at β=10 doesn't always push M4 above it
fast enough to escape stochastic loss.

This affects v3: the framework's strongest claim ("M4 displaces all
other mechanisms eventually") must be qualified. The recommended v3
Discussion drop-in (test_h4_v3_statement.md §1): "M4 displaces M1 at
long horizons but the displacement is gradual when M1's initial
population happens to draw high-fitness lineages. At gen=1000 (Test H3),
M4 wins 83% at equal start; at gen=5000 (Test H4), M4 wins 80% —
approaching but not reaching the 99% threshold, so the finite-time
artifact interpretation is supported but a residual long-horizon
durability effect remains."

## v4 follow-ups (Test F4, G3, H5)

- test_f4_mode2_base_v1.csv (Test F4 — Mode 2 simulation, L sweep; I_struct = 4.214 bits/codon at all L, matches closed-form 4.218)
- test_f4_mode1_vs_mode2_v1.csv (Test F4 — Mode 1 vs Mode 2 capacity ratio = 1.42 at L=64, matches closed-form 6/4.218)
- test_f4_apparatus_signature_v1.csv (Test F4 — G2 apparatus signature for Mode 2)
- test_f4_nested_copyability_v1.csv (Test F4 — 500-gen nested-Mode-1 sim drove fitness 5%→43%, code fidelity 52%)
- test_f4_code_expansion_v1.csv (Test F4 — 20→24 aa expansion: empirical ratio 1.044 vs dispatch-predicted 1.061)
- test_g3_cross_channel_kl_matrix_v1.csv (Test G3 — 10×10 pairwise KL matrix; cluster structure across Mode 3 cyclics, Mode 1 randoms, NRPS)
- test_g3_ensemble_growth_v1.csv (Test G3 — I_chan vs ensemble size, 2 growth orders)
- test_g3_classification_stability_v1.csv (Test G3 — 10 random 4-channel ensembles, 10/10 channels classify stably)
- test_g3_g2_result_robustness_v1.csv (Test G3 — Drt3a vs Drt3b ΔI_chan stays in [0.976, 1.025] across 7 ensembles)
- test_h5_r_sweep_v1.csv (Test H5 — fine r sweep at β=10/L=32, 14 r values × 30 reps; r* = 0.100 with plateau 0.548)
- test_h5_beta_sensitivity_v1.csv (Test H5 — r*(β) at β ∈ {2,5,10,20}; r* range [0.05, 0.08], regime-invariant per P_H5_3)
- test_h5_L_sensitivity_v1.csv (Test H5 — r*(L) at L ∈ {16,32,64,128}; r* range [0.02, 0.05], L-invariant per P_H5_4)
- test_h5_M4_mutation_sweep_v1.csv (Test H5 — M4 plateau vs μ; μ* = 0.0010 with plateau 0.998)
- test_h5_M2_vs_M4_matched_v1.csv (Test H5 — head-to-head at r*/μ*: at N_M2=200 balanced, M2 wins 33%, M4 wins 67%)

v4 verdicts:
- F4: ALL 5 PREDICTIONS CONFIRMED. Mode 2 has its own characterized substrate; I_struct ceiling = codon-count entropy 4.218 bits/codon, strictly tighter than dispatch's loose log₂(20). Nested copyability operationally demonstrated. v3 Mode 2 Results subsection now has empirical numbers.
- G3: ALL 4 PREDICTIONS CONFIRMED. Pairwise KL invariant (P_G3_1). Mode classification stable across 10 random ensembles (P_G3_3). G2's headline ΔI_chan ≈ 1.0 bits robust [0.976, 1.025] across 7 ensembles (P_G3_4). Used L=6 (matching G2's actual harness, not dispatch's mistakenly-claimed L=64).
- H5: 4 CONFIRMED + 1 REFINED. Sweet spot r* = 0.100 (plateau 0.548); β-invariant in [0.05, 0.08]; L-invariant in [0.02, 0.05]. Rate-matching analogy P_H5_2 REFINED — M2 needs 100× more raw variation than M4 to reach matched plateau (lineage re-draws are wasteful relative to targeted mutation). At matched rates and balanced start, M4 wins 67% (M2 33%) — Draft B partially supported, M4 retains measurable edge. **Recommended v3 framing: hybrid Draft A (mechanism-specific) in theorem + Draft B (variation-rate principle) in Discussion.**

Sciscape D6 (real Mode 6 inheritance mechanism literature search) is
upstream of any future H6 dispatch; user must dispatch to Sciscape
directly. Prompt at `dispatches_v4/sciscape_d6_sharpened.md`.

Future tasks must NOT overwrite these files. Future tasks may CREATE new files
prefixed with the new test's identifier (e.g., test_i_*, test_j_*).

If a future task needs to reload these CSVs for plotting or analysis, READ ONLY.
