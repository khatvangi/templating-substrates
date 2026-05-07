# Orchestration: Test F2 + G2 + H2 dispatch plan

## Run order

The three dispatches are independent at the simulation level. Order is
determined by wall-time and by which results unblock which downstream
v3 writing tasks.

1. **Test G2** (apparatus repair) — dispatch first as a regular Boron task.
   Wall-time: <5 min. Unblocks Methods rewrite.

2. **Test F2** (robustness sweep) — dispatch in parallel with G2 as a
   regular task. Wall-time: ~30 min. Unblocks Test F prediction
   writeup.

3. **Test H2** (founder boundary) — dispatch as nohup background task.
   Wall-time: ~9 hours. Unblocks Discussion rewrite.

Strict parallelism: G2 and F2 in foreground, H2 in background. H2
should be started first (so its 9-hour clock starts) but reported on
last.

## H2 nohup pattern

H2 cannot fit in the 15-minute wall-time budget. Dispatch it as:

```
# (Claude Code on Boron)
cd /storage/kiran-stuff/templating_framework
mkdir -p results/test_h2_progress
nohup python code/test_h2_competition_sweep.py \
    --progress-file results/test_h2_progress/progress.txt \
    --completed-csv results/test_h2_progress/completed.csv \
    > results/test_h2_progress/log.txt 2>&1 &
echo $! > results/test_h2_progress/pid.txt
```

The simulation script must:
- Write progress.txt every 50 completed sims with timestamp, n_complete,
  estimated time remaining.
- Append to completed.csv every sim so partial results are usable on
  early termination.
- Catch SIGTERM cleanly; write a summary to progress.txt on exit.

Polling pattern (from the next claude session):
- Read progress.txt every ~30 minutes.
- If no progress for >2 hours, check the log for errors.
- When n_complete = 2200 (or whatever the total is), proceed to Step 4
  (predictions evaluation).

## Result reconciliation

Once all three dispatches finish, the results need to be reconciled into
a v3 plan:

1. **G2 outputs** drive the v3 Methods section apparatus repair. The
   `test_g2_v3_methods_paragraph.md` is the seed.

2. **F2 outputs** drive the v3 Results subsection on family predictions.
   The robust envelope replaces the v1 point estimate. The hypothetical
   universal-gate predictions become a "predictions to be tested by SDM"
   subsection.

3. **H2 outputs** drive the v3 Discussion paragraph on inheritance
   carrier competition. The founder-loss boundary characterization
   replaces v2's published claim with a properly bounded version.

I should hold off on writing v3 until all three finish — partial
reconciliation creates churn. But I can start on the v3 outline now,
leaving placeholders for the F2/G2/H2 inserts.

## Sanity check: what could go wrong

- H2 takes longer than 9 hours. Mitigation: progress.txt polling; if
  N_REPLICATES needs to drop from 30 to 20 to fit budget, that's a
  documented compromise.
- G2 finds that *neither* observable distinguishes Drt3a from Drt3b
  (both give ~0 bits or both give ~12 bits at all conditions). This
  is genuinely possible; the v3 apparatus repair would then need a
  third observable (per-base fidelity as a non-MI measurement). If
  this happens, flag it; don't try to force a resolution.
- F2 finds the envelope is fragile (E26D sensitivity is sharp; ε
  tolerance is tight). This weakens the falsifiability claim and
  needs honest reporting in v3 — the framework's prediction is "narrow
  but valid" rather than "robust."
- All three find their predictions confirmed cleanly. Ideal outcome;
  v3 writes itself.

## File structure on Boron after all three finish

```
/storage/kiran-stuff/templating_framework/
├── code/
│   ├── test_f2_robustness.py
│   ├── test_g2_dual_observable.py
│   └── test_h2_competition_sweep.py
├── results/
│   ├── test_f2_*_v1.csv (5 files)
│   ├── test_f2_robust_falsifiability_v1.md
│   ├── test_g2_*_v1.csv (3 files)
│   ├── test_g2_apparatus_decision_v1.md
│   ├── test_g2_v3_methods_paragraph.md
│   ├── test_h2_*_v1.csv (3 files)
│   ├── test_h2_predictions_v1.md
│   ├── test_h2_v3_discussion_paragraph.md
│   ├── test_h2_progress/ (progress files)
│   └── CANONICAL_RESULTS.md (updated)
├── figures/
│   ├── test_f2_*.png (3 files)
│   ├── test_g2_*.png (2 files)
│   └── test_h2_*.png (3 files)
└── HISTORY.md (entries for F2, G2, H2)
```
