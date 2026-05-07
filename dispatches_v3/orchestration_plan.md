# Orchestration: Test H3 + F3 dispatch plan

## Run order

H3 is the bigger of the two and requires nohup. F3 is short and can
run in foreground. Both should be dispatched promptly so H3's clock
starts.

1. **Test H3** (inheritance landscape) — dispatch first as nohup
   background. Wall-time: ~8 hours. Unblocks v3 inheritance theorem
   reformulation.

2. **Test F3** (gate classification) — dispatch as regular task.
   Wall-time: ~15 min. Unblocks v3 universal-gate prediction
   refinement.

## H3 nohup pattern (same as H2)

```
cd /storage/kiran-stuff/templating_framework
mkdir -p results/test_h3_progress
nohup python code/test_h3_inheritance_landscape.py \
    --progress-file results/test_h3_progress/progress.txt \
    --completed-csv results/test_h3_progress/completed.csv \
    > results/test_h3_progress/log.txt 2>&1 &
echo $! > results/test_h3_progress/pid.txt
```

The simulation script must:
- Write progress.txt every 50 sims with timestamp, n_complete, ETA.
- Append to completed.csv each sim so partial results are usable on
  early termination.
- Catch SIGTERM cleanly; write summary to progress.txt on exit.

Polling: read progress.txt every ~30 minutes. If no progress for >2
hours, check log.txt for errors. When n_complete = 2000 (or whatever
total is), proceed to Step 5 (v3 conceptual reformulation).

## What to do while H3 is running

The two productive things to work on in parallel:

1. **F3 results.** Once F3 finishes (~15 min from dispatch), the
   architecture-vs-selectivity distinction can be folded into v3
   Results. The SDM-ready prediction table is a clean v3 addition.

2. **v3 outline construction.** With G2's Methods paragraph in hand,
   F2's robust falsifiability statement, F3's SDM predictions, and
   the H3 result pending, the v3 outline can be drafted with
   placeholders for the H3 Discussion revision. Outline tasks:
   - Methods: drop in G2 paragraph, update controls subsection
   - Results: extend with F2's robust envelope and F3's gate
     classification
   - Discussion: leave placeholder for H3-driven inheritance
     theorem revision

3. **What NOT to do while H3 is running:** do not write the v3
   Discussion paragraph on Mode 6 / inheritance carrier without
   H3 results. The decision among Drafts A, B, C in
   `test_h3_v3_inheritance_revision.md` requires H3's data, not
   pre-judgment.

## Result reconciliation after both finish

Once F3 and H3 both finish, the v3 picture is:

- **G2 outputs** → v3 Methods (drop-in apparatus paragraph).
- **F2 outputs** → v3 Results subsection on family predictions
  (robust envelope replacing v1 point estimate).
- **F3 outputs** → v3 Results subsection on universal-gate
  architecture-vs-selectivity distinction (refining F2).
- **H3 outputs** → v3 Discussion (inheritance theorem revision,
  Mode 6 reframing, response to the original v2 line 127 claim).

The v3 Methods, Results, and Discussion now have explicit drop-in
paragraphs from each test's output document. The v3 writing is
largely a stitching task — read the drop-ins, integrate them with
v2's existing prose, mark places where v2 prose contradicts v3
results.

## Sanity checks: what could go wrong

- **H3 takes longer than 8 hours.** Mitigation: progress.txt polling.
  If N_REPLICATES needs to drop from 30 to 20 to fit budget, that's
  documented.
- **H3 finds P_H3_10 fails** (some non-M4 mechanism reaches > 0.95
  of max). This would be a deeper challenge to the framework than
  H2's finding. Flag immediately; don't try to interpret. The framework
  would need substantial revision, not just the inheritance condition
  rewrite.
- **F3 finds R253A doesn't actually disrupt the cycle** (P_F3_1 fails,
  I_struct > 0.5 bits). This would mean the framework's identification
  of R253 as architectural is wrong. The SDM prediction set then needs
  redrawing.
- **F3 finds G248A breaks the cycle** (P_F3_2 fails, I_struct < 0.5
  bits). This would mean architectural and selectivity gates can't be
  cleanly separated; the F3 distinction would not hold and v3 should
  treat all primary gates as architectural.

## File structure on Boron after both finish

```
/storage/kiran-stuff/templating_framework/
├── code/
│   ├── test_f3_gate_substitution_sweep.py
│   └── test_h3_inheritance_landscape.py
├── results/
│   ├── test_f3_*_v1.csv (4 files)
│   ├── test_f3_sdm_predictions_v1.md
│   ├── test_f3_v3_results_paragraph.md
│   ├── test_h3_*_v1.csv (4 files)
│   ├── test_h3_v3_inheritance_revision.md
│   ├── test_h3_progress/
│   └── CANONICAL_RESULTS.md (updated)
├── figures/
│   ├── test_f3_*.png (2 files)
│   └── test_h3_*.png (4 files)
└── HISTORY.md (entries for F3, H3)
```

## Reading queue when results land

When F3 finishes, read in order:
1. `test_f3_classification_v1.csv` — does the architecture/selectivity
   distinction hold?
2. `test_f3_sdm_predictions_v1.md` — the v3-ready table.

When H3 finishes, read in order:
1. `test_h3_isolated_dynamics_v1.csv` — per-mechanism plateau heights.
   Verify P_H3_10 across β regimes.
2. `test_h3_pairwise_competition_v1.csv` — does M4 still win pairwise?
3. `test_h3_v3_inheritance_revision.md` — the Drafts A/B/C analysis
   for the revised condition (i).
