# Orchestration: Five-fork parallel exploration plan

## What's being run

Six parallel tasks, dispatched together to Boron and Sciscape:

| Task | Type | Wall-time | Dispatch mode |
|------|------|-----------|---------------|
| H4 long-horizon M4 vs M1 | Defensive | ~2 hours | nohup |
| H5 M2 sweet spot characterization | Generative | ~1.5 hours | nohup |
| F4 Mode 2 simulation | Generative | ~45 min | foreground or nohup |
| G3 ensemble robustness | Defensive | ~10 min | foreground |
| Sciscape D6-sharpened | Literature | async | Sciscape |

H4 and H5 are the long poles. If dispatched in parallel they share
Boron compute resources — confirm Claude Code on Boron can handle
two nohup processes simultaneously. If not, queue: H5 first (more
generative), H4 second.

F4 and G3 are short. Run F4 in foreground after the nohup tasks
launch; run G3 in foreground after F4 finishes.

## Dispatch sequence

```
# (Claude Code on Boron)
cd /storage/kiran-stuff/templating_framework

# 1. Launch H5 nohup (most generative, longest)
mkdir -p results/test_h5_progress
nohup python code/test_h5_sweet_spot.py \
    --progress-file results/test_h5_progress/progress.txt \
    --completed-csv results/test_h5_progress/completed.csv \
    > results/test_h5_progress/log.txt 2>&1 &
echo $! > results/test_h5_progress/pid.txt

# 2. Launch H4 nohup (defensive, comparable length)
mkdir -p results/test_h4_progress
nohup python code/test_h4_long_horizon.py \
    --progress-file results/test_h4_progress/progress.txt \
    --completed-csv results/test_h4_progress/completed.csv \
    > results/test_h4_progress/log.txt 2>&1 &
echo $! > results/test_h4_progress/pid.txt

# 3. F4 in foreground (45 min)
python code/test_f4_mode2_simulation.py
# When F4 finishes:

# 4. G3 in foreground (10 min)
python code/test_g3_ensemble_robustness.py
```

Sciscape dispatch is independent of Boron — paste the D6-sharpened
prompt into Sciscape directly. Results land async (typically <30 min
for a single Sciscape query).

## Polling and status

Every 30 minutes during the long-running phase, poll H4 and H5 progress:

```
cat results/test_h4_progress/progress.txt
cat results/test_h5_progress/progress.txt
```

If either shows no progress for >2 hours, check the log:

```
tail -50 results/test_h4_progress/log.txt
tail -50 results/test_h5_progress/log.txt
```

Expected timeline:
- T+0: H5 + H4 launched, F4 running in foreground
- T+45min: F4 done, G3 starts
- T+55min: G3 done
- T+90min: H5 should be done (if 1.5h estimate holds)
- T+120min: H4 should be done (if 2h estimate holds)
- Total: ~2 hours wall-time on Boron

Sciscape D6 returns asynchronously, probably within an hour.

## What each result unlocks for v3

**G3 unlocks**: A robustness paragraph in v3 Methods next to the G2
apparatus definition. Reviewer-resistant.

**F4 unlocks**: A genuine Mode 2 Results subsection with empirical
information-content numbers, the apparatus signature applied to
translation, and (if Step 4 succeeds) a nested-copyability
demonstration. Mode 2 stops being "Mode 1 plus a lookup" and becomes
its own characterized substrate.

**H4 unlocks**: A finite-time-vs-asymptotic clarification paragraph
in v3 Discussion that defends the Draft A theorem statement against
the P_H3_6 falsification. If H4 confirms M4 wins ≥99% at gen=5000,
the theorem is empirically defended at long horizons.

**H5 unlocks**: The strongest possible support (or refutation) for
the Draft B variation-rate principle. If r* is regime-invariant and
matches μ* at matched rates, Draft B is empirically grounded and v3's
hybrid framing is cleanly defended. If r* shifts with regime, Draft B
is conjectural and v3 should weaken its Discussion framing.

**Sciscape D6 unlocks**: H6 (real Mode 6 simulation), which can't be
written until literature establishes which mechanism class real Mode 6
systems implement. Also feeds v3's Mode 6 paragraph directly: if the
literature shows real Mode 6 = M0, v3 strengthens its claim; if real
Mode 6 = M1, v3 weakens to "Mode 6 has weak inheritance."

## Decision tree after all five finish

When all results land, the v3 writing decisions break down as:

1. **G3 confirms ensemble invariance (P_G3_3)**: v3 Methods adds a
   short note. Continue.

2. **G3 falsifies invariance**: v3 Methods specifies a canonical
   ensemble and reports sensitivity range. Continue.

3. **F4 confirms Mode 2 information bound (P_F4_1) and apparatus
   signature (P_F4_3)**: v3 Results gets a Mode 2 subsection with
   empirical numbers. Continue.

4. **F4 surprises**: depending on which prediction fails, v3 Results
   needs reframing of Mode 2 — possibly Mode 2 reduces to Mode 1
   (taxonomy shrinks), possibly Mode 2 has unexpected properties
   (taxonomy grows).

5. **H4 confirms M4 win rate ≥99% at gen=5000**: v3 Discussion has a
   finite-time clarification footnote, theorem statement holds.

6. **H4 surprises (M1 has structural durability)**: v3 Discussion
   has to acknowledge bounded-horizon dominance, not asymptotic.

7. **H5 supports Draft B (rate-matching wins)**: v3 hybrid framing
   gets empirical teeth. Discussion paragraph on M2_r* sweet spot
   becomes a Results subsection.

8. **H5 supports Draft A (M4 wins even at matched rates)**: v3 leans
   harder on Draft A; M2_r=0.10 sweet spot becomes a footnote
   rather than a principle.

9. **Sciscape D6 returns M0 evidence for biological systems**: v3
   Mode 6 paragraph is unambiguous.

10. **Sciscape D6 returns M1 evidence**: v3 Mode 6 paragraph weakens
    to "Mode 6 has lineage-level inheritance but not M4-class
    inheritance."

11. **Sciscape D6 returns nothing definitive**: v3 flags the
    inheritance question as open and a target for future work.

The most important conjunction: H5 supports Draft B AND Sciscape D6
returns M0 for real systems → v3 has the strongest possible
inheritance theorem framing. Every other combination has caveats but
is still publishable.

## v3 writing readiness after all five finish

v3 Methods: ready (G2 + G3 outputs).
v3 Results: ready (F2 + F3 + F4 outputs).
v3 Theorem: ready (Draft A from H3, robustness-tested by H4).
v3 Discussion: ready (Draft B framing from H3 + H5, M2 sweet spot
  paragraph, Mode 6 reframing from Sciscape D6).

Stitching task remaining:
1. Read all v3-statement and v3-paragraph documents
2. Integrate with existing v2 prose
3. Resolve any v2-vs-v3 contradictions
4. Pre-submission review

## File structure on Boron after all five finish

```
/storage/kiran-stuff/templating_framework/
├── code/
│   ├── test_g3_ensemble_robustness.py
│   ├── test_f4_mode2_simulation.py
│   ├── test_h4_long_horizon.py
│   └── test_h5_sweet_spot.py
├── results/
│   ├── test_g3_*_v1.csv (4 files) + v3_methods_note.md
│   ├── test_f4_*_v1.csv (5 files) + v3_statement.md
│   ├── test_h4_*_v1.csv (3 files) + v3_statement.md
│   ├── test_h5_*_v1.csv (5 files) + v3_statement.md
│   ├── test_h4_progress/, test_h5_progress/
│   ├── lit_search/D6_sharpened_<date>.md
│   └── CANONICAL_RESULTS.md (updated with all five)
├── figures/
│   └── ~12 new PNGs across F4, G3, H4, H5
└── HISTORY.md (5-test consolidated entry)
```

## Risk: simultaneous nohup processes

If H4 and H5 launching together exceed Boron compute capacity
(memory, CPU), one will slow or fail. Mitigation:

- Launch H5 first, wait 60 seconds, check it's running normally
  (progress.txt is being written), then launch H4.
- If both processes show normal progress, proceed.
- If H5 stalls or memory warnings appear, kill H4 and run sequentially.

Both H4 and H5 are pure-Python population-dynamics simulations with
modest memory footprints (K=400 agents, L≤128). Should be fine on
typical Boron setup but worth verifying on launch.
