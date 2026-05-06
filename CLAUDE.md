# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Computational validation suite for the **Templating Substrates Framework** — an
information-theoretic theory that classifies molecular templating systems by
how their structural mutual information `I_struct(X; Y)` scales with output
length `L`. The signature distinction is:

| Mode | Channel                                | I_struct(X; Y) scaling             |
|------|----------------------------------------|------------------------------------|
| 1    | sequence template (DNA-replication)    | linear in `L`                      |
| 3    | cyclic active site (Drt3b, N states)   | saturates at `log_2(N)`, L-indep.  |
| 4    | single conformer                       | bounded                             |
| 5    | conveyor / N modules                   | bounded by templated module count  |
| 6    | 2D surface, lossy inheritance          | bounded                             |

Tests A.1, A.2, B, C, D, E each validate one prediction. Test E is the
biological anchor against the **Drt3 system** (Sharma et al. 2026, *Science* —
supplementary data lives in `science.aed1656_data_s1_to_s8 deng/`).

## Directory layout

```
code/                         one self-contained Python script per test
results/                      CSVs, READMEs, progress.txt, stdout/stderr.log, pid.txt
figures/                      PNG plots written by the test scripts
paper/drafts/                 (currently empty) manuscript drafts
science.aed1656_data_s1_to_s8 deng/   Sharma et al. 2026 Drt3 supp. data
HISTORY.md                    session Q&A log (append-only, see global CLAUDE.md)
recover_results.py            one-off utility: reconstruct CSVs from progress logs
results/CANONICAL_RESULTS.md  sentinel listing CSVs that must NOT be overwritten
```

## Running tests

Each script is independently runnable from the repo root with stdlib + numpy + matplotlib:

```bash
python code/test_a1_mode1_scaling.py        # fast (~seconds)
python code/test_a2_bulk_matched_control.py # fast
python code/test_b_mode3_capacity.py        # ~minutes — long sweeps
python code/test_c_mode5_conveyor.py        # ~minutes
python code/test_d_v2_population_dynamics.py  # slow — 200 agents × 1000 generations × 5 modes
python code/test_e_v2_drt3_classification.py
```

Long-running tests (B, C, D, E) follow a background-execution pattern:
`results/test_<id>_pid.txt`, `progress.txt`, `stdout.log`, `stderr.log`. When
running them as background jobs, write to those same files so the existing
recovery / monitoring scripts keep working.

A test prints a PASS/FAIL report at the end and writes:
- `results/test_<id>_results.csv` — the data
- `results/test_<id>_README.md` — what/why/PASS criteria/interpretation
- `figures/test_<id>_<plotname>.png` — diagnostic plots

Reproducibility: every test seeds with `np.random.seed(42)` and
`np.random.default_rng(42)`. Re-running reproduces tables exactly.

## Hard rules specific to this repo

### 1. CANONICAL_RESULTS.md is a sentinel — do not overwrite listed CSVs

`results/CANONICAL_RESULTS.md` lists CSVs that are *frozen passing results*.
Future work MUST NOT overwrite them. New work must prefix all artifacts with
its own test ID (`test_d_*`, `test_e_*`, …). Re-analysis of a canonical CSV
must be **read-only**.

### 2. Re-runs use `_v2` suffix, never overwrite

When a test is re-run with a corrected setup, write `test_<id>_v2_*` rather
than overwriting `test_<id>_*`. Existing examples: `test_d_v2_*`,
`test_e_v2_*`. This preserves the historical audit trail of every verdict.

### 3. Verdict corrections are pure CSV re-analysis

If only the PASS criterion was wrong (data fine), write a small
`test_<id>_verdict_corrected.py` that uses the `csv` module only, reads the
canonical CSV, applies the corrected criterion, and writes
`results/test_<id>_verdict_corrected.md`. Do NOT re-simulate. Existing
examples: `code/test_e_verdict_corrected.py`, `results/reeval_test_b.py`.

### 4. Test scripts are self-contained on purpose

Each test duplicates its own simulator and estimator rather than importing
from sibling tests. This is deliberate (`test_a2`'s README explicitly notes:
"A.2 cannot drift if A.1 changes"). Do **not** refactor common code into a
shared `utils.py` — the duplication is the point. If you fix a simulator bug
in test X, the equivalent fix must be applied to test Y separately and
intentionally.

### 5. HISTORY.md is append-only

Per the global `~/.claude/CLAUDE.md` rule, append a `## YYYY-MM-DD` block
with `**Q:**` / `**A:**` pairs each session. Don't record code or detailed
work — just the Q&A summary.

## Test script anatomy (the convention to follow)

Every `code/test_*.py` follows the same shape; new tests should mirror it:

1. **Long module docstring** stating: framework section, channel
   parameterization, closed-form theoretical reference, PASS criterion.
2. **Path setup**: `SCRIPT_DIR = Path(__file__).resolve().parent`,
   `PROJECT_ROOT = SCRIPT_DIR.parent`, then `RESULTS_DIR` / `FIGURES_DIR`.
3. **Headless matplotlib**: `matplotlib.use("Agg")` before `import pyplot`.
4. **Simulator** (vectorized over `n_samples` and `L`).
5. **Estimator**: plug-in MI from empirical joints. Per-position vs joint
   choice depends on whether positions are independent given X (Mode 1 yes,
   Mode 3 no — see `test_b_README.md` for the rationale).
6. **Theoretical reference** (closed form when possible).
7. **Sweep** across parameters with seeded RNG.
8. **Write CSV → write figures → write README → print PASS/FAIL report.**

## Visualization caveat

Existing tests use raw `matplotlib.pyplot`. The user's global preference
(`~/.claude/CLAUDE.md`) is **seaborn for static figures, plotly for
interactive**. New plotting code should follow the global preference;
existing matplotlib code should be migrated when touched, but don't refactor
purely-cosmetic plotting.

## Style

- never capitalize comments
- function names like `findMeaning`, not `getMeaning`
- discuss approach before writing code; don't rush to implement
- never fabricate data or "example" outputs — if a sim hasn't run, say so
