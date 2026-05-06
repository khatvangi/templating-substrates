# Templating Substrates Framework

A computational validation suite for an information-theoretic theory of
molecular templating. The framework classifies templating systems by
how their structural mutual information `I_struct(X; Y)` scales with
output length `L`:

| Mode | Channel                                | I_struct(X; Y) scaling             |
|------|----------------------------------------|------------------------------------|
| 1    | sequence template (DNA-replication)    | linear in `L`                      |
| 2    | code-encoded cross-descriptor (translation) | bounded by codon-count entropy |
| 3    | cyclic active site (Drt3b, N states)   | saturates at `log₂(N)`, L-indep.   |
| 4    | single conformer                       | bounded                            |
| 5    | conveyor / N modules                   | bounded by templated module count  |
| 6    | 2D surface, lossy inheritance          | bounded                            |

The biological anchor is the **Drt3 system** (Sharma et al. 2026,
*Science* — DOI [10.1126/science.aed1656](https://doi.org/10.1126/science.aed1656)),
a recently discovered bacterial anti-phage reverse transcriptase that
synthesises alternating poly(AC) DNA without a nucleic acid template.

## Repository layout

```
code/                         self-contained Python script per test
code/figures_v7/              publication figure scripts (one per panel)
results/                      CSVs, READMEs, per-test verdicts
results/CANONICAL_RESULTS.md  sentinel listing CSVs that must NOT be overwritten
figures/                      diagnostic PNGs written by test scripts
paper/figures/v7/             publication-grade panels (PDF + PNG, ready for Inkscape)
templating_substrates_draft_v2.md   manuscript draft
HISTORY.md                    chronological Q&A log of the project
data/                         pointer to third-party data (not redistributed)
recover_results.py            utility: reconstruct CSVs from progress logs
```

## Reproducing the canonical results

The repository ships with the canonical CSVs already generated. Each test
script seeds with `np.random.seed(42)` and `np.random.default_rng(42)`,
so re-running reproduces the tables exactly.

### Dependencies

Stdlib + a few standard scientific Python packages. Tested on Python 3.10+.

```bash
pip install numpy matplotlib seaborn pillow
```

(Some publication figure scripts also use `subprocess` calls to `pdfinfo`
for validation; install via `apt install poppler-utils` on Debian/Ubuntu
or `brew install poppler` on macOS.)

### Quick checks (seconds)

```bash
python code/test_a1_mode1_scaling.py        # Mode 1 length-scaling
python code/test_a2_bulk_matched_control.py # bulk-matched control
```

Both finish in seconds and write CSVs that should match the checked-in
`results/test_a1_results.csv` and `results/test_a2_results.csv` exactly.

### Long-running tests (minutes)

```bash
python code/test_b_mode3_capacity.py        # Mode 3 capacity — minutes
python code/test_c_mode5_conveyor.py        # Mode 5 conveyor   — minutes
python code/test_d_v2_population_dynamics.py  # population dynamics — slow
python code/test_e_v2_drt3_classification.py  # Drt3 anchoring (requires data/)
```

### Family-level tests (require third-party data)

Tests F, F2, F3 read the 1,232-sequence Drt3b alignment from
`science.aed1656_data_s3.fa`. **Download instructions in
[`data/README.md`](data/README.md).** Without the data, these tests
cannot run; the rest of the framework (A.1, A.2, B, C, D, E, G, H, H4,
H5) is fully self-contained.

### Long-horizon population sweeps (background)

Tests H2, H3, H4, H5 are large parameter sweeps written for nohup
dispatch. Their runtime ranges from a few minutes (vectorised at
K=400 agents) to ~30 minutes for H5's 2,130-cell sweep. Invoke pattern:

```bash
mkdir -p results/test_h5_progress
nohup python code/test_h5_sweet_spot.py \
    --progress-file results/test_h5_progress/progress.txt \
    --completed-csv results/test_h5_progress/completed.csv \
    > results/test_h5_progress/log.txt 2>&1 &
```

### Publication figures

```bash
cd code/figures_v7
python fig5_panel_a_m0_m4_plateaus.py   # creates the style module first time
python fig1_panel_a_apparatus_calibration.py
# ... 17 more panel scripts
```

Each script writes one PDF (vector, editable text) + one PNG (600 dpi)
to `paper/figures/v7/figN/figN_panel_X.{pdf,png}`. Composite assembly
into multi-panel figures happens in Inkscape, not Python.

## Hard rules (project conventions)

These are documented in detail in [`CLAUDE.md`](CLAUDE.md). Summary:

1. **`results/CANONICAL_RESULTS.md` is a sentinel.** The CSVs it lists
   are frozen passing results. Future work must not overwrite them; new
   work must prefix all artifacts with its own test ID.
2. **Re-runs use a `_v2` (or `_v3`, `_v4`) suffix, never overwrite.**
3. **Verdict corrections are pure CSV re-analysis.** If only the PASS
   criterion was wrong (data fine), use the `csv` module only and write
   `test_<id>_verdict_corrected.md`. Do not re-simulate.
4. **Test scripts are self-contained on purpose.** Each script
   duplicates its own simulator and estimator. Do not refactor common
   code into a shared `utils.py` — the duplication is the point. If a
   simulator bug is fixed in test X, the equivalent fix must be applied
   to test Y separately and intentionally.
5. **`HISTORY.md` is append-only.** Each session adds a `## YYYY-MM-DD`
   block with `**Q:**` / `**A:**` pairs.

## Manuscript

The active manuscript draft is
[`templating_substrates_draft_v2.md`](templating_substrates_draft_v2.md)
(targeting PRX Life or eLife). The v3 round of work (Tests F2/G2/H2/H3/H4/F3/H5/F4/G3)
produced drop-in paragraphs for the Methods, Results, and Discussion
sections — see the `test_*_v3_statement.md` and `test_*_v3_methods_*.md`
files in `results/`.

## Citing

```
Boggavarapu Kiran (2026). Templating Substrates Framework: a
computational validation suite.
https://github.com/khatvangi/templating-substrates
```

When the manuscript is published, prefer the journal citation.

## License

- Code (`code/`, `code/figures_v7/`, `recover_results.py`):
  [MIT](LICENSE)
- Data, figures, manuscript: [CC-BY-4.0](LICENSE-data.md)
- Third-party Drt3 data (not redistributed): governed by the original
  publisher; see [`data/README.md`](data/README.md).
