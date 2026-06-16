# Supplemental Material — portal field copy

Paste-ready text for the two PRX Life / bioRxiv portal fields about
supplemental material. Each section corresponds to a separate field.

---

## Field 1: Description of Supplemental Material

### Uploaded with this submission

**`supplementary-v3.pdf`** (17 pages) — Mathematical Appendix providing formal definitions and proofs for the apparatus and theorems stated in the main paper. Five sections:

1. **Formal apparatus definitions** — defines the templating event, population mutual information `I_struct^pop`, substrate capacity `C_S`, channel mutual information `I_struct^chan`, per-base mechanism fidelity (phase fidelity `f_phase` and monomer fidelity `f_monomer`), and the joint apparatus signature. Includes the descriptor-relative formulation and the fixed-template degeneracy corollary used to justify the dual `I_pop` / `I_chan` apparatus in the main paper.

2. **Mode-specific information bounds** — theorems and proofs for each mode's capacity scaling: Mode 1 linear (`L · [log₂|A| − h_b(ε) − ε log₂(|A|−1)]`), Mode 3 logarithmic saturation (`I ≤ log₂ N`), Mode 5 module-bounded (`I ≤ min(L,N) · [...]`), Mode 2 code-mediated ceiling (`L · H(ν) = 4.218L` bits/codon for the standard genetic code), Mode 4 length-independence, and the Mode 6 reframing. Each theorem cross-references its empirical validation in the main paper figures.

3. **Descriptor-relativity** — formal statement that mode classification applies to (templating event, descriptor) pairs, not to physical objects. Includes the designed β-sheet peptide-replicator example demonstrating that the same template can classify into different modes under different descriptors, and a remark explaining why amyloid sequence-templating is *not* the right example for this point.

4. **Bulk-matched controls** — formalization of the *k*-th order matched null and the paired-null bulk-control construction used as the operational null hypothesis throughout the main paper.

5. **Capacity and generation theorems** — the two necessity results separating substrate capacity from cumulative generation: the capacity theorem (R1, R2, R3 necessary for unlimited heredity capacity `C_S(L) = Ω(L)`, with input ensembles restricted to the recursively-closed set `R_L`) and the generation theorem (R4 additionally necessary for descendant sets to strictly contain the zero-drift closure `cl_0(S_0)`). Includes the bounded-heredity corollary mapping each mode to its specific R-condition failure, the finite-population finite-horizon proposition that scopes the M0–M4 plateau claim, and a synthesis table mapping ten prior formulations (Schrödinger, Szathmáry, Eigen, Wagner, Pross, Hull, Maynard Smith, Hofmeyr, Vasas, Walker–Davies) to the corresponding R-conditions.

A DOCX version (`supplementary-v3.docx`) with the same content is also available on request and at the deposit DOI below.

### Linked via Data Availability (archived on Zenodo)

All simulation code, raw simulation outputs, the diagnostic apparatus implementation, the family-level analysis pipeline, the pre-specification dispatches, and the figure-generation pipeline are deposited at:

- **`10.5281/zenodo.20060972`** — pre-specification anchor (v0.0-prereg release): the analysis-plan dispatch markdown files at the moment of pre-registration, before empirical work began.
- **`10.5281/zenodo.20272479`** — reproducibility snapshot (v1.0-submission release): the full submission state, containing:
  - **23 Python test scripts** (`code/test_*.py`) — one per apparatus test (A1, A2, B, C, E, F, F2, F3, F4, G, G2, G3, H, H2, H3, H4, H5, plus the v2 re-runs); each is self-contained (deliberate anti-DRY) and reproduces its result CSV under seeded RNG.
  - **57 result CSV files** (`results/*.csv`, ~7.6 MB) — the empirical data underlying every numerical claim in the main paper: per-test sweeps, per-cell summaries, per-clade signatures, per-replicate competition outcomes.
  - **6 figure-dispatch markdown files** (`dispatch_fig1...` through `dispatch_fig6_...`) plus the panel-generation pipeline at `code/figures_v7/` reproducing all 19 manuscript panels.
  - **Phylogenetic input** at `data/` — derived from the published supplementary material of Sharma/Deng et al. (Science 2026, DOI 10.1126/science.aed1656); the third-party data files themselves are not redistributed.
  - **Bibliography** `templating_substrates.bib` and the manuscript source `manuscript-v16.tex` with all figure assets at `paper/figures/v16/`.

The deposit's git history records the commit ordering: Phase 0 retroactively commits the dispatch files with their original filesystem modification times preserved as commit metadata; Phases 1–5 commit the test scripts, result CSVs, and figure pipeline. Pre-specified predictions are reproduced verbatim in the docstring headers of each test script (e.g., `code/test_h_competition.py` lines 18–38).

### Short version (for portals with strict character limits)

> **`supplementary-v3.pdf`** — Mathematical Appendix (17 pp) with formal definitions, mode-specific capacity proofs, descriptor-relativity, bulk-control formalism, and the capacity/generation theorems with bounded-heredity corollary and finite-population proposition. All simulation code, raw outputs, and pre-specification dispatches are deposited on Zenodo at 10.5281/zenodo.20060972 (pre-spec anchor) and 10.5281/zenodo.20272479 (submission snapshot).

---

## Field 2: Special Instructions for Accessing the Files

### Files included with this submission

**`supplementary-v3.pdf`** (Mathematical Appendix) opens in any standard PDF viewer; no special software, plugin, or authentication required. A DOCX companion (`supplementary-v3.docx`) is available at the deposit DOIs below for editors who prefer Word format.

### Files deposited externally and linked via Data Availability

The simulation code, raw outputs, and pre-specification dispatches are deposited on Zenodo with two DOIs (both publicly accessible without account or credentials):

- **Pre-specification anchor**: `10.5281/zenodo.20060972` (release `v0.0-prereg`)
- **Reproducibility snapshot**: `10.5281/zenodo.20272479` (release `v1.0-submission`)

Both DOIs resolve to download pages on Zenodo where the full repository can be retrieved as a `.zip`. The same content is also browsable on GitHub at `https://github.com/khatvangi/templating-substrates` (no login required).

### Recompiling the manuscript from source

`manuscript-v16.tex` requires:

- **pdflatex** + **bibtex** (any modern TeX distribution: TeX Live ≥ 2020, MiKTeX, or MacTeX)
- **REVTeX 4.2** with the `apsrev4-2.bst` bibliography style (install from CTAN: `revtex.tds.zip`)
- The standard packages: `amsmath`, `amssymb`, `subcaption`, `hyperref`, `booktabs`, `xcolor`, `natbib`

Compile sequence: `pdflatex → bibtex → pdflatex → pdflatex` (four passes resolve cross-references and citations). All figure assets live in `paper/figures/v16/` (PDF + PNG dual-format for each of the 19 panels).

### Re-running the simulations

The 23 test scripts in `code/` are self-contained and reproducible from a minimal Python environment:

- Python 3.10
- NumPy 1.24
- Matplotlib 3.7

Each script can be run independently from the repository root (e.g., `python code/test_a1_mode1_scaling.py`). All random seeds are set explicitly (`np.random.seed(42)` plus a per-replicate seed `42 + cell_idx*100 + rep_idx`), so re-running reproduces the result CSVs in `results/` exactly. Long-running tests (B, C, D, E, and the H/H3/H4/H5 family) take minutes; the full apparatus and family-level sweep complete in under an hour on a single 64-core node.

### One file is intentionally NOT redistributed

Tests F, F2, F3, and the Drt3b family-level analysis depend on the 1,232-sequence Drt3b family alignment published as supplementary data by Sharma, Lee, Armijo, Wang, Gao et al. (*Science* 2026, DOI [10.1126/science.aed1656](https://doi.org/10.1126/science.aed1656)). These eight third-party files (`data_s1.xlsx` through `data_s8.sto`) are the publisher's intellectual property and are **not included** in our deposit.

Readers wishing to re-run the family-level analyses must download these files directly from the *Science* article's supplementary materials section and place them in a directory named `science.aed1656_data_s1_to_s8 deng/` at the repository root (note the trailing space and the word "deng" — this matches the layout the test scripts expect). Full instructions are in `data/README.md` within the deposit. Verification: the single most critical file is `science.aed1656_data_s3.fa`, which should contain exactly 1,232 sequences (`grep -c "^>" data_s3.fa` should return 1232).

All other tests in the deposit run without any third-party data and reproduce the corresponding result CSVs from seeded RNG alone.

### Licensing

- All simulation code, analysis scripts, and figure-generation pipeline: **MIT License** (see `LICENSE`)
- Simulation result CSVs, figures, and the manuscript draft: **CC-BY-4.0** (see `LICENSE-data.md`)
- The third-party Sharma/Deng data files retain their original *Science* / AAAS terms; please cite the original *Science* paper if used.
