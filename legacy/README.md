# legacy/

Archived versions of the manuscript, supplementary appendix, and figure
directories. Preserved here for audit-trail and historical reference; the
active submission artifacts live at the repository root.

Per the project convention (`CLAUDE.md`): "move abandoned project material
to legacy rather than treating it as active."

## What is here

### Manuscripts (superseded)

| File | Phase | Notes |
|---|---|---|
| `manuscript-v8.{tex,pdf}` | Phase 6 (PRX Life submission, first compiled) | Initial REVTeX 4.2 conversion from the elife draft. First version with Sreedhar + Boggavarapu co-authorship. 20 pages. |
| `manuscript-v11.{tex,pdf}` | Phase 7 (pre-submission convergence) | First numerical-audit pass (D6/D7/D2 HIGH/MEDIUM fixes applied). 24 pages. |
| `manuscript-v14.{tex,pdf}` | Phase 7 (tone calibration) | Chemical-physics voice pass; capacity/generation theorem separation. 23 pages. |
| `manuscript-v15.{tex,pdf}` | Phase 7 (Discussion scoping) | Added penultimate Discussion paragraph on orthogonal preconditions (metabolism, compartmentalization, sustained driving). 23 pages. |

### Supplementary appendix (superseded)

| File | Notes |
|---|---|
| `supplementary-v8.tex` | First version of the Mathematical Appendix, paired with `manuscript-v8.tex`. Used the original "inheritance theorem proof sketch" formulation of §5, later rewritten in v3 as the separate capacity / generation theorems. |

### Figure directories (superseded)

| Directory | Type | Notes |
|---|---|---|
| `paper/figures/v11/` | Symlink → v7 | Originally a symlink for v11; flattened to a directory of its own here. |
| `paper/figures/v14/` | Flat copy of v7 | Made flat for submission-tarball portability. |
| `paper/figures/v15/` | Flat copy of v7 | Same pattern as v14. |

The canonical figure assets live at `../paper/figures/v7/` (the original
Phase 5 figure pipeline output); the active manuscripts (v16, v18) both
read from `../paper/figures/v16/`, a flat copy of v7.

## What is NOT here (the active artifacts)

At the repository root:

- `manuscript-v18.{tex,pdf,docx}` — current submission manuscript (sole author Boggavarapu Kiran; bioRxiv-ready)
- `manuscript-v16.{tex,pdf,docx}` — immediately prior version (also sole author; the v1.0-submission tagged release on GitHub / Zenodo)
- `supplementary-v3.{tex,pdf,docx}` — current Mathematical Appendix (latest formal proofs; capacity/generation theorem split)
- `supplementary-v14.tex` — immediately prior supplementary (kept as second-latest reference)
- `paper/figures/v7/` (canonical) + `paper/figures/v16/` (active)

## Reproducibility caveat for the archived `.tex` files

The archived manuscripts reference figures via `\graphicspath{{paper/figures/vN/}}`
which assumes the working directory is the repository root. Compiling
`legacy/manuscript-v14.tex` directly from this `legacy/` directory will
not resolve those figure paths. If you need to recompile an archived
version, copy the `.tex` back to the repository root and either:

1. Generate the corresponding flat figure directory at `paper/figures/vN/`
   from `paper/figures/v7/` (`cp -r paper/figures/v7 paper/figures/v14`), or
2. Update the `\graphicspath` line in the copied `.tex` to point at
   `paper/figures/v7/`.

Each archived `.pdf` is the previously-rendered version corresponding
to its `.tex` source at the time of that submission round.
