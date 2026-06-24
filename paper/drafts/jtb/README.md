# JTB submission package

Drafts for the Editorial Manager submission of *"A substrate recursion
principle for biological information, with empirical anchoring through a
templating-mode taxonomy"* to the *Journal of Theoretical Biology*.

The source manuscript is `../jtb-mansucript.tex` (at the repository root),
which compiles to `../jtb-mansucript.pdf` (21 pp, 734 KB) via the standard
`pdflatex → bibtex → pdflatex → pdflatex` sequence with REVTeX 4.2.

SHA256 of the source: `073b99a87cc89af0361eef8e06e6794e76ffc1ad3a57338063c350375b810daf`

## Files in this directory

| File | Purpose | Format for EM |
|---|---|---|
| `presubmit_checklist.md` | End-to-end walk-through before clicking submit. **Read this first.** | (workflow doc) |
| `cover_letter_jtb.md` | Cover-letter prose | Paste into EM cover-letter field |
| `highlights.txt` | 5 bullets ≤85 chars each | Upload as separate file with "highlights" in name |
| `declarations_of_interest.md` | Competing-interests statement | Generate `.docx` via pandoc; upload separately |
| `credit_contributions.md` | CRediT role assignment | **DRAFT — assign real roles** before EM entry |
| `ai_disclosure_for_manuscript.tex` | AI disclosure section | Paste into `jtb-mansucript.tex` before References |
| `data_references_for_bib.bib` | Two `[dataset]` Zenodo refs | Append to `templating_substrates.bib` |
| `suggested_reviewers.md` | Five candidate reviewers with COI guidance | Enter individually in EM reviewer fields |
| `jtb_submission_spec.md` | Full submission spec (verbatim user-provided) | Reference doc |

## Key decisions still pending

1. **Authorship**: JTB manuscript has two authors (Sreedhar + Boggavarapu). PRX/bioRxiv path was sole-author. Resolve before submitting. (See `presubmit_checklist.md`.)
2. **CRediT roles for Sreedhar**: spec leaves these as `[TO BE SET]` — assign per actual contribution.
3. **Funding**: spec assumes unfunded — confirm.
4. **AI disclosure scope**: spec wording covers "language editing and clarity" — broaden if AI was used for analysis, code, or figure design.
5. **Reviewer slate**: five candidates suggested; vet each for current affiliation/email and COI.

## Manuscript-body items NOT yet in `jtb-mansucript.tex`

These need to be pasted into the `.tex` (with confirmation of authorship + scope) before the source-bundle upload:

- `\section*{Declaration of generative AI ...}` block (from `ai_disclosure_for_manuscript.tex`)
- `Declarations of interest: none.` one-liner before References (from `declarations_of_interest.md`'s short form)
- Funding statement before References
- Two `[dataset]` `@misc` entries appended to `templating_substrates.bib` (from `data_references_for_bib.bib`)

After those edits, recompile and re-export the PDF.

## Submission target

| Field | Value |
|---|---|
| Journal | Journal of Theoretical Biology (Elsevier) |
| Article type | Full research article |
| Review type | Single-anonymized (author names stay on manuscript) |
| Open access | Decline (subscription route, $0 to author) |
| Submission system | Editorial Manager |
