# JTB submission package — verbatim spec (reference document)

> **SUPERSEDED 2026-06-26 in part.** The spec below was the planning
> document drafted on 2026-06-24 with placeholders for several pending
> decisions. The resolved live state is now reflected in:
>
> - `README.md` (index + resolved-decision flags)
> - `presubmit_checklist.md` (driver document for EM submission)
> - the individual executable files (`cover_letter_jtb.md`,
>   `credit_contributions.md`, `data_references_for_bib.bib`)
>
> Resolved decisions superseding spec defaults:
>
> | Spec item | Spec default | Resolved (2026-06-26) |
> |---|---|---|
> | Authorship | Two authors (Sreedhar + Boggavarapu) | **Sole author: Boggavarapu Kiran** |
> | CRediT for Sreedhar | `[TO BE ASSIGNED]` | All roles attributed to the sole author |
> | Zenodo data refs | `author = {Kiran, B. and Volet, S.}` | `author = {Kiran, B.}` |
> | Cover-letter voice | "the authors" / "all authors approve" | Singular voice; "the author declares" |
> | Manuscript SHA `073b99...` | (pre-sole-author state) | Post-edit SHA differs; see `presubmit_checklist.md` for the current state |
>
> This document is preserved verbatim below for audit-trail purposes;
> do not consult it as the live submission driver.

---

This is the full submission spec as provided by the corresponding author on
2026-06-24, preserved here as the authoritative reference for the other
drafts in this directory. Each of the section-numbered items below has a
corresponding executable file in `paper/drafts/jtb/`.

---

Target: *Journal of Theoretical Biology* (Elsevier, Editorial Manager). Article type: **Full research article**. Review: single-anonymized (author names stay on the manuscript). Route: **subscription** (decline the OA/APC option at acceptance — no charge).

Manuscript source to compile for the review PDF: `manuscript-v18.tex`, sha256 `073b99a87cc89af0361eef8e06e6794e76ffc1ad3a57338063c350375b810daf`. This is the voice-edited, abstract-tightened, numbers-corrected version — not `prx-life-6-16.pdf`.

> Note (added at packaging time): the sha256 hash actually matches `jtb-mansucript.tex` (the JTB file placed at the repo root), not the current `manuscript-v18.tex`. Compile from the JTB file; treat `jtb-mansucript.tex` as the canonical source.

## 1. Highlights (mandatory)

Upload as a separate file with "highlights" in the filename; each ≤85 chars.

- Four conditions on a template-catalyst pair decide open-ended inheritance
- Six biological templating systems sorted by which condition each fails
- Drt3b makes alternating poly(AC) DNA with no nucleic acid template
- Predicted dG misincorporation matches Drt3b biochemistry within 0.16%
- Two Drt3b active-site mutants give single-experiment tests of the theory

Biology-forward by design — JTB requires highlights to "feature the biological applications as well as any theoretical advancements."

## 2. Keywords (≤5)

Templating; Heredity; Evolvability; Reverse transcriptase; Information theory

## 3. Cover letter

→ `cover_letter_jtb.md`

## 4. Declaration of competing interests

Manuscript line + separate `.docx`.

→ `declarations_of_interest.md` (text source for both)
→ `declarations_of_interest.docx` (generated via pandoc)

## 5. Funding statement (in manuscript)

"This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors."

(Assumes unfunded — confirm.)

## 6. CRediT author contributions

→ `credit_contributions.md` (DRAFT; assign per actual division of labour before submission)

## 7. Generative-AI use disclosure

→ `ai_disclosure_for_manuscript.tex` (LaTeX snippet to paste before References)

Default wording covers language editing only. Broaden scope to match actual AI use across the workflow.

## 8. Data statement and data references (JTB "Option C")

**Data statement (submission form):** "All code, simulation outputs, and pre-specified analysis plans supporting this study are openly available on Zenodo (https://doi.org/10.5281/zenodo.20272479) and GitHub (https://github.com/khatvangi/templating-substrates). The pre-specification anchor is archived at https://doi.org/10.5281/zenodo.20060972."

**Add to the reference list** (JTB tags datasets with a leading `[dataset]`, stripped at proof):

→ `data_references_for_bib.bib` (two `@misc` entries to append to `templating_substrates.bib`)

## 9. Suggested reviewers

→ `suggested_reviewers.md` (5 candidates; vet for current affiliation/email and COI)

## 10. File-upload plan

- **Manuscript (for review)**: author-compiled PDF of `jtb-mansucript.tex` (`jtb-mansucript.pdf`)
- **LaTeX source**: `.tex` + figure files + `.bib` (+ `.bbl`) zipped
- **Highlights**: separate file with "highlights" in name
- **Declaration of competing interest**: separate `.docx`
- **Supplementary mathematical appendix**: upload as supplementary material (`supplementary-v3.pdf`)
- **Figures (high-res)**: required at acceptance as separate EPS/PDF/TIFF ≥300 dpi; embedded figures fine for review

**Deferred to revision/acceptance, not blocking now**: port REVTeX → Elsevier `elsarticle`, switch numeric citations to author–year Harvard, apply numbered sections.

## Pre-submit checklist

→ `presubmit_checklist.md`
