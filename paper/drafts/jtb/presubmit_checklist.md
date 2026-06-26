# JTB Pre-Submit Checklist

Walk through this top to bottom before clicking submit in Editorial Manager.
Each item is either [verify], [confirm], or [paste/upload].

## Authorship and identity

- [x] **Authorship resolved 2026-06-26: sole author (Boggavarapu Kiran).** Applied to `jtb-mansucript.tex`, `manuscript-v18.tex`, `cover_letter_jtb.md`, `credit_contributions.md`, and `data_references_for_bib.bib`. The sole-author state is now consistent across the JTB, PRX/bioRxiv, and Zenodo tracks.
- [ ] Corresponding author email confirmed in EM (`kiran@mcneese.edu`)
- [ ] Corresponding author full postal address entered in EM
- [ ] Corresponding author phone number entered in EM (required field; +1 area code + 10 digits)
- [ ] ORCID linked for each author (CRediT entries map to ORCIDs)

## Manuscript-body items to paste in BEFORE recompiling

Edit `jtb-mansucript.tex` to add the following before `\bibliography{templating_substrates}`:

- [ ] Paste the `\section*{Declaration of generative AI...}` block from `ai_disclosure_for_manuscript.tex`
- [ ] Add the one-line: `\noindent\textbf{Declarations of interest:} none.`
- [ ] Add the funding statement: `\noindent\textbf{Funding:} This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.` (or substitute actual funding source)
- [ ] Add a CRediT section: `\section*{CRediT authorship contribution statement}` with the agreed role assignments from `credit_contributions.md`

Edit `templating_substrates.bib` to append:

- [ ] Append the two `@misc` entries from `data_references_for_bib.bib`

Then cite the Zenodo entries in the Data Availability paragraph:

- [ ] Add `\citep{Zenodo_submission_2026, Zenodo_prespec_2026}` to the Data Availability sentence

## Recompile after the edits

- [ ] `pdflatex jtb-mansucript` (pass 1)
- [ ] `bibtex jtb-mansucript`
- [ ] `pdflatex jtb-mansucript` (pass 2)
- [ ] `pdflatex jtb-mansucript` (pass 3)
- [ ] Visual check: AI disclosure renders before References; data-reference entries show with `[dataset]` prefix; CRediT section renders; declarations + funding render
- [ ] PDF page count: was 21 before additions; expect 22 or 23 after

## EM uploads

- [ ] **Manuscript (main file, for review)**: upload `jtb-mansucript.pdf` (recompiled per above)
- [ ] **Highlights**: upload `highlights.txt` renamed to `highlights_v1.txt` (EM wants "highlights" in the filename)
- [ ] **Declarations of competing interest**: generate `declarations_of_interest.docx` via `pandoc declarations_of_interest.md --to=docx --output=declarations_of_interest.docx`, upload
- [ ] **Source bundle (LaTeX)**: zip `jtb-mansucript.tex` + `templating_substrates.bib` + `paper/figures/v16/` + the `.bbl` from the recompile, upload as `source.zip`
- [ ] **Supplementary material**: upload `supplementary-v3.pdf` (unchanged from PRX path)

## EM form fields

- [ ] **Article type**: Full research article
- [ ] **Title**: A substrate recursion principle for biological information, with empirical anchoring through a templating-mode taxonomy
- [ ] **Abstract**: paste from `jtb-mansucript.tex` lines 95–97 (the rewritten JTB-style compact abstract, ≤280 words)
- [ ] **Keywords**: `Templating; Heredity; Evolvability; Reverse transcriptase; Information theory`
- [ ] **Cover letter**: paste from `cover_letter_jtb.md` (insert bioRxiv DOI when assigned)
- [ ] **Suggested reviewers**: enter 3–5 from `suggested_reviewers.md` (in the recommended order; vet emails first)
- [ ] **Suggested editor**: (optional) leave blank unless you know one
- [ ] **Data statement**: paste the data statement from `presubmit_checklist.md` § "Data statement form text"
- [ ] **Funding source**: declare per the manuscript-body funding statement
- [ ] **Open access**: **DECLINE** at acceptance (subscription route, $0 to author)
- [ ] **Review type**: confirm "single-anonymized" — author names stay on manuscript

## Data statement form text (paste into EM data-statement field)

> All code, simulation outputs, and pre-specified analysis plans supporting this study are openly available on Zenodo (https://doi.org/10.5281/zenodo.20272479) and GitHub (https://github.com/khatvangi/templating-substrates). The pre-specification anchor is archived at https://doi.org/10.5281/zenodo.20060972.

## Pre-submit grep tests

```bash
cd /storage/kiran-stuff/templating_framework
grep -c "TBD\|TODO\|FIXME" jtb-mansucript.tex   # should be 0
grep -c "Anonymous" templating_substrates.bib   # should be 0
grep -c "DOI to be inserted" paper/drafts/jtb/cover_letter_jtb.md  # update with bioRxiv DOI when assigned
```

## Final sanity check before clicking submit

- [ ] PDF opens cleanly in two different readers (Acrobat + Preview / browser)
- [ ] All 6 figures visible with correct panel labels
- [ ] All 6 tables visible and within margins
- [ ] References render with full author lists (no "Author A, et al." truncation unless that's the journal style)
- [ ] No `??` in citation positions (would indicate bib resolution failure)
- [ ] No `[?]` in figure cross-references (would indicate label resolution failure)
- [ ] Page count is reasonable (21–23 pp expected)

## After acceptance (do not do at submission)

- [ ] Decline OA/APC option (subscription route)
- [ ] Submit copyedited proofs within 48-hour window
- [ ] Deposit accepted manuscript to institutional repository if your funder mandates it
