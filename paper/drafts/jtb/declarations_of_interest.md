# Declaration of Competing Interests

**For the manuscript** (insert as a single line before the References section):

> Declarations of interest: none.

**For the separate `.docx` upload** (Elsevier provides a template at https://www.elsevier.com/journals/journal-of-theoretical-biology/0022-5193/guide-for-authors but plain text is also accepted):

> The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

---

## Generation note

This file is the source for both the in-manuscript line and the separate `.docx`. To produce the `.docx`:

```bash
pandoc declarations_of_interest.md --from=markdown --to=docx --output=declarations_of_interest.docx
```

## Confirmation required before submission

- All listed authors agree they have no competing financial interests
- All listed authors agree they have no personal relationships influencing the work
- No undisclosed funding sources
