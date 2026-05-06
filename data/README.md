# Third-party data dependencies

Several tests in this repository (Test E, Test E v2, the Drt3b family
conservation analysis, and Tests F, F2, F3) depend on supplementary
materials from:

> Sharma, P., Lee, H., Armijo, C., Wang, H., Gao, A., et al. (2026).
> Protein-templated synthesis of dinucleotide repeat DNA by an
> antiphage reverse transcriptase. *Science* (First Release, 16 Apr 2026).
> DOI: [10.1126/science.aed1656](https://doi.org/10.1126/science.aed1656)

We do **not** redistribute those files in this repository. They are the
publisher's intellectual property; download them directly from the
*Science* article's "Supplementary materials" section.

## Files needed

The Drt3b conservation analysis (`code/drt3b_conservation_analysis.py`
and Test F's family sweep) expect these files in
`science.aed1656_data_s1_to_s8 deng/` at the repository root:

| Filename | Format | Purpose |
|----------|--------|---------|
| `science.aed1656_data_s1.xlsx` | Excel | Strain metadata |
| `science.aed1656_data_s2.fa` | FASTA | Drt3a sequences |
| `science.aed1656_data_s3.fa` | FASTA | **Drt3b 1,232-sequence family alignment** (used by tests F, F2, F3) |
| `science.aed1656_data_s4.treefile` | Newick | Drt3a phylogeny |
| `science.aed1656_data_s5.treefile` | Newick | Drt3b phylogeny |
| `science.aed1656_data_s6.treefile` | Newick | Combined phylogeny |
| `science.aed1656_data_s7.xlsx` | Excel | Biochemistry data |
| `science.aed1656_data_s8.sto` | Stockholm | Multiple alignment |

The single most important file is `science.aed1656_data_s3.fa` (the
1,232-sequence Drt3b family alignment). Tests F, F2, and F3 cannot
run without it.

## Download instructions

1. Visit the article landing page: <https://doi.org/10.1126/science.aed1656>
2. Find the "Supplementary materials" section (usually a sidebar link
   labelled "PDF + Data" or "Download materials").
3. Download all eight Data S1–S8 files.
4. Place them in a directory at the repository root literally named
   `science.aed1656_data_s1_to_s8 deng/` (note the trailing space and
   the word "deng" — this matches the layout the test scripts expect).

## Verification after download

```bash
ls "science.aed1656_data_s1_to_s8 deng/"
# should list 8 files: data_s1.xlsx through data_s8.sto

# quick sanity check on the family alignment:
grep -c "^>" "science.aed1656_data_s1_to_s8 deng/science.aed1656_data_s3.fa"
# expected: 1232
```

If the count is not 1232, the file is incomplete or the wrong file —
re-download from the *Science* link.

## Citation

If you use these data downstream of this repository, please cite:

> Sharma, P., Lee, H., Armijo, C., Wang, H., Gao, A., et al. (2026).
> Protein-templated synthesis of dinucleotide repeat DNA by an
> antiphage reverse transcriptase. *Science*, DOI: 10.1126/science.aed1656.
