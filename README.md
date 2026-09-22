# Pathology Atlas recount

Reproduction of the prognostic gene analysis in the Human Pathology Atlas v2 (Yuan et al., eBioMedicine 2025;111:105495) in all
21 TCGA cohorts, and what the label "prognostic" (best cut-off log-rank p < 0.001) means once the cut-off search is accounted for.

The note: `docs/05_NOTE_EN.md` (draft, not peer reviewed). Everything else is in Swedish: it is a working repository, not a paper.

## Reproduce

Python 3.13 with numpy, pandas, scipy, statsmodels, matplotlib. Three public inputs, listed with checksums in `data/KALLOR.md`:
HPA v25.1 `rna_cancer_sample.tsv.gz` (1.4 GB) and `cancer_prognostic_data.tsv.zip`, and the TCGA-CDR survival table via UCSC Xena.
Put them under `data/hpa/` and `data/xena/`, then run the numbered scripts in `repro/` in order. `01` builds per-cohort matrices
(about 10 minutes), `02` reproduces (2 minutes), `03` permutes (400,000 draws per cohort, about 25 minutes on 10 cores), `04` to `08`
are the comparisons. All seeds are fixed. Tests: `python3 -m pytest -q repro/`.

`repro/ut/` holds the per-gene outputs (our p, HPA's p, corrected p, q) so the numbers in the note can be checked without rerunning.

## Status

Written by Kim Dahlroth (no affiliation in cancer research) with an AI coding agent, September 2026, as a test of whether that
produces anything a researcher finds useful. Not yet seen by anyone in the field. Corrections welcome as issues.

Licence: MIT for code. HPA data under HPA's own terms, TCGA data open tier.
