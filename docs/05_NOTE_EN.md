Reproduction note · draft, 23 September 2026 · not peer reviewed
  # Human Pathology Atlas v2 prognostic genes: what survives the cut-off search, and what survives stage
  A reproduction of the prognostic gene analysis in Yuan et al. (eBioMedicine 2025;111:105495, online December 2024) in all 21 TCGA cohorts, with two corrections the atlas does not make: for the search over cut-offs, and for tumour stage.
  Written by Kim Dahlroth, with an AI coding agent. No affiliation in cancer research. Everything below is computed from three public files with code that is public at the repository linked at the end, and nothing is a biological claim. The question to the atlas team is at the end.

  98 %of HPA’s 22,946 prognostic gene–cohort pairs recovered by our reimplementation
  12–17×the real false-positive rate of HPA’s p < 0.001, under permutation, in every cohort
  811 → 4Cox-significant genes in lung adenocarcinoma once stage is held fixed; 805 → 0 in papillary kidney, 1,711 → 5 in pancreas
  7,806 → 4,203in clear cell kidney, and 2,377 → 1,117 in liver: signal beyond stage remains, about half of it

## What this is
The Human Pathology Atlas v2 labels a gene prognostic in a cancer when a log-rank test at the best of all expression cut-offs between the 20th and 80th percentile gives p < 0.001, in overall survival, without adjustment for stage, grade or age. The labels are shown per gene on proteinatlas.org. We reimplemented the procedure in Python from the Methods text, ran it on HPA’s public expression file (v25.1) and the TCGA-CDR survival table, and asked two questions: what the threshold means once the cut-off search is accounted for, and how much of the signal is tumour stage.
The first question is not new. That a p-value minimised over cut-offs is not a p-value has been known since Miller and Siegmund (1982), Lausen and Schumacher (1992) and Altman et al. (1994). Gilis et al. (bioRxiv 2020) made this criticism of the Pathology Atlas v1 in breast and liver cancer, adjusting for age. We could not find a journal version of that preprint, and v2 uses the same procedure without mentioning the issue. Section 2 is bookkeeping for v2. The stage analysis in section 4 is the part we have not seen elsewhere. For honesty’s sake: we computed first and found Gilis et al. and the maxstat formula afterwards; had we read them first, section 2 would have been one function call.

## 1. The reproduction holds

| Measure | Ours | HPA | 
| Patients | 6,912 | 6,918 | 
| Gene–cohort pairs with p < 0.001, among the 298,015 pairs both tested | 23,489 | 22,946 | 
| Same gene in same cohort | 22,477 | 98 % of HPA’s | 
| Spearman correlation between our p and HPA’s p | 0.97 to 0.9996 in 20 cohorts; 0.90 in rectal cancer (16 deaths) | 
| Direction agrees (favourable / unfavourable) | 93 to 99.8 % | 

The six missing patients are samples without an overall-survival record in TCGA-CDR, so patient selection is: every HPA sample whose patient has a survival record. Of the pairs we call at p < 0.001, 1,012 are not labelled by HPA, and 469 of HPA’s labels we do not reach; the survival table explains most of the difference. Sections 2 and 3 count 23,506 pairs below 0.001: our raw calls plus 17 genes HPA’s table does not list. Using the exact cut-off rule from the authors’ R code (generateKMplot) instead of our reading of the Methods changes nothing material.

## 2. At p < 0.001 the procedure’s real false-positive rate is 1.25 to 1.73 %
We permuted survival against expression, 400,000 gene-level null draws per cohort, and ran the identical procedure.
Figure 1Actual false-positive rate under permutation against the nominal threshold, one line per cohort. Best cut-off log-rank Univariate Cox. The orange lines that leave the diagonal are TGCT (4 deaths) and THCA (16); the blue line that drops below 1e-4 is KICH (9 deaths), where the permutation distribution has too few distinct values at that tail.
A closed formula gives the same answer. Lausen and Schumacher’s correction, as implemented in R’s maxstat, says that with a 20th-to-80th-percentile search pmin = 0.001 corresponds to 0.0169. Our permutation estimates are 0.74 to 1.02 times that. The formula assumes a continuous statistic; genes with many zeros in pTPM offer fewer distinct cut-offs, which is the likely reason permutation lands below the formula in some cohorts. A real 0.001 needs pmin below about 4e-5.
With Benjamini–Hochberg at 5 % on corrected p-values, 12,409 of the 23,506 labelled pairs remain, in seven cohorts (table in section 3). Five cohorts cannot support any test on overall survival: TGCT has 4 deaths, PRAD 9, KICH 9, THCA 16, READ 16; TCGA-CDR itself (Liu et al., Cell 2018, Table 3) advises against overall survival in TGCT and urges caution in the others. KICH carries 1,284 prognostic labels on nine deaths.

## 3. A calibrated test without a cut-off
Univariate Cox regression on log2(pTPM + 1), likelihood-ratio test, Benjamini–Hochberg 5 %. Under permutation it holds its level (0.06 to 0.17 % at nominal 0.1 %) in the cohorts with at least 28 deaths, though its extreme tail is inflated below about 90 deaths.

| Cohort | Deaths | HPA label | Best cut-off, corrected | Cox | Cox, staged patients | Cox, stratified on stage | 
| KIRC · kidney clear cell | 171 | 7,563 | 7,563 | 7,853 | 7,806 | 4,203 | 
| LIHC · liver | 129 | 3,449 | 3,428 | 3,482 | 2,377 | 1,117 | 
| PAAD · pancreas | 92 | 1,486 | 0 | 2,380 | 1,711 | 5 | 
| KIRP · kidney papillary | 44 | 1,811 | 1,058 | 926 | 805 | 0 | 
| LUAD · lung adenocarcinoma | 180 | 1,507 | 326 | 750 | 811 | 4 | 
| HNSC · head and neck | 212 | 893 | 25 | 170 | 77 | 108 | 
| CESC · cervix | 67 | 927 | 7 | 55 | 61 | 247 | 
| Other fourteen |  | 5,870 | 97 | 194 | see section 4 | 

In seven cohorts Cox finds no gene at all: LUSC (210 deaths), OV (216), STAD (140), COAD (61), SKCM (28), READ (16), THCA (16), where HPA lists 253 to 567 genes each. In the first three the test has the power to find something and does not. That is the plainest evidence that the label in those cohorts is mostly the search over cut-offs. A permutation test on the number of labelled genes per cohort points the same way but is weaker: among the sixteen cohorts with at least 28 deaths, in eight (BLCA, COAD, GBM, LUSC, OV, SKCM, STAD, UCEC) the count is within what permuted data produce. The permuted counts spread widely (Figure 2), so that is absence of evidence for signal, not evidence of its absence.
Figure 2Genes below p < 0.001 in a permuted dataset (bar: median, line: 5th to 95th percentile of 200 permutations, 2,000 in eight cohorts) as a share of the genes HPA labels prognostic in that cohort. Above 100 % the labelled count is within chance. TGCT (22 labelled genes) is clipped.
Two caveats on Cox. It assumes a log-linear effect and proportional hazards, and a gene with genuine threshold behaviour is missed. The likelihood-ratio test on log2(pTPM + 1) is sensitive to heavy tails: on within-gene ranks instead, the counts barely change in KIRC, LIHC and LUAD, but PAAD drops from 1,711 to 657, so the pancreatic count rests partly on extreme values.

## 4. In three cancers the label cannot be told apart from stage
Stage predicts survival, and thousands of genes track stage, so an unadjusted association can be a stage marker rather than a prognostic gene. We refitted Cox stratified on AJCC pathological stage I–IV from TCGA-CDR (FIGO clinical stage for CESC, OV and UCEC), a separate baseline hazard per stage and one coefficient per gene, on the patients with a recorded stage. Dropping unstaged patients already moves the counts (LIHC loses 14 of 129 deaths), so compare the last two columns of the table above.

- Lung adenocarcinoma, papillary kidney and pancreas lose almost everything: 811 → 4, 805 → 0, 1,711 → 5. On ranks: 811 → 3, 867 → 2, 657 → 1.
- Clear cell kidney and liver keep about half: 7,806 → 4,203 and 2,377 → 1,117. Stratifying on stage and grade together leaves 2,407 and 377, a model we have not calibrated.
- The stratified test holds its level: with expression permuted within stage, 0.05 to 0.11 % of genes fall below p = 0.001 in KIRC, LIHC, LUAD, KIRP, HNSC and CESC (20 permutations each). In PAAD it is 0.21 %; stages III and IV hold 3 and 4 patients there. KICH and TGCT cannot be stratified.
- Papillary kidney has 39 staged deaths and pancreas is 84 % stage II, so lost power may contribute there. Lung adenocarcinoma has 178 deaths, and that explanation does not hold.
- In HNSC, CESC and BRCA the count rises with stratification (77 → 108, 61 → 247, 7 → 34), with the test calibrated. We do not have an explanation.

## 5. The labels replicate, and so would a stage marker
HPA reports each gene also in an independent validation cohort for ten cancers. We count a listed gene as validated when the validation p, corrected for the cut-off search with the formula above, is below 0.001 in the same direction. (Using the uncorrected p makes the base rate largely the procedure’s own false-positive rate, which an earlier version of this note did.)

| Cohort | List size | Random gene | HPA list | Cox list, same size | 
| LIHC | 3,061 | 2.9 % | 15.2 % (464) | 14.7 % | 
| KIRC | 7,484 | 1.9 % | 6.3 % (470) | 6.3 % | 
| GBM | 364 | 1.3 % | 4.7 % (17) | 11.0 % | 
| LUAD | 1,462 | 0.4 % | 2.3 % (34) | 3.2 % | 
| PAAD | 1,404 | 0.5 % | 2.3 % (32) | 1.8 % | 
| BRCA, COAD, LUSC, OV, READ | 280–571 | 0.1–0.4 % | 0 to 3 genes each | 0 to 4 | 

Where the numbers are large enough to say anything, labelled genes replicate at three to five times the rate of a random gene. The labels carry information. But the validation cohorts are not stage-adjusted either, so a gene that marks stage replicates as well as one that does not: lung adenocarcinoma replicates at five times the base rate while stage-stratified Cox keeps four genes. Replication cannot separate a prognostic gene from a stage marker. At equal list size a Cox ranking replicates better in GBM (11.0 against 4.7 %) and LUAD, equally in KIRC, slightly worse in LIHC and worse in PAAD.

## What we think follows, and what we do not know

### Follows

- The per-gene label “potential prognostic, p < 0.001” carries very different weight by cohort. In clear cell kidney and liver there is association with survival beyond stage. In lung adenocarcinoma, papillary kidney and pancreas the label cannot, on this analysis, be told apart from a stage proxy. In LUSC, OV and STAD a calibrated test finds nothing. The page shows the same label in all of them.
- The smallest fix keeps HPA’s procedure and replaces the p-value with maxstat’s corrected one. The next one stratifies on stage, which is in the same TCGA table. Both are a few lines of R.

### Unknown

- How the labels are used. If they are only a starting point for wet-lab work, false positives and stage markers may be an accepted cost. We cannot judge that.
- Why stratifying on stage increases the count in HNSC, CESC and BRCA.
- Not done: adjustment for age, sex or treatment; calibration of stage-and-grade stratification; stage as a covariate instead of strata; Breslow ties instead of R’s Efron; the authors’ R code was read, not run; the whole-cohort permutation test in the strongest cohorts has 200 permutations (resolution 0.005); Altman 1994 and Lausen & Schumacher 1992 were checked against the maxstat source, not read in the original.

## The question
Why does proteinatlas.org show an uncorrected, unadjusted p-value and label, when the corrected value is one function call and stage is in the same table? In three of the five cancers with the most labels, the label appears to mark stage. If there is a reason to show it anyway, or a flaw in what we did, we would like to know which.

Every number above is produced by a numbered script with a fixed seed from three public files: HPA v25.1 rna_cancer_sample.tsv.gz and cancer_prognostic_data.tsv, and TCGA-CDR (Liu et al., Cell 2018) via UCSC Xena, listed with checksums. Stage analysis: repro/11_stadium.py, calibration 12, validation 13. Code, per-gene tables and this note: github.com/Pochtli137/pathology-atlas-recount. Contact: kim.dahlroth@gmail.com.
