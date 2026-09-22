# Reproduction note · draft, 22 September 2026 · not peer reviewed
  How many Human Pathology Atlas v2 prognostic genes survive correction for the cut-off search?
  A reproduction of the prognostic gene analysis in Yuan et al. (eBioMedicine 2025;111:105495, online December 2024) in all 21 TCGA cohorts, and what the label “prognostic” means once the cut-off search is accounted for.
  Written by Kim Dahlroth, with an AI coding agent. No affiliation in cancer research. Everything below is computed from three public files with code that is public at the repository linked at the end, and nothing is a biological claim. The statistics are not new; the coverage is. We are asking one question: is this useful to you, and if not, why not?

  98 %of HPA’s 22,946 prognostic gene–cohort pairs recovered by our reimplementation
  12–17×the real false-positive rate of HPA’s p < 0.001, under permutation, in every cohort
  10 of 21cohorts where the number of “prognostic” genes is within what chance alone produces
  1.8–2.7 %of labelled genes validate in LUSC, BRCA, COAD and READ, against 1.0–1.7 % for a random gene

What this is
The Human Pathology Atlas v2 labels a gene prognostic in a cancer when a log-rank test at the best of all expression cut-offs between the 20th and 80th percentile gives p < 0.001. The labels are shown per gene on proteinatlas.org. We reimplemented the procedure in Python from the Methods text, ran it on HPA’s public expression file (v25.1) and the TCGA-CDR survival table, and then asked what the threshold means once the cut-off search is accounted for.
None of the statistics is new. That a p-value minimised over cut-offs is not a p-value has been known since Miller and Siegmund (1982), Lausen and Schumacher (1992) and Altman et al. (1994). Gilis et al. (bioRxiv 2020) made this criticism of the Pathology Atlas v1 in breast and liver cancer, including a permutation showing 1.88 % of genes below 0.001 under the null. We could not find a journal version of that preprint, and v2 uses the same procedure without mentioning the issue. What we add is coverage: all 21 cohorts of v2, per-cohort numbers, per-gene corrected values, and three comparisons. Of these, only one is a result we have not seen elsewhere: the whole-cohort permutation test in section 3, which shows that in ten of the twenty-one cohorts the number of labelled genes is within what chance alone produces. The rest is bookkeeping that, as far as we can tell, nobody had done for v2. For honesty's sake: we computed first and found Gilis et al. and the maxstat formula afterwards; both confirmed what we had, and had we read them first, section 2 would have been one function call.

1. The reproduction holds

MeasureOursHPA
Patients6,9126,918
Gene–cohort pairs with p < 0.001, among the 298,015 pairs both tested23,48922,946
Same gene in same cohort22,47798 % of HPA’s
Spearman correlation between our p and HPA’s p0.97 to 0.9996 in 20 cohorts; 0.90 in rectal cancer (16 deaths)
Direction agrees (favourable / unfavourable)93 to 99.8 %

The six missing patients are samples without an overall-survival record in TCGA-CDR (BRCA 2, OV, SKCM, STAD, UCEC one each), so patient selection is: every HPA sample whose patient has a survival record. Of the pairs we call at p < 0.001, 1,012 are not labelled by HPA, and 469 of HPA’s labels we do not reach; an expression filter on median pTPM explains fewer than a tenth of the former, the rest is the survival table. Sections 2 to 4 count 23,506 pairs below 0.001 in our run: the 23,489 above plus 17 genes we test that HPA’s table does not list. These are our raw calls, not the labels shown on proteinatlas.org, which is why the totals differ. Figure 3 and the validation counts use only genes present in both the TCGA and the validation cohort, so its n is smaller (KIRC 7,484 of 7,563; LIHC 3,061 of 3,449).
Using the exact cut-off rule from the authors’ R code (generateKMplot) instead of our reading of the Methods changes nothing material. The residual difference in rectum and colon is most likely patient selection, which we could not determine.

2. At p < 0.001 the procedure’s real false-positive rate is 1.25 to 1.73 %
We permuted survival against expression, 400,000 null draws per cohort, and ran the identical procedure.
Figure 1Actual false-positive rate under permutation against the nominal threshold, one line per cohort. Blue HPA’s best cut-off log-rank. Orange univariate Cox on log2(pTPM + 1), likelihood-ratio test. The orange lines that leave the diagonal are TGCT (4 deaths) and THCA (16); the blue line that drops below 1e-4 is KICH (9 deaths), where the permutation distribution has too few distinct values at that tail.

The real false-positive rate at the nominal 0.001 is 12 to 17 times too high in every cohort.
A closed formula gives the same answer. Lausen and Schumacher’s correction, as implemented in R’s maxstat, says that with a 20th-to-80th-percentile search pmin = 0.001 corresponds to 0.0169. Our permutation estimates are 0.74 to 1.02 times that value. By the same formula, a real 0.001 needs pmin below about 4e-5. No permutation is needed to fix this; it is one function call.
With Benjamini–Hochberg at 5 % on corrected p-values, 12,409 (permutation) or 12,356 (formula) of the 23,506 labelled genes remain, in seven cohorts: KIRC 7,563, LIHC 3,428, KIRP 1,058, LUAD 326, HNSC 25, CESC 7, BRCA 2. A further 95 pass in KICH, which we exclude because it has nine deaths (next point). In thirteen of twenty-one cohorts no individual gene passes this correction chain; that is a statement about the labels, not evidence that no prognostic association exists there. Gilis et al. argue that permutation does not address gene–gene correlation or unmeasured confounding, so these counts are upper bounds.
Several cohorts cannot support any test: TGCT has 4 deaths, PRAD 9, KICH 9, THCA 16, READ 16. TCGA-CDR itself (Liu et al., Cell 2018, Table 3) advises against overall survival in TGCT and urges caution in PRAD, KICH, THCA, READ and BRCA. KICH carries 1,284 prognostic labels on nine deaths, and its 95 corrected survivors rest on the same nine events.

3. How many “prognostic” genes does chance alone produce in one dataset?
Section 2 gives the mean. Genes are correlated, so the count in a single dataset scatters far more than a binomial would. We permuted survival 200 times per cohort and ran the procedure on all genes each time.
Figure 2Genes below p < 0.001 in a permuted dataset (bar: median, line: 5th to 95th percentile of 200 permutations) as a share of the genes HPA labels prognostic in that cohort. Above 100 % the labelled count is within chance. TGCT (22 labelled genes) runs off the axis and is clipped at the top.

In a typical permuted dataset 100 to 190 genes come out “prognostic” (median). In one dataset in twenty it is 450 to 870. The worst of 200 draws gave 3,900 in THCA, 3,163 in KIRC and 2,676 in PAAD.
In ten of twenty-one cohorts the observed number of prognostic genes is within what permutation produces (one-sided permutation p > 0.05 for the count; with 200 permutations the resolution is 0.005): COAD, GBM, LUSC, OV, PRAD, SKCM, STAD, TGCT, THCA, UCEC. READ (10 of 200 permutations at or above the observed count) and BLCA (11 of 200) are borderline: with 200 permutations the 95 % Monte Carlo interval around p = 0.05 is roughly 0.02 to 0.08, so these two cannot be placed. The other nine carry more signal than chance (at most 5 of 200), BRCA the weakest.
The same correlation affects any per-gene multiple-testing correction. Under permutation, Cox with Benjamini–Hochberg at 5 % returned at least one “discovery” in 1.5 to 13.5 % of permuted datasets in cohorts with enough deaths, and when it did, often hundreds or thousands. Per-cohort gene counts in this note are therefore upper bounds; the whole-cohort permutation test above does not have this problem.

4. A test without a cut-off finds the same cancers, and more genes in pancreatic cancer
Univariate Cox regression on log2(pTPM + 1), likelihood-ratio test, Benjamini–Hochberg 5 %. Under the same permutations this test holds its level (0.06 to 0.17 % at nominal 0.1 % in the seventeen cohorts with at least 28 deaths).

CohortDeathsHPA labelBest cut-off, corrected, BH 5 %Cox, BH 5 %
KIRC · kidney clear cell1717,5637,5637,853
LIHC · liver1293,4493,4283,482
PAAD · pancreas921,48602,380
KIRP · kidney papillary441,8111,058926
LUAD · lung adenocarcinoma1801,507326750
HNSC · head and neck21289325170
CESC · cervix67927755
Other fourteen5,87097194
Total23,50612,50415,810

Two caveats on Cox as the reference. The likelihood-ratio test on log2(pTPM + 1) is sensitive to heavy tails, so a few extreme values can drive a p-value; and a gene with genuine threshold behaviour is missed by a log-linear model. Figure 3 is the empirical answer to the second: at equal list size the Cox ranking validates as well or better. Cox also assumes a log-linear effect and proportional hazards, and more Cox hits do not mean more biology. Neither test adjusts for stage or grade, so 7,853 significant genes in KIRC (about 40 % of the coding genome) is a statement about association with survival in an unadjusted cohort, not about prognostic value on top of stage.
In seven cohorts Cox finds no gene at all (LUSC, OV, READ, SKCM, STAD, COAD, THCA), where HPA lists 253 to 567 each. In pancreatic cancer the corrected best cut-off test cannot single out any gene while Cox finds 2,380: correcting the cut-off procedure keeps it valid but costs power compared with not dichotomising in the first place. Do not read the Cox counts in TGCT (22), KICH (154) or PRAD (2) as findings.
Validation at equal list size
Take as many genes as HPA labels prognostic, but rank by Cox p. The share that is prognostic in the same direction in HPA’s own validation cohort is equal or higher in eight of ten cohorts (GBM 29.4 vs 19.2 %, LUSC 3.2 vs 1.8 %), lower in PAAD (8.5 vs 9.8 %) and marginally lower in LIHC (34.2 vs 34.5 %). The outcome is defined with HPA’s own procedure, which favours HPA’s list, and is too generous in absolute terms for every list. In LUSC, BRCA, COAD and READ the labelled genes validate in 1.8 to 2.7 % of cases against 1.0 to 1.7 % for a random gene: the label carries almost no information there. HPA’s own analysis (hypergeometric test, Figure 4) finds a significant overlap between TCGA and validation prognostic genes in four of ten cancers, GBM, KIRC, LIHC and LUAD, which is the same observation.
Figure 3Share of listed genes that are prognostic (p < 0.001, same direction) in HPA’s validation cohort. HPA list against a Cox-ranked list of the same size.

What we think follows, and what we do not know

Follows

The per-gene label on proteinatlas.org separates “validated prognostic” from “potential prognostic”, so a user can see whether a gene replicated. What the page does not show is that “potential prognostic, p < 0.001” carries very different weight by cohort: unadjusted association is present in kidney (KIRC, KIRP), liver, lung adenocarcinoma and pancreatic cancer, weak in HNSC, CESC and BRCA, and indistinguishable from noise in ten cohorts, where the same label is shown. None of this is prognostic value on top of stage and grade, which neither HPA nor we adjust for.
The smallest fix keeps HPA’s procedure and replaces the p-value with maxstat’s corrected one, plus a multiple-testing adjustment. A larger fix ranks by a continuous-expression Cox model and keeps the cut-off only for drawing the Kaplan–Meier plot.

Unknown

Whether any of this matters for how the atlas is used. If the labels are only ever a starting point for wet-lab work, a high false-positive rate may be an accepted cost. That is the question we cannot answer ourselves.
Not done: no adjustment for age, stage, grade or sex (Gilis et al. adjust for age), which is the main threat to every gene count here; Breslow ties instead of R’s Efron; the authors’ R code was read, not run; Altman 1994 and Lausen & Schumacher 1992 were checked against the maxstat source, not read in the original.

Every number above is produced by a numbered script with a fixed seed from three public files: HPA v25.1 rna_cancer_sample.tsv.gz and cancer_prognostic_data.tsv, and TCGA-CDR (Liu et al., Cell 2018) via UCSC Xena, listed with checksums. Code, per-gene tables and this note: github.com/Pochtli137/pathology-atlas-recount. Contact: kim.dahlroth@gmail.com.
