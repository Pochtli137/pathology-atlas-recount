# Mail till Adil Mardinoglu (utkast, Kim skriver om och skickar själv)

Mottagare: adilm@scilifelab.se (korresponderande författare enligt artikeln, PMC11683280, primär). Science for Life Laboratory,
KTH, Stockholm, och King's College London. Kopia: ingen. Bilaga: länk till artefakten (delad först när Kim valt det) och, om han
vill, repot.

---

Subject: Reproduction of the Pathology Atlas v2 prognostic genes in all 21 cohorts, and a question

Dear Adil,

I reproduced the prognostic gene analysis in your eBioMedicine 2024 paper (Human Pathology Atlas v2) from HPA's public expression
file and the TCGA-CDR survival table, in all 21 TCGA cohorts. The reproduction holds: 98 % of your 22,946 prognostic gene–cohort
pairs come back, with Spearman 0.97 to 0.9996 against your published p-values.

Then I asked what p < 0.001 means once the best cut-off search is accounted for. Under permutation the procedure's real
false-positive rate is 1.25 to 1.73 % in every cohort, which Lausen and Schumacher's 1992 formula also gives (0.0169). After
correction and Benjamini–Hochberg at 5 %, 12,504 of the 23,506 labelled genes remain, in eight cohorts. In ten cohorts the number of
labelled genes is within what permutation alone produces. A Cox model on continuous expression finds the same cancers and more genes
in pancreatic cancer, where the corrected cut-off test finds none.

I know this is not new: Gilis et al. raised it against v1 in 2020, and the statistics go back to Altman 1994. What I add is the
coverage, per-cohort numbers, per-gene corrected values, and a validation comparison in your own validation cohorts. A two-page note
with three figures is here: [länk]. Every number comes from a numbered script with a fixed seed on public files.

I am not a cancer researcher. I work with experimentation and statistics in e-commerce, and this was a test of whether that skill
set produces anything a researcher would find useful. So my question is simple: is this useful to you, and if not, why not? If the
labels are only ever a starting point for wet-lab work, a high false-positive rate may be an accepted cost, and I would like to
understand that.

Happy to share the code, and happy to be told where I have gone wrong.

Best regards,
Kim Dahlroth
Stockholm
kim.dahlroth@gmail.com

---

Anteckningar:
- Kort, fakta först, frågan sist. Ingen anklagelse: "not new", "where I have gone wrong".
- Artefakten är privat tills Kim delar den. Skicka inte länken innan den är delad, annars får han 404.
- Svarar han inte inom två veckor: Lieven Clement (Gent) är den som redan bryr sig om frågan.
