# Kandidater A: sökning på reproducerbarhet (agent, 2026-09-21)

Geografin bortsågs från. Artikelsidor hos Nature, Springer och PubMed gick inte att hämta, texterna är lästa via PMC, bioRxiv eller
universitetsportaler, README och Zenodo direkt. Tal märkta "cirka" är inte kontrollerade mot källtexten. Betyg 1 till 5 på kriterium
2 (öppna data), 3 (körbar på bärbar dator), 5 (kod), 6 (nästa steg).

| # | Studie | Huvudresultat att reproducera | Data | Kod | 2 / 3 / 5 / 6 |
|---|---|---|---|---|---|
| 1 | **DrEval.** Bernett … List (TUM), Nat Commun 2026, doi 10.1038/s41467-026-72903-w | Djupa modeller för läkemedelsrespons slår knappt en naiv modell (medel per läkemedel och cellinje). Inget slår trädensembler. Preprintens tal, CTRPv2: slumpvis undanhållna par R² cirka 0,83 naiv mot 0,84 random forest. Osedda cellinjer, normaliserat R²: naiv cirka 0,08, random forest cirka 0,19, bästa djupa modell cirka 0,11. Osedda läkemedel: cirka 0. | Zenodo 10.5281/zenodo.12633909, v47 (2026-08-05), 5,5 GB. CTRPv2 519 MB, GDSC2 522 MB, författarnas egna resultat 2,7 GB (facit) | github.com/daisybio/drevalpy, Python, `pip install drevalpy`, poetry.lock, Dockerfile, MIT, push 2026-09-19 | 5 / 4 / 5 / 5 |
| 2 | **Chen … Ideker, Mattson** (UCSD), Bioinformatics 2026, doi 10.1093/bioinformatics/btag382 | AUC på trunkerat dosintervall (0,03–10 µM) ger bäst överensstämmelse mellan databaser, Spearman i medel 0,63 | CTRP, GDSC1/2, PRISM. Zenodo 20620409 (562 MB) | github.com/digitaltumors/Pharmacogenomics_Screens_Harmonization, Python, environment.yml, ingen licens | 5 / 5 / 4 / 5 |
| 3 | **Crawford … Greene**, Patterns 2024, doi 10.1016/j.patter.2024.101115 | TCGA till CCLE: bästa modellen vinner för 37 av 71 gener (52,1 %), minsta tillräckliga för 24 av 71 (figur 3B) | TCGA Pan-Cancer Atlas, CCLE 22Q2 | github.com/greenelab/pancancer-evaluation, Python, BSD, flera artiklar i samma repo | 5 / 4 / 3 / 4 |
| 4 | **Kebabci … Ryan** (UCD), Genome Medicine 2026, doi 10.1186/s13073-026-01759-y | Random forest för syntetisk letalitet mellan paraloger: ROC AUC 0,92, extern validering 0,87 och 0,81 (tal ur preprinten) | DepMap 22Q4, Ito-screenen, BioGRID, STRING | github.com/cancergenetics/context_specific_paralog_SL, notebooks, MIT | 4 / 4 / 4 / 5 |
| 5 | **Gross … Romagnoni** (Owkin), Sci Rep 2024, doi 10.1038/s41598-024-67023-8 | PCA och identitet matchar djupa representationer för överlevnad på TCGA (C-index skiljer högst 1,5 %) | zip på Owkins S3 | github.com/owkin/drl-evaluation, Python 3.9, poetry | 4 / 4 / 4 / 3 |
| 6 | Branson … Bessant (QMUL), Bioinformatics 2025 | Ingen prestanda kommer ur läkemedelsfeatures. BinaryET AUC 0,771 | GDSC2, manuell hämtning | github.com/Nik-BB/Understanding_DRP_models, ingen miljöfil | 4 / 3 / 2 / 3 |
| 7 | Codicè m.fl., J Cheminformatics 2025 | | | kräver CUDA | 4 / 2 / 2 / 3 |

**Agentens val: DrEval.** Ensam om pip-paket, låsta beroenden, versionerad data och publicerat facit. Huvudpåståendet är statistiskt,
inte biologiskt. Förslag på ett steg till: ersätt den naiva modellen med en mixed model med korsade slumpeffekter (läkemedel × cellinje
plus vävnad) och bootstrappa random forest-vinsten per läkemedel. Reserv: Chen, som bara kräver scipy.

Risker: Python 3.13 mot äldre låsta miljöer (kandidat 3, 5, 6). DepMap slutade enligt en enda sekundär källa publicera på Figshare efter 24Q4.
