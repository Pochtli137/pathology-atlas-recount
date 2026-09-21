# Kandidater B: sökning på nordiska forskargrupper (agent, 2026-09-21)

Artiklarna är lästa via ett hämtverktyg som sammanfattar sidan. Tal ska kontrolleras mot PDF före användning. Betyg 1 till 5 på
kriterium 2 (öppna data), 3 (bärbar dator), 5 (kod), 6 (nästa steg).

| # | Studie | Huvudresultat | Data | Kod | 2 / 3 / 5 / 6 |
|---|---|---|---|---|---|
| 1 | **Human Pathology Atlas v2.** Yuan … Uhlén, Mardinoglu (sist och korresponderande, KTH/SciLifeLab). eBioMedicine 2024, doi 10.1016/j.ebiom.2024.105495 | 6 918 patienter, 21 cancertyper, 19 652 gener. Prognostisk gen = log-rank p < 0,001 med den cutoff mellan 20:e och 80:e percentilen som ger lägst p. KIRC och LIHC flest. Överensstämmelse mot valideringskohorter r = 0,64 och 0,66 | proteinatlas.org: TCGA-uttryck per prov, prognostabell per gen (facit), valideringstabell | github.com/cellur-m/pathology_atlas, R, Apache 2.0 | 5 / 5 / 3 / 5 |
| 2 | **Spatial artifact detection (NRFE).** Ianevski … Aittokallio (FIMM, Oslo). iScience 2025, doi 10.1016/j.isci.2025.113470 | GDSC1 mot GDSC2 Spearman 0,59, 0,51 på dåliga plattor, 0,65 på de bästa | PRISM och GDSC på plattnivå | github.com/IanevskiAleksandr/plateQC, R | 5 / 4 / 3 / 5 |
| 3 | **PDAC, ofullständiga multi-omics.** Paja-García … López (Oslo). PLOS Comput Biol 2026-09-10, doi 10.1371/journal.pcbi.1014735 | 154 patienter, två kluster, 4,5 gånger högre dödsrisk vid sex månader | TCGA-PAAD, CPTAC-PDA | github.com/ocbe-uio/imoc_pdac, Python, R, Octave | 5 / 4 / 4 / 4 |
| 4 | Drug response profiles. Abdel-Rehim … King (Cambridge, Chalmers). Bioinformatics 2026 | Pearson 0,904 mot 0,834 på GDSC | GDSC, i repot | Python, scikit-learn | 5 / 5 / 5 / 4 |
| 5 | CTLA4 i basal-lik bröstcancer. Røssevold m.fl. (Oslo). Commun Med 2025 | 97,2 % utan fjärrrecidiv vid fem år vid högt CTLA4 | GEO GSE96058 | ingen kod | 4 / 5 / 1 / 4 |
| 6 | hEFS i pankreascancer. Zobolas … Aittokallio. BioData Mining 2026 | C-index 0,54 till 0,64 | TCGA-PAAD, CPTAC | R, 40 kärnor och 40 GB | 4 / 2 / 2 / 3 |
| 7 | Pan-cancer emQTL. Ankill … Fleischer (Oslo). PLOS Comput Biol 2024 | 772 730 emQTL i sex kluster | Xena PANCAN | R | 5 / 3 / 2 / 3 |

Agentens val: kandidat 1.
