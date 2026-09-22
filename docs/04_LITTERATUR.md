# Litteraturläsning (2026-09-21)

Punkt 2 i `03_RESULTAT.md` avsnitt 5. Bara öppna källor. Varje påstående har URL och, där det går, avsnitt, tabell eller figur.
Märkning: **[original]** = jag har läst fulltexten, **[abstract]** = bara abstract, **[sekundär]** = annan källa som återger originalet,
**[sammanfattad]** = läst via sammanfattande hämtning och därför mindre pålitlig, **[egen räkning]** = räknat här, inte hämtat ur en källa.
Formelkontrollen ligger i sessionens scratchpad och är inte sparad som script. Ska talen användas i något som visas måste de in i `repro/`.

---

## Vad detta ändrar i 03_RESULTAT.md avsnitt 4

| # | Står nu | Ska stå | Källa |
|---|---|---|---|
| 1 | "Gilis m.fl. gjorde inget av det" (om verklig felfrekvens) | Fel för bröstcancer. Version 2 av preprinten (19/8 2020) redovisar felfrekvensen: efter en omkastning av BRCA gav HPA:s procedur 320 av 17 040 gener med p < 0,001 (1,88 %, 19 gånger nominellt) och 6 410 med p < 0,05 (37,6 %). Vårt tillskott är alla 21 cancerformer, 400 000 dragningar i stället för en omkastning, och ett kalibrerat p per gen. | Gilis v2, huvudtext s. 2 och Suppl. 2.2.1 |
| 2 | "visade ... p ner till 2,6e-7" | Rätt (2,62E-07, Suppl. figur S3), men det är det svagare av deras två tal. Lyft 320 av 17 040 i stället. | samma |
| 3 | "fann inga säkerställda gener i bröstcancer efter korrektion" | Rätt men ofullständigt. Med Cox, ålder som confounder och BH 5 %: BRCA 0 gener, LIHC 2 812 (mot HPA:s 3 501). Med Efrons empiriska nollfördelning (locfdr, lfdr < 0,2): **0 gener i både BRCA och LIHC**. Det står mot våra 3 428 robusta i LIHC och måste bemötas, se punkt 1.4 nedan. | Gilis v2, s. 3 och Suppl. avsnitt 4 och 5 |
| 4 | "Preprinten ser inte ut att ha publicerats" | Bekräftat. bioRxiv:s API anger `published: NA`, Crossref har ingen "is-preprint-of"-relation, två versioner (18/3 och 19/8 2020), artikeltyp "contradictory results". Citerad 1 till 3 gånger. | se 1.5 |
| 5 | Altman och Lausen & Schumacher "inte lästa i original" | Altman 1994 och Lausen & Schumacher 1992 ligger bakom betalvägg och är fortfarande **inte lästa i original**. Talen (cirka 40 % falskt positiva vid nominella 5 %, P = 0,002 motsvarar verkliga 0,05) är bekräftade i Altmans eget brev i Br J Cancer 1998, som är läst i fulltext. Formeln är bekräftad i källkoden till R-paketet maxstat (Hothorn och Lausen). | se avsnitt 2 |
| 6 | "Det som är nytt här: ... den verkliga felfrekvensen" | **Felfrekvensen går att räkna ut med en formel från 1982.** Med HPA:s intervall (20:e till 80:e percentilen) ger Miller och Siegmunds formel p_korr = 0,0169 vid p_min = 0,001, alltså 17 gånger. Vår permutation gav 12 till 17 gånger. Permutationen tillför avvikelsen från asymptotiken i små kohorter, inte storleksordningen. Skriv det så. | se 2.3, egen räkning |
| 7 | (saknas) | TCGA-CDR, den överlevnadstabell HPA v2 själva anger att de använt, avråder från OS i testikelcancer och manar till försiktighet i KICH, PRAD, READ, THCA och BRCA. Vår poäng om för få dödsfall är alltså redan gjord av datakällan. | Liu m.fl. 2018, tabell 3, se 4.1 |
| 8 | "HPA version 2 (2024) använder samma procedur" | Rätt. Tillägg: artikeln nämner varken Gilis, Altman, multipel testning eller cutoff-problemet, och kallar p < 0,001 "stringent" och "rigorous". Korrekt referens är Yuan m.fl., "The Human Pathology Atlas for deciphering the prognostic features of human cancers", eBioMedicine 2025;111:105495, online 10/12 2024. | se avsnitt 3 |

Följdrättelse i **avsnitt 1** i 03_RESULTAT.md: "vilken tabell HPA använde är OKÄNT" stämmer inte längre. Artikeln anger TCGA-CDR (Liu m.fl. 2018)
och endpoint overall survival. Vi använder alltså samma tabell, och avvikelserna i ändtarm och tjocktarm behöver en annan förklaring.
Två kandidater ur källorna: (a) genfiltret skiljer sig mellan artikel (medeluttryck > 1 TPM) och webbplats (median under 1 pTPM märks som
icke-prognostisk), (b) proteinatlas.org säger fortfarande att kliniska data kommer från GDC:s första release 6/6 2016, vilket motsäger artikeln.

---

## 1. Gilis, Taelman, Davey, Martens, Clement (bioRxiv 2020) [original, v2, hela PDF:en inklusive supplement]

Källa: https://www.biorxiv.org/content/10.1101/2020.03.16.994038v2.full.pdf (27 sidor, version 2 postad 19/8 2020, CC BY 4.0).

### 1.1 Vad de granskade
- Uhlén m.fl., Science 2017, alltså HPA:s patologiatlas **version 1** (FPKM, 17 cancerformer). Supplementet kallar texten en
  "technical comment on" Science-artikeln (Suppl. s. 7).
- Två cancerformer: bröst (BRCA) och lever (LIHC), TCGA via GDC. "the procedures for other datasets ... are analogous" (Suppl. s. 7).
- HPA:s egna tal som de utgår från: BRCA 847 gener med p < 0,001 och 7 749 med p < 0,05 av 17 040. LIHC 3 501 respektive 10 632 av 16 397
  (Suppl. 2.2).

### 1.2 Permutationerna
- **Nollanalysen av HPA:s procedur bygger på en enda omkastning** av överlevnadsdata i BRCA (Suppl. 2.2.1, script
  `/BRCA/originalBRCA_problems.Rmd`). Resultat: **320 gener med p < 0,001 och 6 410 med p < 0,05**, mot väntade 17 och 852. Lägsta p 2,62E-07
  (figur S3). Största p-värdet i hela nollkörningen var 0,65 (figur S4).
  [egen räkning] 320/17 040 = 1,88 %, 6 410/17 040 = 37,6 %.
- Sex upprepade omkastningar av överlevnad per dataset används senare, men för en annan fråga: om Cox-modellens teoretiska nollfördelning
  håller (figur S11 och S12). En omkastning av uttrycket per gen används i figur S9 och S10.
- De räknar ingen felfrekvens per cancerform, inget kalibrerat p-värde per gen och ingen validering mot oberoende kohorter.

### 1.3 Deras korrektion
De korrigerar inte HPA:s procedur, de byter ut den: Cox-regression på kontinuerligt uttryck (log-CPM) med splineterm, ålder som confounder
i BRCA, ålder och BMI i LIHC, likelihood ratio-test, därefter Benjamini-Hochberg 5 %. Resultat: **BRCA 0 gener, LIHC 2 812 gener**
(huvudtext s. 3). 14 573 gener testade i BRCA, 13 376 i LIHC (tabell S1 och S3).

### 1.4 Den del vi inte kände till, och som träffar vår metod
- Även deras Cox-modell har fel nollfördelning. Med Efrons empiriska noll (paketet locfdr, lfdr < 0,2) återstår **ingen gen i vare sig BRCA
  eller LIHC** (Suppl. avsnitt 4 och 5, tabell S1 till S4 där lfdr = 1 för alla topp 10). De påpekar själva att den empiriska nollan passar
  data dåligt (figur S8, S13).
- De hänvisar till Efron (Large-Scale Inference, kap. 6): nollfördelningen kan fallera av fyra skäl, (I) matematiska antaganden,
  (II) korrelation mellan gener, (III) korrelation mellan patienter, (IV) omätta confounders, och **"the issue of failed mathematical
  assumptions is the only type of failure that can be addressed via permutations"** (Suppl. s. 16).
- Konsekvens för oss [extrapolation]: vår kalibrering är en permutation. Den rättar cutoff-valet (skäl I) men inte confounding (ålder, stadium)
  och inte korrelationen mellan gener. En läsare i Gent kommer att säga att våra 12 504 "robusta" gener är en övre gräns, och att LIHC 3 428
  står mot deras 0. Det ska stå i vår text innan någon annan säger det.

### 1.5 Kod, publicering, svar
- Kod: https://github.com/statOmics/pitfallsOfHumanPathologyAtlas (R Markdown, mappar BRCA och LIHC, Dockerfile, ingen licensfil, senaste
  push 30/11 2020, 2 stjärnor, inga öppna ärenden). Kontrollerat via GitHubs API 21/9 2026.
- Tidskrift: **nej**. `https://api.biorxiv.org/details/biorxiv/10.1101/2020.03.16.994038` ger `published: NA` för båda versionerna.
  Crossref listar den som preprint utan relation till någon artikel.
- Citeringar: Semantic Scholar 1, OpenAlex 3. Den enda tidskriftsartikeln är Jeuken, Tobin och Käll, PLoS Comput Biol 2022
  (doi 10.1371/journal.pcbi.1010020), som citerar Gilis i förbigående om confounders, inte om cutoff-problemet [original, PMC8989354].
  Lukas Käll sitter på KTH och SciLifeLab, samma miljö som HPA.
- Svar från HPA-gruppen (Uhlén, Mardinoglu): **inget funnet**. Sökt på webben, i HPA v2-artikelns text och referenslista. bioRxiv:s
  kommentarsfält och Science eLetters gick inte att läsa (OKÄNT).

---

## 2. Altman m.fl. 1994 och Lausen & Schumacher 1992

### 2.1 Vad som är läst
- Altman, Lausen, Sauerbrei, Schumacher, JNCI 1994;86(11):829-35, doi 10.1093/jnci/86.11.829: **inte läst, betalvägg**. PubMed-posten
  (PMID 8182763) saknar abstract.
- Lausen och Schumacher, Biometrics 1992;48:73-85, doi 10.2307/2532740: **[abstract]** via https://repository.essex.ac.uk/2515/.
  Abstractet säger att den asymptotiska nollfördelningen för en maximalt vald rangstatistik, med valet begränsat till en inre del av
  variabelns spann, är fördelningen för supremum av absolutbeloppet av en standardiserad Gaussprocess, att det gäller även vid censurerade
  observationer, och att de jämför Monte Carlo med en approximation. Inga tal i abstractet.
- Altman, "Suboptimal analysis using 'optimal' cutpoints", Br J Cancer 1998;78(4):556-7 **[original]**,
  https://pmc.ncbi.nlm.nih.gov/articles/PMC2063091/. Altman sammanfattar där sin egen artikel från 1994.

### 2.2 Kvantitativt (ur Altman 1998, som återger Altman 1994)
- "because of multiple testing the false-positive rate is around 40% rather than the nominal 5%"
- "the P-value is far too small (P = 0.002 corresponds to a genuine P = 0.05)"
- Räknade exempel i brevet: rapporterade p = 0,0224 och 0,0028 blir korrigerade 0,29 och 0,06.
- Vilket percentilintervall 40 % avser står inte i brevet. [egen räkning] Formeln nedan med ε = 10 % ger 0,045 för p_min = 0,002 och
  återger 0,29 och 0,06 exakt, så talen avser med stor sannolikhet valet 10:e till 90:e percentilen.

### 2.3 Formeln
Källa för exakt form: källkoden till R-paketet **maxstat** (Hothorn och Lausen, CRAN 0.7-26), funktionen `pLausen92` i `R/maxstat.R`:

```r
pLausen92 <- function(b, minprop=0.1, maxprop=0.9) {
  if (b < 1) return(1)
  db <- dnorm(b)
  p <- 4*db/b + db*(b - 1/b)*log((maxprop*(1 - minprop))/((1-maxprop)*minprop))
  max(p,0)
}
```

Alltså, med z = (1 − p_min/2)-kvantilen i standardnormalfördelningen, φ dess täthet, och cutoffs tillåtna mellan kvantilerna ε1 och ε2:

**p_korr = φ(z) · (z − 1/z) · ln[ ε2(1 − ε1) / ((1 − ε2) ε1) ] + 4 φ(z) / z**

Paketets vinjett (https://cran.r-project.org/web/packages/maxstat/vignettes/maxstat.pdf, avsnitt 2) säger att Lausen och Schumacher 1992
visade att gränsfördelningen är supremum av en standardiserad brownsk brygga och att Miller och Siegmunds approximation (Biometrics 1982)
därför kan användas. Formeln är alltså Miller och Siegmunds, tillämpad på log-rank av Lausen och Schumacher.

Samma formel, plus Altmans förenklingar, står i manualen till R-paketet OptimalCutpoints, funktionen `control.cutpoints`, Details **[sekundär]**
(https://cran.r-project.org/web/packages/OptimalCutpoints/OptimalCutpoints.pdf):
- ε = 5 %: p_korr = −3,13 · p_min · (1 + 1,65 · ln p_min)
- ε = 10 %: p_korr = −1,63 · p_min · (1 + 2,35 · ln p_min)
- "work well for low pmin values (0.0001 < pmin < 0.1)". Manualens beskrivning av ε_high är tvetydig, maxstat-koden är entydig. Använd den.

**Ingen förenkling finns publicerad för HPA:s ε = 20 %. Använd den allmänna formeln**, där ln-termen blir ln(0,8·0,8 / (0,2·0,2)) = ln 16 = 2,773.

[egen räkning, scipy] Formeln med ε1 = 0,2 och ε2 = 0,8:

| p_min | p_korr | Uppblåsning |
|---:|---:|---:|
| 0,05 | 0,354 | 7 ggr |
| 0,01 | 0,110 | 11 ggr |
| **0,001** | **0,0169** | **17 ggr** |
| 1e-4 | 2,29e-3 | 23 ggr |
| 1e-5 | 2,90e-4 | 29 ggr |
| 1e-6 | 3,5e-5 | 35 ggr |

- För verkliga 0,001 per gen krävs p_min under cirka 3,9e-5. För verkliga 0,05 krävs p_min under cirka 3,7e-3.
- Jämförelse: vår permutation 1,25 till 1,73 % vid 0,001, formeln 1,69 %. Gilis en omkastning i BRCA: 1,88 % vid 0,001 och 37,6 % vid 0,05,
  formeln 1,69 % och 35,4 %. Tre oberoende vägar ger samma tal. Det stärker resultatet och minskar nyhetsvärdet.
- Förbehåll: formeln är asymptotisk. Vinjetten anger att Lausen m.fl. 1994 (förbättrad Bonferroni) och Hothorn och Lausen 2003 (exakt, små
  stickprov) finns för små n. Att våra permutationstal ligger under formeln i kohorter med få dödsfall är väntat [extrapolation].
- Formeln rättar bara cutoff-valet per gen. Multipel testning över 20 000 gener är ett separat steg, precis som hos oss.

---

## 3. HPA v2-artikeln [original, fulltext via Europe PMC, PMC11683280]

Yuan M, Zhang C, von Feilitzen K, Zwahlen M, Shi M, Li X, Yang H, Song X, Turkez H, Uhlén M, Mardinoglu A. "The Human Pathology Atlas for
deciphering the prognostic features of human cancers." eBioMedicine 2025;111:105495. doi 10.1016/j.ebiom.2024.105495.
https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(24)00531-0/fulltext

- **Nämner den problemet? Nej.** Orden Gilis, Altman, Lausen, multiple testing, false positive, false discovery, Bonferroni och adjusted
  förekommer inte i artikeltexten eller referenslistan (textsökning i hela XML-filen). Inte heller Uhlén 2017 Science gör det (fulltext via
  https://research.chalmers.se/publication/251875/file/251875_Fulltext.pdf).
- Diskussionen säger i stället: "our methodology for selecting PGs was stringent, utilizing a p-value threshold of less than 0.001 to ensure
  robust statistical significance. This rigorous cut-off minimizes the influence of potential gene expression fluctuations". Referenserna
  15 till 17 som följer är gruppens egna artiklar om läkemedelsrepositionering.
- Metod, avsnitt "Statistics": KM och log-rank på **overall survival**, alla TPM-värden från 20:e till 80:e percentilen prövas, lägsta p
  väljs, prognostisk om p < 0,001. R 4.2.3, paketet survival 3.5.5. Ogynnsam om högt uttryckt grupp har fler observerade händelser än väntat.
- Överlevnadstabell, avsnitt "Data collection": "This clinical information was sourced from the TCGA Pan-Cancer Clinical Data Resource
  (TCGA-CDR)" (ref. 18, Liu m.fl. 2018). 6 918 givare, 21 kohorter, rå-BAM via GDC, Kallisto, GRCh38, Ensembl 103.
- Filter: gener med medeluttryck > 1 inom cancerformen, patienter med överlevnadstid > 0 dagar.
- **Minsta gruppstorlek: anges inte.** Den enda begränsningen är percentilintervallet. I koden (`toolbox_pathology.R`, rad 233 till 263)
  är cutoffs de unika värdena strikt över 20:e och till och med 80:e percentilen, hög grupp är uttryck ≥ cutoff, p räknas med
  `survdiff(..., rho=0)` och `pchisq`.
- Kod: https://github.com/cellur-m/pathology_atlas (Apache 2.0, skapad 8/5 2024, senaste push 5/6 2024). Innehåller
  `code/survival_analysis.R`, `code/toolbox_pathology.R` och en exempelfil med kolumnerna `status` och `survival_time`. Inga patientlistor
  per kohort, så urvalet i urinblåsecancer går inte att läsa ur repot.
- Validering: en gen är "confidence prognostic gene" om den är prognostisk i både TCGA och uppföljningskohorten. Artikelns egen figur 4d:
  överlappet är signifikant (hypergeometriskt test, p < 0,05) i **4 av 10** cancerformer (GBM, KIRC, LIHC, LUAD). Spearman för KM-koefficienter
  0,64 i KIRC och 0,66 i LIHC. Det föregriper delvis vårt avsnitt 3: HPA visar själva att valideringen bär i ungefär samma cancerformer.
- Spår om urinblåsan [gissning]: TCGA-CDR tabell 3 har BLCA med 412 patienter och 181 dödsfall. HPA:s 169, alla avlidna, liknar ett urval
  på händelse snarare än hela kohorten. Artikeln förklarar det inte (den nämner bara att UCEC krympte 67,5 % för att BAM-filer saknades).

---

## 4. Andra kritiker och metodartiklar (sex stycken)

### 4.1 Liu m.fl., TCGA-CDR, Cell 2018;173:400-416 [original, PMC6066282]
https://pmc.ncbi.nlm.nih.gov/articles/PMC6066282/ . Detta är tabellen HPA v2 använder. Tabell 3 bedömer varje endpoint per cancerform:

| Cancerform | N | OS-händelser | OS-bedömning | Kommentar i tabell 3 |
|---|---:|---:|---|---|
| TGCT | 134 | 4 | avrådes (×) | "number of events is small for OS and DSS; need a longer follow-up" |
| PRAD | 500 | 10 | med försiktighet (✓*) | "need a longer follow-up for OS and DSS" |
| KICH | 113 | 13 | med försiktighet | "number of events is too small, need a longer follow-up" |
| THCA | 507 | 16 | med försiktighet | "number of events is small for OS and DSS" |
| READ | 170 | 26 | med försiktighet | "need a longer follow-up for OS, DSS, and DFI" |
| BRCA | 1 097 | 151 | med försiktighet | "need a longer follow-up for OS and DSS" |

Löptexten: "For an even less aggressive TCGA cancer type like PRAD, where there were only 10 OS events out 500 cases, OS is clearly not a
suitable study endpoint (Table 3)." Kravet för godkänd endpoint är bland annat fler än 20 händelser. HPA v2 citerar källan och använder
ändå OS i alla 21. **Vår observation om för få dödsfall är därmed inte ny, men kopplingen till HPA har vi inte sett någon göra.**

### 4.2 Smith och Sheltzer, Cell Reports 2022;38:110569 [original via PMC BioC, PMC9042322] och eLife 2018;7:e39217 [original, PMC6289580]
10 884 patienter, 33 cancerformer, 3 091 782 univariata Cox-modeller, 112 303 signifikanta par vid BH-FDR 1 %. Väljer Cox uttryckligen för
att slippa cutoffs: "unlike Kaplan-Meier analysis, Cox models do not require the selection of threshold values, so continuous data like gene
expression measurements do not need to be dichotomized" (STAR Methods). eLife 2018 citerar Uhlén 2017 som tidigare arbete utan att kritisera
det. Begränsning de själva anger: en kohort per cancerform, ingen oberoende validering. Relevans: en färdig, publicerad jämförelsetabell
för vår planerade känslighetsanalys med Cox (punkt 4 i avsnitt 5), och den naturliga invändningen "varför inte bara använda deras Z-värden?".

### 4.3 Lánczky och Győrffy, KMplot, J Med Internet Res 2021;23(7):e27633 [original, PMC8367126]
KM plotter har samma "best cutoff" (mellan nedre och övre kvartilen, Cox per cutoff, mest signifikanta väljs) men **räknar BH-FDR över
cutoff-värdena som standard**: "During the computation of multiple cut-off values, multiple hypotheses are generated. Therefore, the false
discovery rate (FDR) is computed by default". [egen bedömning] HPA är alltså strängt taget sämre än det mest använda verktyget i genren, som
åtminstone visar en korrektion. Någon publicerad kritik riktad specifikt mot KM plotters cutoff-val hittade jag inte.

### 4.4 Budczies m.fl., Cutoff Finder, PLoS ONE 2012;7(12):e51862 [original, PMC3522617]
Verktygsartikel som själv säger: "Cutoff optimization was demonstrated to contribute to overestimation of results when multiple cutoff points
are investigated and the problem of multiple hypothesis testing is ignored", och hänvisar till Altman 1994, Lausen och Schumacher 1992 samt
Miller och Siegmund 1982. Visar att problemet var allmänt känt i biomarkörfältet fem år före HPA v1.

### 4.5 Hothorn och Lausen, maxstat (R-paket och vinjett; metodartikel Comput Stat Data Anal 2003) [vinjett och källkod lästa, artikeln inte]
https://cran.r-project.org/package=maxstat . Färdig implementation av maximalt vald log-rank med korrigerat p (`smethod="LogRank"`,
`pmethod` Lau92, Lau94, exactGauss, HL, condMC). HPA:s R-pipeline kunde ha bytt `survdiff`-slingan mot ett anrop. Det är den minsta
möjliga rättelsen att föreslå dem.

### 4.6 Altman, Br J Cancer 1998 [original], se 2.1
Tre exempel ur samma tidskrift där korrigerade p-värden blev 0,29 och 0,06. Användbar som mall för tonen: saklig, kort, med räkneexempel.

Inte funnet: någon publicerad kritik av HPA:s prognostiska gener utöver Gilis m.fl. "Evaluate Cutpoints" (Ogłuszka m.fl., Comput Methods
Programs Biomed 2019) ligger bakom betalvägg och är inte läst.

---

## 5. Proteinatlas.org i dag (version 25.1, läst 21/9 2026)

- Metodsida: https://www.proteinatlas.org/humanproteome/cancer/method . Versionsnumret 25.1 står i gensidornas metadata.
- Metodtexten ("How has the data been analyzed?"): cutoffs från 20:e till 80:e percentilen, lägsta log-rank P väljs, "genes with log rank P
  values less than 0.001 in maximally separated Kaplan-Meier analysis were defined as prognostic genes". Gener med median under 1 pTPM
  klassas som icke-prognostiska oavsett p.
- **Ingen varningstext** om cutoff-optimering, multipel testning eller falska fynd. Samma sida anger Bonferroni-justering för
  proteomikjämförelserna (CPTAC), så frånvaron gäller just överlevnadsanalysen.
- Etiketter för användaren: sökfiltret "Prognosis" har fyra värden, "Favorable/Unfavorable - potential prognostic" och
  "Favorable/Unfavorable - validated prognostic". På en gensida, exempel https://www.proteinatlas.org/ENSG00000160211-G6PD/cancer :
  rubriken "PROGNOSTIC SUMMARY", texten "Kaplan-Meier plots for all cancers where high expression of this gene has significant (p<0.001)
  association with patient survival are shown", därefter "G6PD is a prognostic marker in:" med "Validated prognostic - unfavorable" för
  lever och "Potential prognostic - unfavorable" för KIRC och KIRP. Ordet "potential" är den enda reservationen.
- Definitionen av "validated" står inte på metodsidan. Artikeln definierar den som prognostisk i både TCGA och uppföljningskohort (se 3).
- Metodsidan säger att TCGA-data kommer från "the initial release of Genomic Data Commons (GDC) on June 6, 2016", vilket motsäger artikelns
  uppgift om TCGA-CDR. Vilket som gäller för 25.1 är OKÄNT.
- Webbplatsen låter användaren välja egen cutoff och subgrupp (stadium med mera) och ritar om KM-kurvan direkt, vilket öppnar för ännu en
  nivå av sökande efter lågt p [egen bedömning].

---

## Källor som inte gick att läsa
Altman m.fl. JNCI 1994 (betalvägg). Lausen och Schumacher Biometrics 1992 (bara abstract). Miller och Siegmund Biometrics 1982 (inte sökt i
original). Hothorn och Lausen CSDA 2003. Ogłuszka m.fl. 2019. bioRxiv-kommentarer och Science eLetters.
