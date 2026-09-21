# Vald studie: Human Pathology Atlas v2

Beslut 2026-09-21 (Claude, enligt kriterierna i `01_URVAL.md`; Kim kan riva upp det).

**Studien:** Yuan m.fl., "The Human Pathology Atlas for deciphering the prognostic features of human cancers", eBioMedicine 2024,
doi 10.1016/j.ebiom.2024.105495. Sista och korresponderande författare Adil Mardinoglu, KTH och SciLifeLab, Stockholm. Mathias Uhlén
medförfattare. Resultaten visas per gen på proteinatlas.org, där forskare över hela världen slår upp om en gen är "prognostisk".

## Varför den, och inte DrEval

- **Forskaren finns i Stockholm.** Steg 4 är att visa resultatet för en forskare, och en författare har mest skäl att titta. DrEval
  (TU München) var bäst paketerad och är reserv.
- **Facit finns per gen.** HPA publicerar p-värdet för varje gen och cancerform (`cancer_prognostic_data.tsv`, 610 303 rader). En
  reproduktion kan kontrolleras på hundratusentals tal, inte på en figur.
- **Nästa steg är statistik, inte biologi, och det är Kims hemmaplan.** Metoden (kontrollerad mot artikelns metodavsnitt via PMC,
  PMC11683280): för varje gen prövas alla cutoffs mellan 20:e och 80:e percentilen, och den som ger lägst log-rank-p väljs. En gen
  kallas prognostisk vid p < 0,001. Ingen korrektion för multipla test nämns. Att välja den bästa cutoffen i efterhand är samma fel som
  att stoppa ett A/B-test när det ser bäst ut.

## Vad som ska reproduceras

1. p-värdet per gen och cancerform, mot HPA:s facit (HPA version 25.1). Mål: rangkorrelation nära 1 och samma antal prognostiska gener
   per cancerform inom några procent. Avvikelser redovisas.
2. Antalet prognostiska gener per cancerform (se `repro/ut/00_facit_per_cancer.csv`).

## Det enda steget till (OPRÖVAT, en hypotes tills det är räknat och en forskare har sett det)

**Hur många av de prognostiska generna överlever när p-värdet kalibreras för cutoff-sökningen?** Permutera överlevnadstiderna inom
varje cancerform, kör samma min-p-procedur, och läs av vad p < 0,001 betyder under nollhypotesen. Sedan: validerar den överlevande
delmängden bättre i HPA:s egna valideringskohorter än den bortfallna? Det sista är det som gör svaret användbart och inte bara kritiskt.

Första läsningen av facit (ingen ny beräkning, bara HPA:s egna tal): medianen av p-värdena är 0,06 till 0,11 i de flesta cancerformer.
Ett giltigt test utan signal har median 0,5. I prostata, livmoder, melanom, skivepitelcancer i lunga och tjocktarm kallas 170 till
290 gener prognostiska, där ett giltigt test skulle ge cirka 14 av slumpen. I njurcancer (KIRC) är det 7 536 av 14 201, och där
finns uppenbart verklig signal. Hypotesen är alltså att listorna är mest brus i de svaga cancerformerna och mest signal i de starka.
Det är inte visat förrän permutationen är körd.

## Risker

- Författarnas kod är R och R saknas här. Metoden är enkel nog att skriva om i Python, men då är det en omimplementering, och
  skillnader mot facit kan bero på min kod.
- Överlevnadsdata ligger inte i HPA:s filer. De hämtas från TCGA (GDC eller UCSC Xena), och HPA kan ha använt en annan version.
- Resultatet är i praktiken en kritik av en Stockholmsgrupps metod. Det ska presenteras som "den robusta delmängden av era gener och hur
  den validerar", inte som ett avslöjande. Optimerade cutoffs är ett känt problem i litteraturen (allmän kunskap, källa ska läggas till).
