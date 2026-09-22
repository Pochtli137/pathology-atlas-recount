# Resultat, dag ett och natten efter (2026-09-21 till 22)

Status: reproduktionen är klar, utökningen är räknad med 400 000 nolldragningar per cancerform. Allt i avsnitt 2 och 3 är **oprövat**
tills en forskare har sett det.

## 1. Reproduktionen: håller

Studien: Yuan m.fl., Human Pathology Atlas v2, eBioMedicine 2024. Metoden skrevs om i Python ur artikelns metodavsnitt
(`repro/logrank.py`, testad mot `scipy.stats.logrank` och mot en rå sökning över alla cutoffs) och kördes mot HPA:s egna data
(`repro/02_reproducera.py`). Facit är HPA:s publicerade p-värde per gen och cancerform, version 25.1.

| Mått | Vårt | HPA:s |
|---|---:|---:|
| Patienter | 6 912 | 6 918 (artikeln) |
| Prognostiska gener (p < 0,001), alla 21 cancerformer | 23 489 | 22 946 |
| Varav samma gen i samma cancerform | 22 477 (98 % av HPA:s) | |
| Rangkorrelation mellan våra och HPA:s p-värden | 0,97 till 0,9996 i 20 cancerformer, 0,90 i ändtarmscancer | |
| Samma riktning (gynnsam eller ogynnsam) | 93 till 99,8 % | |
| Körtid, alla 21 | 113 sekunder | |

Avvikelserna är störst där antalet dödsfall är litet (ändtarm 16, tjocktarm 61). **Rättat natten 21–22/9:** artikeln anger TCGA-CDR och
overall survival, alltså samma överlevnadstabell som vår (`docs/04_LITTERATUR.md` avsnitt 3). Författarnas R-kod är läst
(github.com/cellur-m/pathology_atlas, `toolbox_pathology.R`, funktionen generateKMplot) men inte körd: de tar cutoffs strikt över 20:e och
till och med 80:e percentilen, hög grupp = uttryck större än eller lika med cutoffen, ingen minsta gruppstorlek. Den regeln finns nu som
`rule="authors"` i `repro/logrank.py`. Den ändrar inget: rangkorrelationen i ändtarm går från 0,9005 till 0,9023, i tjocktarm från 0,9735
till 0,9751. Avvikelsen sitter alltså inte i cutoff-regeln. Kvar som förklaring är patienturvalet, och vilka patienter HPA använde per
cancerform är OKÄNT.

En iakttagelse vi inte kan förklara: i HPA:s urinblåsekohort har alla 169 patienter avlidit enligt överlevnadstabellen. TCGA har över
400 patienter med urinblåsecancer. Hur urvalet gjordes är OKÄNT.

## 2. Det enda steget till: vad betyder p < 0,001 när cutoffen väljs för att minimera p?

Metod (`repro/03_nollfordelning.py`): patienternas överlevnad kastas om slumpmässigt mot genuttrycket, så att inget samband kan
finnas, och exakt samma procedur körs. Det ger den verkliga felfrekvensen för HPA:s tröskel, per cancerform.
**Slutliga tal: 400 000 nolldragningar per cancerform** (400 omkastningar gånger 1 000 slumpvalda gener, fast slumpfrö). En första
körning med 40 000 gav samma bild.

- Utan något samband alls är medianen av p-värdena 0,085 till 0,125. Ett giltigt test har median 0,5. I HPA:s publicerade tabell är
  medianen 0,06 till 0,11 i de flesta cancerformer, alltså nära det rena bruset.
- Tröskeln p < 0,001 har en verklig felfrekvens på **1,25 till 1,73 %, alltså 12 till 17 gånger** den nominella.
- Ren slump ger därmed **175 till 250 "prognostiska" gener per cancerform**, sammanlagt cirka 4 400 av de 23 506 (19 %).

| Grupp | Cancerformer | Kallade prognostiska | Andel som väntas av ren slump |
|---|---|---:|---:|
| Tydlig signal | njure KIRC, lever LIHC, njure KIRP | 7 563, 3 449, 1 811 | 3 till 11 % |
| Signal i stort, osäkert gen för gen | njure KICH, lunga LUAD, bukspottkörtel PAAD, livmoderhals CESC, huvud-hals HNSC | 893 till 1 507 | 14 till 24 % |
| Svårt att skilja från brus | bröst, ändtarm, urinblåsa, sköldkörtel, äggstock, hjärna GBM, magsäck, tjocktarm, melanom | 253 till 724 | 35 till 74 % |
| I nivå med brus | lunga LUSC, livmoder UCEC, prostata PRAD, testikel TGCT | 286, 227, 177, 22 | 77 till 100 % |

Efter kalibrering och kontroll av andelen falska fynd (Benjamini-Hochberg, 5 %) återstår **12 504 av 23 506 gener**, i åtta
cancerformer: KIRC 7 563, LIHC 3 428, KIRP 1 058, LUAD 326, KICH 95, HNSC 25, CESC 7 och BRCA 2. **I tretton av tjugoen cancerformer
går ingen enskild gen att säkerställa.** Vid 1 % återstår 8 190. Bukspottkörtelcancer är talande: 1 486 gener mot 232 väntade av
slump, så signal finns i stort, men ingen enskild gen går att peka ut på 5 %-nivån.

Flera kohorter har så få dödsfall att inget test kan bära slutsatser: testikel 4, prostata 9, kromofob njurcancer 9, sköldkörtel 16,
ändtarm 16. Kromofob njurcancer har ändå 1 203 prognostiska gener i HPA:s tabell, på nio dödsfall bland 64 patienter. Att 95 av dem
klarar kalibreringen här visar mest att även permutationstestet vilar på nio händelser: de ska inte tas som fynd.

## 3. Är det användbart och inte bara kritiskt? De robusta generna validerar bättre

HPA:s egen valideringstabell (oberoende kohorter, samma metod) används som utfall (`repro/04_validering.py`): är genen prognostisk
även där (p < 0,001), åt samma håll? Basnivån är vad en slumpvis vald gen uppnår.

| Cancerform | Robusta (klarar 5 %) | Validerar | Sköra | Validerar | Inte prognostiska | Basnivå |
|---|---:|---:|---:|---:|---:|---:|
| Lunga LUAD | 316 | **21,8 %** | 1 146 | 12,2 % | 4,2 % | 2,9 % |
| Lever LIHC | 3 043 | **34,6 %** | 18 | 16,7 % | 9,2 % | 8,4 % |
| Njure KIRC | 7 484 | **25,5 %** | 0 | | 5,2 % | 8,3 % |
| Hjärna GBM | 0 | | 364 | 19,2 % | 8,9 % | 5,8 % |
| Bukspottkörtel PAAD | 0 | | 1 404 | 9,8 % | 4,4 % | 3,7 % |
| Äggstock OV | 0 | | 333 | 5,4 % | 3,5 % | 2,2 % |
| Ändtarm READ | 0 | | 566 | 2,7 % | 2,2 % | 1,7 % |
| Tjocktarm COAD | 0 | | 334 | 2,7 % | 1,5 % | 1,1 % |
| Bröst BRCA | 1 | | 570 | 1,9 % | 1,3 % | 1,0 % |
| Lunga LUSC | 0 | | 280 | 1,8 % | 1,5 % | 1,0 % |

- Där båda grupperna finns (lungadenokarcinom) validerar de robusta generna nästan dubbelt så ofta som de sköra: 21,8 mot 12,2 %.
- I skivepitelcancer i lunga, bröst, tjocktarm och ändtarm validerar de "prognostiska" generna i 1,8 till 2,7 % av fallen. Gener som
  HPA INTE kallar prognostiska validerar i 1,3 till 2,2 %. Etiketten tillför alltså nästan ingenting där, vilket är vad bruskalkylen
  förutsäger.
- Förbehåll: valideringskohorterna är analyserade med samma uppblåsta procedur, så även "validerad" är för generöst. Jämförelsen mellan
  grupperna håller ändå, eftersom felet är detsamma i alla grupper.

## 4. Det här är inte nytt i sak, och det måste sägas

- Att välja cutoff efter lägsta p blåser upp felfrekvensen är känt sedan länge: Altman m.fl., "Dangers of using 'optimal' cutpoints",
  JNCI 1994, och korrektionen hos Lausen och Schumacher 1992 (källor funna via sökning, inte lästa i original).
- **Kritiken är redan riktad mot HPA:** Gilis, Taelman, Davey, Martens och Clement (Gent), "Pitfalls in re-analysis of observational omics
  studies: a post-mortem of the human pathology atlas", bioRxiv 2020, doi 10.1101/2020.03.16.994038 (läst via sammanfattande hämtning).
  De granskade version 1 i två cancerformer (bröst och lever), visade med omkastade data att proceduren ger p ner till 2,6e-7 utan
  något samband, och fann inga säkerställda gener i bröstcancer efter korrektion. Preprinten ser inte ut att ha publicerats i tidskrift.
- HPA version 2 (2024) använder samma procedur, och proteinatlas.org visar etiketterna per gen.

**Rättat natten 21–22/9 efter läsning i original (`docs/04_LITTERATUR.md`):**

- Gilis m.fl. mätte felfrekvensen i bröstcancer: en omkastning gav 320 av 17 040 gener med p < 0,001, alltså 1,88 % eller 19 gånger
  nominellt. Meningen "Gilis m.fl. gjorde inget av det" som stod här var fel. Med Efrons empiriska nollfördelning fick de noll gener i
  både bröst och lever. Deras invändning, att permutation inte rättar för korrelation mellan gener eller omätta confounders, träffar
  också vår kalibrering. Våra "robusta" gener ska därför läsas som en övre gräns.
- **Felfrekvensen går att räkna ut med en formel från 1992** (Lausen och Schumacher, efter Miller och Siegmund 1982, finns i R-paketet
  maxstat). `repro/08_formel.py`: för HPA:s intervall betyder p = 0,001 i verkligheten 0,0169, alltså 17 gånger. Våra permutationer gav
  0,74 till 1,02 gånger formelns värde, lägst i kohorter med få dödsfall. Med formeln och Benjamini-Hochberg 5 % återstår 12 451 gener,
  med permutationerna 12 504. Permutationerna behövdes alltså inte: den minsta rättelsen för HPA är ett funktionsanrop.
- TCGA-CDR, den överlevnadstabell HPA själva använder, avråder från overall survival i testikelcancer och manar till försiktighet i
  KICH, PRAD, READ, THCA och BRCA (Liu m.fl., Cell 2018, tabell 3). Poängen om för få dödsfall är redan gjord av datakällan.
- HPA version 2 nämner varken Gilis, Altman, multipel testning eller cutoff-problemet. Proteinatlas.org version 25.1 har ingen
  varningstext. Korrekt referens: Yuan m.fl., eBioMedicine 2025;111:105495.

**Det som återstår som vårt:** alla 21 cancerformer i version 2 med tal per cancerform, visningen att formeln från 1992 räcker,
jämförelsen mot ett test utan cutoff (avsnitt 5), spridningen i antalet slumpgener per dataset (avsnitt 6) och valideringsjämförelsen.
Det är ett hantverksbidrag, inte en upptäckt, och mindre nytt än det såg ut i går.

## 5. Känslighetsanalys utan cutoff: Cox-regression (natten 21–22/9, OPRÖVAT)

Samma fråga ställd med ett test som inte väljer någon cutoff: Cox-regression per gen på log2(pTPM + 1), sannolikhetskvotstest,
Benjamini-Hochberg 5 % (`repro/cox.py`, `repro/05_cox.py`). Koden är testad mot statsmodels PHReg (Breslow) och mot en rå optimering
(`repro/test_cox.py`). R:s coxph använder Efron som standard, skillnaden är OKÄND här men liten vid tider i dagar.

**A. Testet håller sin felfrekvens.** Samma omkastning som i avsnitt 2, 100 000 nolldragningar per cancerform: medianen av p är 0,49 till
0,52 (ska vara 0,5), felfrekvensen vid 0,05 är 4,4 till 5,9 %, och vid 0,001 är den 0,06 till 0,17 % i de sjutton cancerformer som har
minst 28 dödsfall. HPA:s procedur ligger på 1,25 till 1,73 %. Undantag: testikel (4 dödsfall) och sköldkörtel (16) ligger på 0,33 %,
där bär inget test. Se `docs/figurer/fig1_calibration.png`.

**B. Bilden blir densamma, med ett viktigt undantag.**

| Cancerform | Dödsfall | HPA:s etikett | Klarar kalibreringen (avsnitt 2) | Cox, 5 % falska fynd |
|---|---:|---:|---:|---:|
| Njure KIRC | 171 | 7 563 | 7 563 | 7 853 |
| Lever LIHC | 129 | 3 449 | 3 428 | 3 482 |
| Bukspottkörtel PAAD | 92 | 1 486 | 0 | **2 380** |
| Njure KIRP | 44 | 1 811 | 1 058 | 926 |
| Lunga LUAD | 180 | 1 507 | 326 | 750 |
| Huvud-hals HNSC | 212 | 893 | 25 | 170 |
| Livmoderhals CESC | 67 | 927 | 7 | 55 |
| Övriga fjorton | | 5 870 | 97 | 194 |
| **Summa** | | **23 506** | **12 504** | **15 810** |

- I tolv cancerformer hittar Cox högst fem gener, i sju av dem ingen alls: lunga LUSC, äggstock, ändtarm, melanom, magsäck, tjocktarm
  och sköldkörtel. HPA anger 253 till 567 prognostiska gener i var och en av dem.
- **Bukspottkörtelcancer:** kalibreringen i avsnitt 2 kunde inte peka ut någon enskild gen, Cox pekar ut 2 380. Avsnitt 2 sa "signal i
  stort", och det här visar att den går att lokalisera med ett test som har bättre styrka. Det kalibrerade p-värdet är alltså giltigt
  men svagt: att rädda HPA:s procedur med kalibrering kostar styrka jämfört med att byta test.
- Av HPA:s 23 506 gener bekräftas 12 519 av Cox. Cox hittar 3 291 som HPA inte har. Riktningen (gynnsam eller ogynnsam) är densamma
  i 99 till 100 % av HPA:s gener.
- Räkna inte testikel (22 gener på 4 dödsfall, modellen konvergerar inte för 22 gener), kromofob njurcancer (154 på 9 dödsfall) och
  prostata (2 på 9) som fynd.

**C. Validering vid samma listlängd** (`repro/ut/05_validering.csv`, `docs/figurer/fig3_validation.png`): ta lika många gener som HPA
kallar prognostiska, men rangordna efter Cox-p. Andelen som är prognostisk åt samma håll i HPA:s valideringskohort är lika eller
högre i åtta av tio cancerformer (glioblastom 29,4 mot 19,2 %, LUSC 3,2 mot 1,8 %, LUAD 15,6 mot 14,3 %), lägre i bukspottkörtel (8,5
mot 9,8 %) och marginellt lägre i lever (34,2 mot 34,5 %). Skillnaderna är små utom i glioblastom, och utfallet är definierat med HPA:s egen procedur, vilket
gynnar HPA:s lista. Slutsats: det giltiga testet förlorar inget i validering.

## 6. Hur många slumpgener ger ETT dataset? (natten 21–22/9, OPRÖVAT)

Avsnitt 2 ger medelvärdet. Gener är korrelerade, så antalet i ett enskilt dataset sprider mycket mer än binomialt. 200 omkastningar per
cancerform på alla gener (`repro/07_antal_under_noll.py`, `repro/ut/07_sammanfattning.csv`), figur 2 omritad.

- Ett typiskt omkastat dataset ger 100 till 190 "prognostiska" gener (median). Vart tjugonde ger 450 till 870. Värsta av 200: sköldkörtel
  3 900, njure KIRC 3 163, bukspottkörtel 2 676.
- **I tio av tjugoen cancerformer ligger HPA:s antal prognostiska gener inom vad slumpen ger** (ensidigt permutationstest på antalet,
  p > 0,05): COAD, GBM, LUSC, OV, PRAD, SKCM, STAD, TGCT, THCA, UCEC. Ändtarm (0,050) och urinblåsa (0,055) på gränsen. Nio bär signal
  (p ≤ 0,025), bröst svagast (724 mot 95:e percentilen 490). Det här testet är giltigt oavsett korrelation mellan gener.
- Samma korrelation slår mot varje per-gen-korrektion: Cox med Benjamini-Hochberg 5 % gav under omkastning minst ett "fynd" i 1,5 till
  13,5 % av dataseten där dödsfallen räcker, och då ofta hundratals eller tusentals (HNSC upp till 2 299, KIRC 3 534). Det är Gilis
  invändning, och den gäller alla gental i avsnitt 2 och 5. Talen är övre gränser.

## 7. Kvar innan något visas för en forskare

1. ~~Körningen med 400 000 dragningar~~ klar 21/9.
2. ~~Läs Gilis m.fl. i original, sök fler kritiker~~ `docs/04_LITTERATUR.md`. Altman 1994 och Lausen & Schumacher 1992 kvar (betalvägg).
3. Författarnas R-kod läst, cutoff-regeln utesluten som förklaring (avsnitt 1). Kvar: patienturvalet. Kräver R eller kontakt.
4. ~~Räkna om basnivån~~ klar 21/9. ~~Känslighetsanalys med Cox-regression~~ klar natten 21–22/9, avsnitt 5.
5. ~~Skriv två sidor på engelska med tre figurer~~ utkast `docs/05_NOTE_EN.md` natten 21–22/9, Kim läser före allt annat. Kim avgör vem som får dem: Adil Mardinoglu (korresponderande författare, KTH och
   SciLifeLab) är den naturliga, Lieven Clement i Gent den som redan bryr sig om frågan.
