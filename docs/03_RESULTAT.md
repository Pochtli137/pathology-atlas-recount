# Resultat, dag ett (2026-09-21)

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

Avvikelserna är störst där antalet dödsfall är litet (ändtarm 16, tjocktarm 61). Trolig orsak: vi använder TCGA:s kurerade
överlevnadstabell (TCGA-CDR via UCSC Xena), och vilken tabell HPA använde är OKÄNT. Författarnas R-kod är inte körd, R saknas här.

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

**Det som är nytt här:** alla 21 cancerformer i version 2, den verkliga felfrekvensen och det väntade antalet falska fynd per cancerform,
ett kalibrerat p-värde per gen som behåller HPA:s egen procedur (så att det går att lägga bredvid deras tabell), och en första
visning av att robustheten förutsäger validering i HPA:s egna valideringskohorter. Gilis m.fl. gjorde inget av det. Det är ett
hantverksbidrag, inte en upptäckt.

## 5. Kvar innan något visas för en forskare

1. ~~Körningen med 400 000 dragningar~~ klar 21/9.
2. Läs Gilis m.fl. och Altman i original. Sök efter fler kritiker och efter svar från HPA.
3. Kör författarnas R-kod (kräver R) eller förklara avvikelserna i ändtarm och tjocktarm på annat sätt.
4. ~~Räkna om basnivån~~ klar 21/9. Lägg till en känslighetsanalys med Cox-regression på kontinuerligt uttryck.
5. Skriv två sidor på engelska med tre figurer. Kim avgör vem som får dem: Adil Mardinoglu (korresponderande författare, KTH och
   SciLifeLab) är den naturliga, Lieven Clement i Gent den som redan bryr sig om frågan.
