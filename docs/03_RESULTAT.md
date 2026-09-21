# Resultat, dag ett (2026-09-21)

Status: reproduktionen är klar. Utökningen är räknad med 40 000 nolldragningar per cancerform, en körning med 400 000 pågår och
ersätter talen i avsnitt 2 när den är klar. Allt i avsnitt 2 och 3 är **oprövat** tills en forskare har sett det.

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

**Preliminära tal (40 000 nolldragningar per cancerform):**

- Utan något samband alls är medianen av p-värdena 0,085 till 0,125. Ett giltigt test har median 0,5. I HPA:s publicerade tabell är
  medianen 0,06 till 0,11 i de flesta cancerformer, alltså nära det rena bruset.
- Tröskeln p < 0,001 har en verklig felfrekvens på 1,1 till 2,0 %, **11 till 20 gånger** den nominella.
- Ren slump ger därmed **165 till 290 "prognostiska" gener per cancerform**, cirka 4 500 av de 23 500.

| Grupp | Cancerformer | Kallade prognostiska | Väntade av ren slump |
|---|---|---:|---:|
| Tydlig signal | njure KIRC, lever LIHC, njure KIRP | 7 563, 3 449, 1 811 | 291, 167, 159 (4 till 9 %) |
| Signal i stort, osäkert gen för gen | lunga LUAD, bukspottkörtel PAAD, njure KICH, livmoderhals CESC, huvud-hals HNSC, bröst BRCA | 724 till 1 507 | 15 till 30 % |
| Svårt att skilja från brus | ändtarm, urinblåsa, sköldkörtel, tjocktarm, äggstock, hjärna GBM, melanom, livmoder, magsäck | 227 till 567 | 40 till 74 % |
| I nivå med brus | lunga LUSC, prostata PRAD, testikel TGCT | 286, 177, 22 | 268, 173, 206 (94 till 100 %) |

Efter kalibrering och kontroll av andelen falska fynd (Benjamini-Hochberg, 5 %) återstår 12 901 av 23 506 gener, alla i sju
cancerformer: KIRC 7 563, LIHC 3 449, KIRP 1 264, LUAD 409, KICH 140, CESC 51, HNSC 25. I de övriga fjorton går ingen enskild gen att
säkerställa på den nivån. Med 40 000 dragningar är upplösningen för grov för att hitta ett fåtal mycket starka gener, därför körs
400 000.

Flera kohorter har så få dödsfall att inget test kan bära slutsatser: testikel 4, prostata 9, kromofob njurcancer 9, sköldkörtel 16,
ändtarm 16. Kromofob njurcancer har ändå 1 203 prognostiska gener i HPA:s tabell, på nio dödsfall bland 64 patienter.

## 3. Är det användbart och inte bara kritiskt? De robusta generna validerar bättre

HPA:s egen valideringstabell (oberoende kohorter, samma metod) används som utfall: är genen prognostisk även där, åt samma håll?

| Cancerform | Robusta (klarar 5 %) | Validerar | Sköra | Validerar | Basnivå |
|---|---:|---:|---:|---:|---:|
| Lunga LUAD | 398 | 20,9 % | 1 064 | 11,8 % | 2,9 % |
| Njure KIRC | 7 484 | 25,5 % | 0 | | 8,3 % |
| Lever LIHC | 3 061 | 34,5 % | 0 | | 8,4 % |
| Lunga LUSC | 0 | | 280 | 1,8 % | 1,0 % |
| Bröst BRCA | 0 | | 571 | 1,9 % | 1,0 % |
| Tjocktarm COAD | 0 | | 334 | 2,7 % | 1,1 % |
| Ändtarm READ | 0 | | 566 | 2,7 % | 1,7 % |
| Bukspottkörtel PAAD | 0 | | 1 404 | 9,8 % | 3,7 % |
| Hjärna GBM | 0 | | 364 | 19,2 % | 5,8 % |

I den enda cancerform där båda grupperna finns (lungadenokarcinom) validerar de robusta generna nästan dubbelt så ofta som de sköra.
I skivepitelcancer i lunga, bröst, tjocktarm och ändtarm validerar de "prognostiska" generna knappt över basnivån, vilket är vad
bruskalkylen förutsäger. Basnivån är grov (andelen prognostiska i valideringskohorten delat med två) och ska räknas om ordentligt.

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

1. Körningen med 400 000 dragningar, och talen ovan uppdaterade.
2. Läs Gilis m.fl. och Altman i original. Sök efter fler kritiker och efter svar från HPA.
3. Kör författarnas R-kod (kräver R) eller förklara avvikelserna i ändtarm och tjocktarm på annat sätt.
4. Räkna om basnivån i avsnitt 3 ordentligt, och lägg till en känslighetsanalys med Cox-regression på kontinuerligt uttryck.
5. Skriv två sidor på engelska med tre figurer. Kim avgör vem som får dem: Adil Mardinoglu (korresponderande författare, KTH och
   SciLifeLab) är den naturliga, Lieven Clement i Gent den som redan bryr sig om frågan.
