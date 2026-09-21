# cancer-prov

Ett prov, inte ett projekt att bota cancer. Startat 2026-09-21 av Kim efter frågan "lets cure cancer".

## Frågan provet svarar på

**Kan Kim och ett gäng agenter producera något som en cancerforskare tycker är användbart?** Inget annat. Är svaret nej vet vi det
billigt. Är svaret ja finns något att visa den partner som riktigt arbete kräver: utan någon som kan pröva en hypotes i labb eller
klinik är allt vi tar fram brus.

## Upplägget

1. **Välj** en nyligen publicerad beräkningsstudie om cancer med öppna data och öppen kod (`docs/01_URVAL.md`, kandidater i `kandidater/`).
2. **Reproducera** huvudresultatet här, från rådata till figur (`repro/`). Avvikelser redovisas, de göms inte.
3. **Utöka ett steg.** Ett, inte fem: en fråga studien inte ställde men som dess data kan svara på.
4. **Visa** resultatet för en forskare, helst någon av studiens författare. Kim tar kontakten, aldrig Claude.
5. **Döm:** användbart eller inte, med forskarens ord. Skrivs i `docs/BESLUT.md`.

Ram: en till två veckor av Kims lediga dagar. Går tiden ut utan reproduktion är det också ett svar.

## Regler

- **Det här är Kims enda spekulativa satsning** (WIP-tak: en). Inget annat spekulativt startas medan den pågår.
- **Ingen hypotes presenteras som fynd.** Vi kan ta fram hundra rimliga påståenden på en eftermiddag och kan inte skilja de sanna
  från de falska. Allt som inte är en reproduktion av något publicerat märks "oprövat" tills en forskare sagt något om det.
- **Reproduktionen går före utökningen.** Stämmer inte talen med artikeln stannar vi där och reder ut varför.
- **Källa på varje påstående** om biologi eller om studien: artikel, tabell, figur, datasetversion. Skriv OKÄNT hellre än något rimligt.
- **Bara öppna data.** Inga patientdata som kräver ansökan, inga inloggade tjänster, inget som kräver etikprövning.
- **Inga körningar till API-pris** (`claude -p`, extra förbrukning) utan att Kim sagt ja med belopp.
- Data ligger i `data/` och är gitignorerad. Varje dataset får en rad i `data/KALLOR.md`: varifrån, version, datum, licens, kontrollsumma.
- Varje analys är ett körbart script med fast slumpfrö. Inga resultat som bara finns i en session.
- Lokalt git utan remote tills det finns något att visa. Då publikt på GitHub (Kims beslut 21/9): öppen kod på öppna data är det som
  gör resultatet granskningsbart.
- Svenska i dokumenten, engelska i kod och i det som ska visas för en forskare.

## Maskinen

Apple M3 Max, 64 GB minne, cirka 70 GB ledig disk (21/9). Python 3.13 med numpy, pandas, scipy, scikit-learn, matplotlib och torch.
R, conda och docker saknas. Urvalet ska därför gynna studier vars kod är Python, eller vars analys är enkel nog att skriva om.

## Grannprojekt

- `~/Projects/jakt/`: jobbjakten. Spår E där (uppdrag inom screening) är parkerat och är en annan fråga än den här.
- `~/Projects/Hsng/`: rör inte. Klientrepo.
