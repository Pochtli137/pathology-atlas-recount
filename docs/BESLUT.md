# Beslut

Append-only. Datum, vem, varför.

## 2026-09-21 (Kim, i session)

- Projektet startat som ett prov: reproducera en publicerad beräkningsstudie på öppna data, utöka ett steg, visa för en forskare.
- Git: lokalt nu, publikt på GitHub när det finns något att visa.
- Det här är den enda spekulativa satsningen (WIP-tak: en).

## 2026-09-21 (Claude, i session, enligt kriterierna)

- Vald studie: Human Pathology Atlas v2 (Yuan … Mardinoglu, eBioMedicine 2024). Skäl och risker i `02_VAL.md`. Reserv: DrEval.
- Reproduktionen höll samma dag: 22 477 av HPA:s 22 946 prognostiska gener återfunna, rangkorrelation 0,97 till 0,9996.
- Det enda steget till: permutationskalibrering av bästa-cutoff-proceduren. Preliminärt: tröskeln p < 0,001 har 11 till 20 gånger för
  hög felfrekvens, cirka 4 500 av 23 500 gener väntas av ren slump, och i fjorton av tjugoen cancerformer går ingen enskild gen att
  säkerställa. OPRÖVAT.
- Kritiken är inte ny: Gilis m.fl. 2020 (Gent, preprint) gjorde den mot HPA version 1 i två cancerformer. Vårt tillägg är bredden,
  talen per cancerform, en kalibrerad lista och valideringen. Det ska stå överst i allt som visas för en forskare.

## 2026-09-22, natten (Claude ensam, Kim sov; mandat: "bygg något som kommer mänskligheten till nytta")

- Litteraturen läst i original där det gick (`04_LITTERATUR.md`). Bidraget är mindre nytt än det såg ut: Gilis mätte felfrekvensen i
  bröstcancer 2020, och Lausen & Schumachers formel från 1992 ger samma tal som våra permutationer. Skrivet in i 03 avsnitt 4.
- Fyra tillägg körda med fast frö: Cox utan cutoff (05), Lausen-formeln (08), spridning per dataset och helkohorttest (07), författarnas
  exakta cutoff-regel (logrank.py). Alla OPRÖVADE.
- Utkast till tvåsidaren på engelska: `05_NOTE_EN.md`. Ingen har sett den. Kim läser och avgör om, och till vem.
- Inget skickat, inget publicerat, ingen remote.

## 2026-09-22 (Kim, i session)

- Mottagare: Adil Mardinoglu, adilm@scilifelab.se (ur artikeln). Mailutkast i `06_MAIL_MARDINOGLU.md`, artefakt publicerad (privat).
- Granskning av en annan agent inarbetad: tre totaler förklarade, KICH struken ur överlevarna (12 409 i sju kohorter), Cox-förbehåll
  (svansar, tröskelgener, ingen stadiejustering), rubriken "survive correction for the cut-off search", valideringstalet in i talrutan,
  permutationsupplösning 0,005 utskriven. Kontrollerat: proteinatlas.org skiljer "validated" från "potential" på gensidan, så
  meningen "a user cannot see the difference" är omskriven. HPA:s överlapp i fyra av tio (GBM, KIRC, LIHC, LUAD) verifierat i artikeln.
- Repot publikt: github.com/Pochtli137/pathology-atlas-recount (MIT). Granskarens skäl: "on request" läses som misstro mot koden.
- Andra granskning (ChatGPT, på version 1): överlappar med den första. Nytt och inarbetat: urvalsregeln för patienter utskriven,
  rå-träffar mot HPA-etiketter åtskilda (1 012 våra utan HPA-etikett, 469 tvärtom; medianfilter förklarar under en tiondel),
  "13 av 21" mjukat till "passes this correction chain", Monte Carlo-intervall för READ/BLCA, Cox-antaganden. Granskarens uppgift att
  HPA filtrerar på median pTPM < 1 gick inte att verifiera: metodtexten laddas inte i sidkällan på proteinatlas.org.

## 2026-09-22, kväll (tredje granskningen, okänd modell, inarbetad)

- Granskningen är den bästa hittills. Rätt om: rubriksiffran inverterad (rättad), avsnitt 3 för tunt (200 permutationer, ingen
  korrektion för 21 kohorter), "indistinguishable from noise" för starkt, permutationsnollan testar global nolla och inte
  genspecifikt värde, gemensam korrigering för alla gener är ungefärlig, valideringen saknar osäkerhetsmått, konkreta gener saknas.
- **Testat och bekräftat delvis: BH-felet är Cox-svansen, inte korrelationen** (`repro/09_cox_svans.py`, 2,8 miljoner nolldragningar
  per kohort). Cox-LRT är kalibrerad vid 3,6e-6 bara i kohorter med cirka 90 dödsfall eller fler. Under det blåser svansen upp:
  KIRP 15x, UCEC 12x, KICH 73x, TGCT 822x. Cox-kolumnen i avsnitt 4 är inte pålitlig för KIRP, CESC, KICH, PRAD, TGCT. Inskrivet.
- Fyra exempelgener inlagda: CHEK2 (LUSC), MGMT (GBM), MYC och ERBB2 (OV), alla "potential prognostic" i HPA, q 0,33 till 0,48.
- 2 000 permutationer körs på de åtta kohorter där svaret kan ändras (p 0,01 till 0,20). Alla 21 hade tagit fem timmar.
  Cox-delen hoppas över, den är 200 gånger dyrare. Bonferroni och BH över 21 kohorter räknas när körningen är klar.
- Inte gjort: PFI i stället för OS, justering för ålder och stadium, genspecifik korrigering. Står som "not done" i noten.
- Granskarens sista råd följs redan: noten går till HPA-teamet, inte till publicering.

## 2026-09-23, kväll (Claude, i session, Kims "kör")

- 2 000-permutationskörningen (klar 10:41) inarbetad. BLCA hamnar inom slumpen (p 0,062), READ över (0,042). Elva av 21 kohorter
  inom slumpen, inte tio: rubriktalet, noten, mailet, figur 2 och artefakten rättade.
- Korrektion över 21 kohorter i nytt script `repro/10_kohorttest.py`: BH 5 % behåller sju (KIRC, KIRP, LIHC, LUAD, KICH, CESC, PAAD).
  HNSC, BRCA och READ faller. Bonferroni går inte att avgöra för de sex starkaste med 200 permutationer; står som "not done" i noten.
- Artefaktlänken insatt i mailutkastet. Inget skickat.

## 2026-09-23, natt (fjärde granskningen, Kim klistrade in)

- Granskningens viktigaste punkt testad: stadiumstratifierad Cox (`repro/11_stadium.py`, kalibrering `12`). LUAD 811 → 4, KIRP 805 → 0,
  PAAD 1 711 → 5, medan KIRC (7 806 → 4 203) och LIHC (2 377 → 1 117) står kvar. Stratifieringen är kalibrerad i de kohorterna.
  Detta ändrar "Follows" i noten: ojusterad association i LUAD, KIRP och PAAD går inte att skilja från stadium.
- Stratifierad Cox tillagd i `cox.py`, testad mot brute force (statsmodels saknas, de testerna hoppas över).
- Noten är inte omskriven än. Övriga punkter i granskningen bedömda i sessionen, se svaret till Kim.
