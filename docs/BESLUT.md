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
