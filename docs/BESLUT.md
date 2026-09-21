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
