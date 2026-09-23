# Mail till Adil Mardinoglu (utkast v3, 23/9 kväll, Kim skriver om och skickar själv)

Mottagare: adilm@scilifelab.se (korresponderande författare, PMC11683280). Ramen är Kims (22/9): testet är om en vanlig människa
med en AI kan göra forskning utan att förstå fältet, och det ska sägas rakt, inte som ursäkt.

---

Subject: Pathology Atlas v2: in LUAD, KIRP and PAAD the prognostic label seems to track stage (an AI-assisted reproduction)

Dear Adil,

I need to start with a disclaimer, because it is the point of this mail.

I am not a researcher and I do not understand the biology in what I am about to send you. I work with A/B testing and statistics
in e-commerce. I wanted to test one question: can an ordinary person, using an AI agent, produce research work that a researcher
finds useful, without understanding the field? To test it I asked the agent to reproduce a recent computational study on open data
and extend it one step. It chose your Human Pathology Atlas v2 paper.

What it did: reproduced your prognostic gene analysis from HPA's public expression file and the TCGA-CDR survival table in all 21
TCGA cohorts (98 % of your 22,946 prognostic gene–cohort pairs come back), then made two corrections the atlas does not make.

For the search over cut-offs: the real false-positive rate at p < 0.001 is 1.25 to 1.73 %, as Lausen and Schumacher's formula
predicts. Gilis et al. made this point against v1 in 2020, so that part is not new.

For stage: a Cox model stratified on AJCC stage keeps about half the signal in KIRC and LIHC, but in lung adenocarcinoma it goes from
811 significant genes to 4, in papillary kidney from 805 to 0, and in pancreas from 1,711 to 5. The stratified test is calibrated.
We have not seen this reported.

So the question is sharper than whether this is useful: why does proteinatlas.org show an uncorrected, unadjusted label when the
corrected p-value is one function call and stage is in the same table? If there is a good reason, or a flaw in what we did, I would
like to know which.

The note with two figures is here: https://claude.ai/artifact/ACSJ9kjTutmbP5FKxvnqKV. The code and per-gene tables are at
github.com/Pochtli137/pathology-atlas-recount, every number from a numbered script with a fixed seed.

Four AI reviews have criticised the draft and the agent worked their points in, but no human who knows the field has read it. I
cannot judge whether it is sound myself, and any answer, including "this is known" or "this is wrong", is a result for my experiment.

Thank you for your time, and apologies if this lands as noise in a busy inbox.

Best regards,
Kim Dahlroth
Stockholm
kim.dahlroth@gmail.com

---

Anteckningar:
- Rubriken säger vad det är innan han öppnar. Ingen forskare vill bli lurad att läsa fem stycken innan disclaimern.
- Artefakten är redan delad med länk (version 5, 23/9).
- Svarar han inte inom två veckor: Lieven Clement (Gent), som redan bryr sig om frågan.
