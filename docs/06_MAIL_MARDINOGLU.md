# Mail till Adil Mardinoglu (utkast v2, Kim skriver om och skickar själv)

Mottagare: adilm@scilifelab.se (korresponderande författare, PMC11683280). Ramen är Kims (22/9): testet är om en vanlig människa
med en AI kan göra forskning utan att förstå fältet, och det ska sägas rakt, inte som ursäkt.

---

Subject: An experiment: an AI reproduced your Pathology Atlas prognostic genes. Is the result useful to you?

Dear Adil,

I need to start with a disclaimer, because it is the point of this mail.

I am not a researcher and I do not understand the biology in what I am about to send you. I work with A/B testing and statistics
in e-commerce. I wanted to test one question: can an ordinary person, using an AI agent, produce research work that a researcher
finds useful, without understanding the field? To test it I asked the agent to reproduce a recent computational study on open data
and extend it one step. It chose your Human Pathology Atlas v2 paper.

What it did: reproduced your prognostic gene analysis from HPA's public expression file and the TCGA-CDR survival table in all 21
TCGA cohorts (98 % of your 22,946 prognostic gene–cohort pairs come back). Then it asked what p < 0.001 means once the best cut-off
search is accounted for: under permutation the real false-positive rate is 1.25 to 1.73 % in every cohort, which matches Lausen and
Schumacher's 1992 formula, and in ten cohorts the number of labelled genes is within what permutation alone produces. It also found
that Gilis et al. raised the same point against v1 in 2020, so the statistics are not new; the coverage across v2 is.

The two-page note with three figures is here: [artefaktlänk]. The code and per-gene tables are at
github.com/Pochtli137/pathology-atlas-recount, every number from a numbered script with a fixed seed.

I had three other AI models criticise the draft and the agent worked their points in, but no human who knows the field has read
it. I cannot judge whether it is sound. So my question is simple, and either answer is a result for my experiment: is this useful to
you? If it is wrong or beside the point, where?

Thank you for your time, and apologies if this lands as noise in a busy inbox.

Best regards,
Kim Dahlroth
Stockholm
kim.dahlroth@gmail.com

---

Anteckningar:
- Rubriken säger vad det är innan han öppnar. Ingen forskare vill bli lurad att läsa fem stycken innan disclaimern.
- Dela artefakten innan länken går ut.
- Svarar han inte inom två veckor: Lieven Clement (Gent), som redan bryr sig om frågan.
