#!/usr/bin/env python3
"""Does the stage-stratified Cox of 11 hold its level? Expression is permuted WITHIN stage strata (keeps the stage-survival link, breaks
gene-survival given stage), all genes, N_PERM times per cohort. Reports the false-positive rate at p < 0.001 and how often BH q < 0.05
returns anything. The same for the unstratified Cox on the same patients (permuted over all patients) as reference.
Usage: python3 repro/12_stadium_kalibrering.py [N_PERM] [cohorts]. Seed fixed. Output: repro/ut/12_kalibrering.csv"""
import sys, importlib.util
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
spec = importlib.util.spec_from_file_location("s11", HERE / "11_stadium.py")
src = (HERE / "11_stadium.py").read_text().split("\nrows = []")[0]; s11 = {"__file__": str(HERE / "11_stadium.py")}; exec(compile(src, "11_stadium.py", "exec"), s11)
UT = s11["UT"]; NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 20; KOH = sys.argv[2:] or ["BRCA", "CESC", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "PAAD", "TGCT"]
rows = []
for abbr in KOH:
    genes, m, s = s11["ladda_klinik"](abbr); keep = np.nanmean(m, axis=1) > 1; expr = np.log2(np.nan_to_num(m[keep].astype(np.float64)) + 1)
    stg = s["ajcc_pathologic_tumor_stage"].map(s11["stadium"])
    if stg.isna().mean() > 0.5: stg = s["clinical_stage"].map(s11["stadium"])
    har = stg.notna().to_numpy(); X = expr[:, har]; t = s["OS.time"].to_numpy(float)[har]; e = s["OS"].to_numpy(int).astype(bool)[har]; sg = stg.to_numpy()[har]
    rng = np.random.default_rng([20260923, sum(map(ord, abbr))]); fp_s, fp_u, bh_s, bh_u = [], [], [], []
    for i in range(NPERM):
        perm = np.arange(len(sg))
        for g in np.unique(sg): idx = np.flatnonzero(sg == g); perm[idx] = rng.permutation(idx)
        ps = s11["kor"](X[:, perm], t, e, sg); pu = s11["kor"](X[:, rng.permutation(len(sg))], t, e)
        fp_s.append((ps < 1e-3).mean()); fp_u.append((pu < 1e-3).mean()); bh_s.append(int((s11["bh"](ps) < .05).sum())); bh_u.append(int((s11["bh"](pu) < .05).sum()))
    r = dict(cohort=abbr, events=int(e.sum()), strata=len(np.unique(sg)), permutations=NPERM, fpr001_stage=round(float(np.mean(fp_s)), 5), fpr001_unadj=round(float(np.mean(fp_u)), 5),
             bh_any_stage=round(float(np.mean(np.array(bh_s) > 0)), 2), bh_max_stage=max(bh_s), bh_any_unadj=round(float(np.mean(np.array(bh_u) > 0)), 2), bh_max_unadj=max(bh_u))
    rows.append(r); print(r, flush=True)
pd.DataFrame(rows).to_csv(UT / "12_kalibrering.csv", index=False); print(pd.DataFrame(rows).to_string(index=False))
