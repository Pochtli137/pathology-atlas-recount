#!/usr/bin/env python3
"""Do the genes that survive calibration validate better? Outcome = HPA's own validation table (independent cohorts, same best-cutoff
method): prognostic there too (p < 0.001) AND same direction. Base rate = share of all tested genes that are prognostic in the validation
cohort with a given direction, averaged over the two directions (what a random gene list would achieve).
Input: repro/ut/03_<ABBR>_gener.csv, data/hpa/cancer_prognostic_data.tsv. Output: repro/ut/04_validering.csv"""
import importlib.util
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent; spec = importlib.util.spec_from_file_location("rep", HERE / "02_reproducera.py"); rep = importlib.util.module_from_spec(spec); spec.loader.exec_module(rep)
fac = pd.read_csv(rep.ROOT / "data" / "hpa" / "cancer_prognostic_data.tsv", sep="\t"); kol = list(fac.columns[3:]); fac["p"] = fac[kol].bfill(axis=1).iloc[:, 0]
fac["dir"] = np.where(fac[kol].notna().idxmax(axis=1).str.contains("unfavorable"), 1, -1); fac = fac.dropna(subset=["p"])
rows = []
for a, namn in rep.NAMN.items():
    v = fac[fac.Cancer == namn + " (validation)"]
    if v.empty: continue
    d = pd.read_csv(rep.UT / f"03_{a}_gener.csv").merge(v[["Gene", "p", "dir"]].rename(columns={"p": "p_val", "dir": "dir_val"}), on="Gene")
    d["validated"] = (d.p_val < 1e-3) & (d.dir_val == d.direction); pr = d[d.p_ours < 1e-3]; rob, fra = pr[pr.q < 0.05], pr[pr.q >= 0.05]
    base = np.mean([((d.p_val < 1e-3) & (d.dir_val == s)).mean() for s in (1, -1)])
    rows.append(dict(cohort=a, genes=len(d), prognostic=len(pr), robust=len(rob), robust_validated=round(rob.validated.mean(), 3) if len(rob) else np.nan,
                     fragile=len(fra), fragile_validated=round(fra.validated.mean(), 3) if len(fra) else np.nan, not_prognostic_validated=round(d[d.p_ours >= 1e-3].validated.mean(), 3), base_rate=round(base, 3)))
out = pd.DataFrame(rows); out.to_csv(rep.UT / "04_validering.csv", index=False); print(out.to_string(index=False))
