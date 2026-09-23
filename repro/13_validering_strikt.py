#!/usr/bin/env python3
"""Validation with an outcome that does not share the false-positive rate of the procedure being tested (fourth review).
05 counts a gene as validated when HPA's validation cohort gives best cut-off p < 0.001 in the same direction, so the base rate (1.0-1.7 %)
is largely the procedure's own false-positive rate. Here the validation p is first corrected for the cut-off search with Lausen and
Schumacher's formula (08), and a gene validates when the corrected p < 0.001 in the same direction. Lists of equal size: HPA's (ranked by
our best cut-off p), Cox, and Cox stratified on stage (11). Lift = validated share / base rate.
Input: repro/ut/05_<ABBR>_gener.csv, 11_<ABBR>_gener.csv, data/hpa/cancer_prognostic_data.tsv. No randomness. Output: repro/ut/13_validering.csv"""
import importlib.util
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("rep", HERE / "02_reproducera.py"); rep = importlib.util.module_from_spec(spec); spec.loader.exec_module(rep)
src = (HERE / "08_formel.py").read_text(); ns = {"__file__": str(HERE / "08_formel.py")}; exec(src.split("\ndef bh")[0], ns); lausen = ns["p_lausen92"]; UT = rep.UT

fac = pd.read_csv(rep.ROOT / "data" / "hpa" / "cancer_prognostic_data.tsv", sep="\t"); kol = list(fac.columns[3:]); fac["p"] = fac[kol].bfill(axis=1).iloc[:, 0]
fac["dir"] = np.where(fac[kol].notna().idxmax(axis=1).str.contains("unfavorable"), 1, -1); fac = fac.dropna(subset=["p"]); rows = []
for a, namn in rep.NAMN.items():
    v = fac[fac.Cancer == namn + " (validation)"]
    if v.empty: continue
    d = pd.read_csv(UT / f"05_{a}_gener.csv").dropna(subset=["p_ours"]).merge(v[["Gene", "p", "dir"]].rename(columns={"p": "p_val", "dir": "dir_val"}), on="Gene")
    f11 = UT / f"11_{a}_gener.csv"; d = d.merge(pd.read_csv(f11), on="Gene", how="left") if f11.exists() else d.assign(cox_p_stage_strata=np.nan)
    d["pc_val"] = lausen(d.p_val.to_numpy()); n = int((d.p_ours < 1e-3).sum())
    for navn, ok in (("raw", d.p_val < 1e-3), ("corrected", d.pc_val < 1e-3)):
        base = float(np.mean([(ok & (d.dir_val == s)).mean() for s in (1, -1)]))
        val = lambda x, col: float((ok[x.index] & (x.dir_val == x[col])).mean()) if len(x) else np.nan
        r = dict(cohort=a, outcome=navn, list_size=n, base_rate=round(base, 4), hpa_list=round(val(d.nsmallest(n, "p_ours"), "minp_direction"), 4), cox_list=round(val(d.nsmallest(n, "cox_p"), "cox_direction"), 4))
        if d.cox_p_stage_strata.notna().any():
            # direction for the stratified model is not stored; use the unstratified Cox direction (changes sign only for near-null genes)
            r["cox_stage_list"] = round(val(d.dropna(subset=["cox_p_stage_strata"]).nsmallest(n, "cox_p_stage_strata"), "cox_direction"), 4)
        r.update(hpa_lift=round(r["hpa_list"] / base, 1) if base else np.nan, cox_lift=round(r["cox_list"] / base, 1) if base else np.nan); rows.append(r)
out = pd.DataFrame(rows); out.to_csv(UT / "13_validering.csv", index=False); print(out.to_string(index=False))
