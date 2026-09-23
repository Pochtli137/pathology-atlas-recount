#!/usr/bin/env python3
"""Is the Cox signal stage? Univariate Cox (as 05) against Cox stratified on tumour stage, and on stage x grade, per cohort.
Stage from TCGA-CDR ajcc_pathologic_tumor_stage, collapsed to I-IV (clinical_stage, FIGO, when ajcc is mostly missing, as in CESC).
Grade from histological_grade (G1-G4). Patients without stage (or grade) are dropped, so the unstratified model is also rerun on
exactly the patients the stratified model uses: the difference between those two columns is the effect of stage alone.
Counts: genes with BH q < 0.05, and how many of the genes our best cut-off run labels prognostic (p < 0.001) keep Cox q < 0.05.
No randomness. Output: repro/ut/11_stadium.csv"""
import sys, re, importlib.util
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from logrank import prepare  # noqa: E402
from cox import cox_for_genes  # noqa: E402
spec = importlib.util.spec_from_file_location("rep", HERE / "02_reproducera.py"); rep = importlib.util.module_from_spec(spec); spec.loader.exec_module(rep)
UT = rep.UT; ROOT = HERE.parent
KOHORTER = sys.argv[1:] or list(rep.NAMN)


def bh(p):
    n = len(p); o = np.argsort(p); q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]; return np.minimum(q, 1)


def stadium(x):
    m = re.match(r"Stage (IV|III|II|I)", str(x)); return m.group(1) if m else None


def ladda_klinik(abbr):
    """Same patient selection as rep.ladda, plus the clinical columns."""
    c = rep.NAMN[abbr] + " (TCGA)"; z = np.load(ROOT / "data" / "hpa" / "per_cancer" / f"{rep.slug(c)}.npz", allow_pickle=False)
    genes, samples, m = z["genes"], z["samples"], z["ptpm"]
    surv = pd.read_csv(ROOT / "data" / "xena" / "Survival_SupplementalTable_S1_20171025_xena_sp", sep="\t").drop_duplicates("_PATIENT").set_index("_PATIENT")
    pat = pd.Series([s[:12] for s in samples]); keep = ~pat.duplicated() & pat.isin(surv.index)
    s = surv.loc[pat[keep]]; ok = (s["OS"].notna() & s["OS.time"].notna() & (s["OS.time"] > 0)).to_numpy()
    return genes, m[:, np.flatnonzero(keep.to_numpy())[ok]], s[ok]


def kor(expr, tid, ev, grupp=None):
    prep = prepare(tid, ev)
    st = None if grupp is None else [(np.flatnonzero(grupp == g), prepare(tid[grupp == g], ev[grupp == g])) for g in np.unique(grupp)]
    return cox_for_genes(expr, prep, strata=st)["p_lrt"]


rows = []
for abbr in KOHORTER:
    genes, m, s = ladda_klinik(abbr); keep = np.nanmean(m, axis=1) > 1; expr = np.log2(np.nan_to_num(m[keep].astype(np.float64)) + 1); g = genes[keep]
    tid = s["OS.time"].to_numpy(float); ev = s["OS"].to_numpy(int).astype(bool)
    stg = s["ajcc_pathologic_tumor_stage"].map(stadium); kalla = "ajcc"
    if stg.isna().mean() > 0.5: stg = s["clinical_stage"].map(stadium); kalla = "clinical"
    grd = s["histological_grade"].where(s["histological_grade"].isin(["G1", "G2", "G3", "G4"]))
    labels = pd.read_csv(UT / f"02_{abbr}_gener.csv").set_index("Gene")["p_ours"].reindex(g).to_numpy() < 1e-3
    if stg.notna().mean() < 0.6: print(abbr, "no usable stage, skipped", flush=True); continue
    har = stg.notna().to_numpy(); X = expr[:, har]; t, e = tid[har], ev[har]; sg = stg.to_numpy()[har]
    p_all = kor(expr, tid, ev); p_sub = kor(X, t, e); p_st = kor(X, t, e, sg)
    r = dict(cohort=abbr, stage_source=kalla, patients=len(tid), events=int(ev.sum()), patients_staged=int(har.sum()), events_staged=int(e.sum()),
             stage_counts=" ".join(f"{k}:{v}" for k, v in pd.Series(sg).value_counts().sort_index().items()),
             cox_all=int((bh(p_all) < .05).sum()), cox_staged_patients=int((bh(p_sub) < .05).sum()), cox_stage_strata=int((bh(p_st) < .05).sum()),
             labelled=int(labels.sum()), labelled_keep_unadj=int((labels[:] & (bh(p_sub) < .05)).sum()), labelled_keep_stage=int((labels & (bh(p_st) < .05)).sum()))
    hg = har & grd.notna().to_numpy()
    if hg.sum() > 0.6 * len(tid):
        sgg = (stg.astype(str) + "_" + grd.astype(str)).to_numpy()[hg]; p_sub2 = kor(expr[:, hg], tid[hg], ev[hg]); p_sg = kor(expr[:, hg], tid[hg], ev[hg], sgg)
        r.update(patients_stage_grade=int(hg.sum()), cox_sg_patients=int((bh(p_sub2) < .05).sum()), cox_stage_grade_strata=int((bh(p_sg) < .05).sum()))
    rows.append(r); print(r, flush=True)
    pd.DataFrame(dict(Gene=g, cox_p_staged_patients=p_sub, cox_p_stage_strata=p_st)).to_csv(UT / f"11_{abbr}_gener.csv", index=False)
out = pd.DataFrame(rows); out.to_csv(UT / "11_stadium.csv", index=False); print(out.to_string(index=False))
