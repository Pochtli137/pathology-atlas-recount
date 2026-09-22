#!/usr/bin/env python3
"""Sensitivity analysis: the same question (is this gene associated with overall survival?) asked with a test that has no cutoff.

Per cohort and gene: univariate Cox regression on standardised log2(pTPM + 1), likelihood-ratio p, Benjamini-Hochberg q (cox.py).
Three things are reported.
  A. Validity: the same permutation scheme as 03 (patients permuted against expression), so the Cox test's real false-positive rate at
     nominal 0.001 and 0.05 can be laid next to the best-cutoff procedure's. A valid test gives 0.001 and 0.05.
  B. Agreement: how many genes Cox finds at q < 0.05, and the overlap with HPA's label (min-p < 0.001) and with the permutation-calibrated
     list from 03 (q < 0.05).
  C. Validation at equal list size: take as many genes as HPA calls prognostic, but ranked by Cox p instead of min-p. Outcome as in 04:
     prognostic in HPA's validation cohort (p < 0.001 there) with the same direction. Same caveat as 04: the outcome is itself produced
     by the inflated procedure, so rates are comparable between lists but too generous in absolute terms.
Usage: python3 repro/05_cox.py [N_PERM] [GENES_PER_PERM]   (default 40 x 1000). Seed fixed.
Output: repro/ut/05_<ABBR>_gener.csv, repro/ut/05_null_<ABBR>.npy, repro/ut/05_sammanfattning.csv, repro/ut/05_validering.csv"""
import sys, time, importlib.util
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from logrank import prepare  # noqa: E402
from cox import cox_for_genes  # noqa: E402
spec = importlib.util.spec_from_file_location("rep", HERE / "02_reproducera.py"); rep = importlib.util.module_from_spec(spec); spec.loader.exec_module(rep)
UT = rep.UT; NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 40; NGEN = int(sys.argv[2]) if len(sys.argv) > 2 else 1000


def bh(p):
    n = len(p); o = np.argsort(p); q = p[o] * n / np.arange(1, n + 1); q = np.minimum.accumulate(q[::-1])[::-1]; ut = np.empty(n); ut[o] = np.minimum(q, 1); return ut


def en_kohort(abbr):
    t0 = time.time(); rng = np.random.default_rng(20260922 + sum(map(ord, abbr)))
    c, genes, m, tid, handelse, _ = rep.ladda(abbr); keep = np.nanmean(m, axis=1) > 1
    expr = np.log2(np.nan_to_num(m[keep].astype(np.float64)) + 1); prep = prepare(tid, handelse)
    r = cox_for_genes(expr, prep)
    d = pd.DataFrame(dict(Gene=genes[keep], cox_beta=r["beta"], cox_se=r["se"], cox_p=r["p_lrt"], cox_p_wald=r["p_wald"], cox_p_score=r["p_score"], cox_converged=r["converged"]))
    d["cox_q"] = bh(d["cox_p"].to_numpy()); d["cox_direction"] = np.where(d.cox_beta > 0, 1, -1)
    d = d.merge(pd.read_csv(UT / f"03_{abbr}_gener.csv")[["Gene", "p_ours", "direction", "p_hpa", "p_cal", "q"]].rename(columns={"q": "minp_q", "direction": "minp_direction"}), on="Gene", how="left")
    d.to_csv(UT / f"05_{abbr}_gener.csv", index=False)
    null = []
    for _ in range(NPERM):
        g = rng.choice(expr.shape[0], size=min(NGEN, expr.shape[0]), replace=False); perm = rng.permutation(expr.shape[1]); null.append(cox_for_genes(expr[g][:, perm], prep)["p_lrt"])
    null = np.sort(np.concatenate(null)); np.save(UT / f"05_null_{abbr}.npy", null); hpa = d.p_ours < 1e-3; cox = d.cox_q < 0.05; rob = hpa & (d.minp_q < 0.05)
    return dict(cohort=abbr, patients=len(tid), events=int(handelse.sum()), genes=len(d), not_converged=int((~d.cox_converged).sum()), null_draws=len(null),
                cox_null_median_p=round(float(np.median(null)), 3), cox_alpha_at_001=round(float((null < 1e-3).mean()), 5), cox_alpha_at_05=round(float((null < 0.05).mean()), 4),
                hpa_prognostic=int(hpa.sum()), minp_robust=int(rob.sum()), cox_q05=int(cox.sum()), cox_p001=int((d.cox_p < 1e-3).sum()),
                cox_q05_and_hpa=int((cox & hpa).sum()), cox_q05_not_hpa=int((cox & ~hpa).sum()), hpa_not_cox_q05=int((hpa & ~cox).sum()),
                robust_and_cox=int((rob & cox).sum()), same_direction_in_hpa_set=round(float((d.cox_direction[hpa] == d.minp_direction[hpa]).mean()), 4) if hpa.any() else np.nan,
                seconds=round(time.time() - t0, 1))


def validering():
    fac = pd.read_csv(rep.ROOT / "data" / "hpa" / "cancer_prognostic_data.tsv", sep="\t"); kol = list(fac.columns[3:]); fac["p"] = fac[kol].bfill(axis=1).iloc[:, 0]
    fac["dir"] = np.where(fac[kol].notna().idxmax(axis=1).str.contains("unfavorable"), 1, -1); fac = fac.dropna(subset=["p"]); rows = []
    for a, namn in rep.NAMN.items():
        v = fac[fac.Cancer == namn + " (validation)"]
        if v.empty: continue
        d = pd.read_csv(UT / f"05_{a}_gener.csv").dropna(subset=["p_ours"]).merge(v[["Gene", "p", "dir"]].rename(columns={"p": "p_val", "dir": "dir_val"}), on="Gene")
        ok = d.p_val < 1e-3; n = int((d.p_ours < 1e-3).sum())
        top_minp = d.nsmallest(n, "p_ours"); top_cox = d.nsmallest(n, "cox_p"); cq = d[d.cox_q < 0.05]
        val = lambda x, col: round(float((ok[x.index] & (x.dir_val == x[col])).mean()), 3) if len(x) else np.nan
        rows.append(dict(cohort=a, genes=len(d), list_size=n, minp_list_validated=val(top_minp, "minp_direction"), cox_list_validated=val(top_cox, "cox_direction"),
                         overlap=len(set(top_minp.Gene) & set(top_cox.Gene)), cox_q05=len(cq), cox_q05_validated=val(cq, "cox_direction"),
                         base_rate=round(float(np.mean([(ok & (d.dir_val == s)).mean() for s in (1, -1)])), 3)))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    from multiprocessing import Pool
    with Pool(7) as pool:
        ordning = sorted(rep.NAMN, key=lambda a: -{"BRCA": 9, "HNSC": 8, "LUSC": 8, "LUAD": 7, "KIRC": 7, "OV": 6, "STAD": 5, "LIHC": 4}.get(a, 1)); rader = []
        for r in pool.imap_unordered(en_kohort, ordning): rader.append(r); print(r, flush=True)
    pd.DataFrame(rader).sort_values("cohort").to_csv(UT / "05_sammanfattning.csv", index=False)
    v = validering(); v.to_csv(UT / "05_validering.csv", index=False); print(v.to_string(index=False))
