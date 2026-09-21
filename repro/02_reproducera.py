#!/usr/bin/env python3
"""Reproduces HPA's per-gene "best cut-off" log-rank p-values for one TCGA cohort and compares them with HPA's published values.

Usage: python3 repro/02_reproducera.py KIRC [LIHC ...]      (TCGA abbreviations; "alla" runs every cohort)
Input:  data/hpa/per_cancer/<slug>.npz (from 01_matriser.py), data/xena/Survival_SupplementalTable_S1_20171025_xena_sp (TCGA-CDR, overall
        survival), data/hpa/cancer_prognostic_data.tsv (HPA's answer key)
Output: repro/ut/02_<ABBR>_gener.csv (one row per gene: our p, HPA's p, cutoff, direction) and a line in repro/ut/02_sammanfattning.csv

Method as stated in Yuan et al. 2024: genes with mean expression > 1 within the cohort; cutoffs = every value between the 20th and 80th
percentile; keep the cutoff with the lowest log-rank p. Known unknowns, to be resolved against the authors' R code: which survival
table and endpoint they used, how duplicate samples per patient are handled, and whether a minimum group size applies."""
import re, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from logrank import prepare, minp_for_genes  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]; UT = ROOT / "repro" / "ut"
NAMN = {"BLCA": "Bladder Urothelial Carcinoma", "BRCA": "Breast Invasive Carcinoma", "CESC": "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma",
        "COAD": "Colon Adenocarcinoma", "GBM": "Glioblastoma Multiforme", "HNSC": "Head and Neck Squamous Cell Carcinoma", "KICH": "Kidney Chromophobe",
        "KIRC": "Kidney Renal Clear Cell Carcinoma", "KIRP": "Kidney Renal Papillary Cell Carcinoma", "LIHC": "Liver Hepatocellular Carcinoma",
        "LUAD": "Lung Adenocarcinoma", "LUSC": "Lung Squamous Cell Carcinoma", "OV": "Ovary Serous Cystadenocarcinoma", "PAAD": "Pancreatic Adenocarcinoma",
        "PRAD": "Prostate Adenocarcinoma", "READ": "Rectum Adenocarcinoma", "SKCM": "Skin Cuteneous Melanoma", "STAD": "Stomach Adenocarcinoma",
        "TGCT": "Testicular Germ Cell Tumor", "THCA": "Thyroid Carcinoma", "UCEC": "Uterine Corpus Endometrial Carcinoma"}
slug = lambda c: re.sub(r"[^a-z0-9]+", "-", c.lower()).strip("-")


def ladda(abbr):
    c = NAMN[abbr] + " (TCGA)"; z = np.load(ROOT / "data" / "hpa" / "per_cancer" / f"{slug(c)}.npz", allow_pickle=False)
    genes, samples, m = z["genes"], z["samples"], z["ptpm"]
    surv = pd.read_csv(ROOT / "data" / "xena" / "Survival_SupplementalTable_S1_20171025_xena_sp", sep="\t").drop_duplicates("_PATIENT").set_index("_PATIENT")
    pat = pd.Series([s[:12] for s in samples]); keep = ~pat.duplicated() & pat.isin(surv.index)
    s = surv.loc[pat[keep]]; ok = (s["OS"].notna() & s["OS.time"].notna() & (s["OS.time"] > 0)).to_numpy()
    idx = np.flatnonzero(keep.to_numpy())[ok]; s = s[ok]
    return c, genes, m[:, idx], s["OS.time"].to_numpy(float), s["OS"].to_numpy(int).astype(bool), len(samples)


def kor(abbr, facit):
    t0 = time.time(); c, genes, m, tid, handelse, n_prov = ladda(abbr)
    uttryckt = np.nanmean(m, axis=1) > 1; m2 = m[uttryckt]; g2 = genes[uttryckt]
    p, cut, nh, riktning = minp_for_genes(m2.astype(np.float64), prepare(tid, handelse))
    d = pd.DataFrame(dict(Gene=g2, p_ours=p, cutoff=cut, n_high=nh, direction=riktning))
    f = facit[facit["Cancer"] == c][["Gene", "p", "label"]].rename(columns={"p": "p_hpa"}); d = d.merge(f, on="Gene", how="outer", indicator=True)
    both = d[d["_merge"] == "both"]; rho = spearmanr(both["p_ours"], both["p_hpa"]).statistic
    hpa_dir = np.where(both["label"].str.contains("unfavorable"), 1, -1)
    res = dict(cohort=abbr, samples_hpa=n_prov, patients_used=len(tid), events=int(handelse.sum()), genes_ours=int(uttryckt.sum()), genes_hpa=len(f), genes_both=len(both),
               spearman_p=round(rho, 4), median_abs_log10_diff=round(float(np.median(np.abs(np.log10(both["p_ours"].clip(1e-300)) - np.log10(both["p_hpa"].clip(1e-300))))), 3),
               prognostic_ours=int((both["p_ours"] < 1e-3).sum()), prognostic_hpa=int((both["p_hpa"] < 1e-3).sum()),
               prognostic_both=int(((both["p_ours"] < 1e-3) & (both["p_hpa"] < 1e-3)).sum()), same_direction=round(float((both["direction"].to_numpy() == hpa_dir).mean()), 4),
               seconds=round(time.time() - t0, 1))
    d.drop(columns="_merge").to_csv(UT / f"02_{abbr}_gener.csv", index=False)
    return res


def main():
    fac = pd.read_csv(ROOT / "data" / "hpa" / "cancer_prognostic_data.tsv", sep="\t"); kol = list(fac.columns[3:])
    fac["p"] = fac[kol].bfill(axis=1).iloc[:, 0]; fac["label"] = fac[kol].notna().idxmax(axis=1); fac = fac.dropna(subset=["p"])
    vilka = list(NAMN) if sys.argv[1:] == ["alla"] else [a.upper() for a in sys.argv[1:]] or ["KIRC"]
    rader = []
    for a in vilka:
        r = kor(a, fac); rader.append(r); print(r, flush=True)
    ny = pd.DataFrame(rader); f = UT / "02_sammanfattning.csv"
    if f.exists(): ny = pd.concat([pd.read_csv(f).query("cohort not in @vilka"), ny])
    ny.to_csv(f, index=False)


if __name__ == "__main__":
    main()
