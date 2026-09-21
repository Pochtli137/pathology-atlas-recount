#!/usr/bin/env python3
"""Streams HPA's long-format TCGA expression file in chunks and writes one genes x samples matrix per cancer cohort.
Input:  data/hpa/rna_cancer_sample.tsv.gz  (columns Gene, Sample, Cancer, pTPM; HPA v25.1)
Output: data/hpa/per_cancer/<slug>.npz with arrays genes, samples, ptpm (float32)
Memory: 12 bytes per row (two int32 codes and a float32), about 2 GB for the whole file."""
import re, sys
from collections import defaultdict
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]; UT = ROOT / "data" / "hpa" / "per_cancer"; UT.mkdir(parents=True, exist_ok=True)
gene_ix, samp_ix = defaultdict(dict), defaultdict(dict); parts = defaultdict(list); n = 0
for ch in pd.read_csv(ROOT / "data" / "hpa" / "rna_cancer_sample.tsv.gz", sep="\t", chunksize=5_000_000, dtype={"Gene": "string", "Sample": "string", "Cancer": "string", "pTPM": "float32"}):
    for c, d in ch.groupby("Cancer", sort=False):
        gi, si = gene_ix[c], samp_ix[c]
        for g in d["Gene"].unique(): gi.setdefault(g, len(gi))
        for s_ in d["Sample"].unique(): si.setdefault(s_, len(si))
        parts[c].append((d["Gene"].map(gi).to_numpy(np.int32), d["Sample"].map(si).to_numpy(np.int32), d["pTPM"].to_numpy(np.float32)))
    n += len(ch); print(f"{n / 1e6:.0f} M rader", file=sys.stderr, flush=True)
for c, p in parts.items():
    g = np.concatenate([x[0] for x in p]); s_ = np.concatenate([x[1] for x in p]); v = np.concatenate([x[2] for x in p])
    m = np.full((len(gene_ix[c]), len(samp_ix[c])), np.nan, dtype=np.float32); m[g, s_] = v
    slug = re.sub(r"[^a-z0-9]+", "-", c.lower()).strip("-")
    np.savez_compressed(UT / f"{slug}.npz", genes=np.array(list(gene_ix[c])), samples=np.array(list(samp_ix[c])), ptpm=m, cancer=c)
    print(f"{c}: {m.shape[0]} gener x {m.shape[1]} prov, saknade värden {int(np.isnan(m).sum())}", flush=True)
