#!/usr/bin/env python3
"""Three figures for the two-page note (English, for a researcher). Reads only files written by 03, 05 and 07.
  fig1_calibration.png   actual vs nominal false-positive rate under permutation: best-cutoff log-rank vs univariate Cox, one line per cohort
  fig2_chance_share.png  genes below p < 0.001 in a permuted dataset (median and 5th to 95th percentile) vs genes HPA labels prognostic
  fig3_validation.png    validation rate at equal list size: ranked by best-cutoff p vs ranked by Cox p, with the base rate
Output: docs/figurer/*.png (300 dpi). No randomness."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]; UT = ROOT / "repro" / "ut"; FIG = ROOT / "docs" / "figurer"; FIG.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, GREY, INK, MUTED, SURFACE = "#2a78d6", "#eb6834", "#8a8984", "#0b0b0b", "#52514e", "#fcfcfb"
plt.rcParams.update({"font.family": "Helvetica Neue", "font.size": 9, "axes.edgecolor": "#c9c8c2", "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": MUTED, "text.color": INK, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE, "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.5, "xtick.major.size": 0, "ytick.major.size": 0})
s3 = pd.read_csv(UT / "03_sammanfattning.csv"); s5 = pd.read_csv(UT / "05_sammanfattning.csv")


def fig1():
    fig, ax = plt.subplots(figsize=(4.6, 4.2)); grid = np.logspace(-4, -1, 60)
    for k, (pre, col) in enumerate((("03", BLUE), ("05", ORANGE))):
        for a in s3.cohort:
            null = np.load(UT / f"{pre}_null_{a}.npy"); ax.plot(grid, np.searchsorted(null, grid, side="left") / len(null), color=col, lw=0.8, alpha=0.55, zorder=3 - k)
    ax.plot(grid, grid, color=INK, lw=1, ls=(0, (4, 3)), zorder=4); ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e-4, 1e-1); ax.set_ylim(3e-5, 1)
    ax.axvline(1e-3, color=MUTED, lw=0.6); ax.text(1.08e-3, 4e-5, "HPA threshold\np < 0.001", color=MUTED, fontsize=7.5, va="bottom")
    ax.text(1.9e-4, 6.5e-2, "Best cut-off log-rank\n(HPA procedure)", color=INK, fontsize=8, ha="left"); ax.plot([1.35e-4, 1.8e-4], [9.5e-2, 9.5e-2], color=BLUE, lw=2)
    ax.text(6e-3, 1.1e-3, "Univariate Cox\n(no cut-off)", color=INK, fontsize=8, ha="left"); ax.plot([4.4e-3, 5.7e-3], [1.6e-3, 1.6e-3], color=ORANGE, lw=2)
    ax.text(4.2e-2, 1.35e-2, "dashed: a valid test", color=MUTED, fontsize=7.5, ha="center")
    ax.set_xlabel("Nominal p-value threshold"); ax.set_ylabel("Actual false-positive rate under permutation")
    ax.set_title("One line per TCGA cohort (21). Survival permuted against expression.", fontsize=8, color=MUTED, loc="left")
    fig.tight_layout(); fig.savefig(FIG / "fig1_calibration.png", dpi=300); plt.close(fig)


def fig2():
    s7 = pd.read_csv(UT / "07_sammanfattning.csv").set_index("cohort"); d = s3.set_index("cohort").join(s7[["null_median", "null_p05", "null_p95", "p_global_minp"]])
    d["share_med"] = d.null_median / d.prognostic * 100; d["share_p95"] = d.null_p95 / d.prognostic * 100; d["share_p05"] = d.null_p05 / d.prognostic * 100
    d = d.sort_values("share_med"); fig, ax = plt.subplots(figsize=(5.6, 5.0)); y = np.arange(len(d)); cap = 160
    ax.barh(y, d.share_med.clip(upper=cap), height=0.62, color=BLUE, zorder=3)
    for yi, r in zip(y, d.itertuples()):
        ax.plot([min(r.share_p05, cap), min(r.share_p95, cap)], [yi, yi], color=INK, lw=1.1, zorder=4); ax.plot([min(r.share_p95, cap)], [yi], marker="|", color=INK, ms=7, mew=1.1, zorder=4)
        ax.text(min(r.share_p95, cap) + 3, yi, f"{r.null_median:.0f} ({r.null_p05:.0f} to {r.null_p95:.0f}) / {r.prognostic:,}   ·   {r.events} deaths", va="center", fontsize=6.6, color=MUTED)
    ax.set_yticks(y); ax.set_yticklabels(d.index, color=INK, fontsize=8); ax.set_xlim(0, 260); ax.set_xticks([0, 25, 50, 75, 100, 125, 150]); ax.axvline(100, color=MUTED, lw=0.6, zorder=2)
    ax.grid(axis="y", visible=False); ax.set_xlabel("Genes below p < 0.001 under permutation, as % of genes labelled prognostic (capped at 160)")
    ax.set_title("Bar: median over 200 permuted datasets (2,000 in eight cohorts). Line: 5th to 95th percentile.\nLabel: median (5th to 95th) / labelled", fontsize=8, color=MUTED, loc="left")
    fig.tight_layout(); fig.savefig(FIG / "fig2_chance_share.png", dpi=300); plt.close(fig)


def fig3():
    v = pd.read_csv(UT / "05_validering.csv").sort_values("minp_list_validated"); fig, ax = plt.subplots(figsize=(5.4, 3.9)); y = np.arange(len(v))
    for yi, r in zip(y, v.itertuples()): ax.plot([r.base_rate * 100, max(r.minp_list_validated, r.cox_list_validated) * 100], [yi, yi], color="#d6d5cf", lw=1.2, zorder=2)
    ax.scatter(v.base_rate * 100, y, s=26, color=GREY, marker="|", linewidths=1.6, zorder=3, label="Random gene (base rate)")
    ax.scatter(v.minp_list_validated * 100, y + 0.13, s=34, color=BLUE, edgecolor=SURFACE, linewidths=1, zorder=4, label="Ranked by best cut-off p (HPA)")
    ax.scatter(v.cox_list_validated * 100, y - 0.13, s=30, color=ORANGE, marker="D", edgecolor=SURFACE, linewidths=1, zorder=5, label="Same number of genes, ranked by Cox p")
    ax.set_yticks(y); ax.set_yticklabels([f"{c}  (n = {n:,})" for c, n in zip(v.cohort, v.list_size)], color=INK, fontsize=8); ax.grid(axis="y", visible=False)
    ax.set_xlabel("% of listed genes prognostic in HPA's validation cohort, same direction"); ax.set_xlim(0, 40)
    ax.legend(loc="lower right", frameon=False, fontsize=7.5, labelcolor=INK, handletextpad=0.3)
    fig.tight_layout(); fig.savefig(FIG / "fig3_validation.png", dpi=300); plt.close(fig)


if __name__ == "__main__":
    fig1(); fig2(); fig3(); print("wrote", *sorted(p.name for p in FIG.glob("*.png")))
