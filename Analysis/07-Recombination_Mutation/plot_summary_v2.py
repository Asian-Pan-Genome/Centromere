#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plot_summary.py — 24 条染色体着丝粒左右两侧 phylogeny 一致性总表可视化。

读 compare_all_summary.tsv（complete 列 = 去掉 filterflag!=0 不完整 node 的树），
自然染色体顺序 chr1..chr22, chrX, chrY。三个指标：
  - IBS Spearman ρ          (IBSrho_complete)   红
  - Quartet similarity @70  (Q70_complete)      蓝
  - MCI similarity          (MCI70 或 MCI0_complete) 绿

输出两张图：
  1) compare_all_summary_grouped.{pdf,png}  分组柱状图（每条染色体 3 根柱，MCI@70）
  2) compare_all_summary_panels.{pdf,png}   方案1 三分面柱状图（MCI@0 / Q70 / IBS ρ）
"""
import argparse
import csv
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ===== 设置全局字体为 Arial，PDF 可编辑 =====
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42   # TrueType，Illustrator 可编辑
plt.rcParams['ps.fonttype'] = 42

RED, BLUE, GREEN = "#d62728", "#1f77b4", "#2ca02c"


def chr_key(c):
    s = c[3:]
    return {"X": 100, "Y": 101}[s] if s in ("X", "Y") else int(s)


def load(path, mci70_path=None):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append(r)
    rows.sort(key=lambda r: chr_key(r["chrom"]))
    chroms = [r["chrom"] for r in rows]
    ibs = [float(r["IBSrho_complete"]) for r in rows]
    q70 = [float(r["Q70_complete"]) for r in rows]
    mci0 = [float(r["MCI0_complete"]) for r in rows]
    # MCI@70（去噪、与 Quartet@70 同口径）：从各 chrN/compare 的 treecong_complete 取
    mci70 = load_mci70(chroms) if mci70_path is None else None
    return chroms, np.array(ibs), np.array(q70), np.array(mci0), mci70


def load_mci70(chroms):
    """从 chrN/compare/chrN_treecong_complete.txt (chrY 特殊) 取 thresh==70 行的 MCI_normSim。"""
    base = os.path.join(os.path.dirname(__file__), "..")
    out = []
    for c in chroms:
        if c == "chrY":
            p = os.path.join(base, "chrY", "compare", "tree_congruence_LvR.txt")
        else:
            p = os.path.join(base, c, "compare", f"{c}_treecong_complete.txt")
        val = np.nan
        try:
            with open(p) as f:
                for line in f:
                    t = line.split()
                    if t and t[0] == "70":
                        val = float(t[3])   # MCI_normSim
                        break
        except FileNotFoundError:
            pass
        out.append(val)
    return np.array(out)


def grouped(chroms, ibs, q70, mci70, prefix):
    x = np.arange(len(chroms))
    w = 0.27
    fig, ax = plt.subplots(figsize=(15, 4.6))
    ax.bar(x - w, ibs, w, color=RED, label="IBS Spearman ρ")
    ax.bar(x, q70, w, color=BLUE, label="Quartet similarity")
    ax.bar(x + w, mci70, w, color=GREEN, label="MCI similarity")
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("chr", "") for c in chroms], fontsize=9, fontfamily='Arial')
    ax.set_xlabel("Chromosome", fontfamily='Arial')
    ax.set_ylabel("Left–right congruence", fontfamily='Arial')
    ax.set_ylim(0, 1.0)
    ax.set_yticks(np.arange(0, 1.01, 0.1))
    ax.grid(axis="y", ls=":", lw=0.6, alpha=0.6)
    ax.set_axisbelow(True)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3,
              frameon=False, fontsize=10, prop={'family': 'Arial'})
    
    # 设置刻度标签字体
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontfamily('Arial')
    
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{prefix}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[done] {prefix}.pdf / .png")


def panels(chroms, mci0, q70, ibs, prefix, stat="median"):
    x = np.arange(len(chroms))
    P1, P2, P3 = "#956BBF", "#F2BF27", "#5D87B3"
    specs = [("Mutual Clustering Information similarity", mci0, P1),
             ("Quartet similarity", q70, P2),
             ("IBS (Spearman ρ)", ibs, P3)]
    fig, axes = plt.subplots(3, 1, figsize=(14, 9), sharex=True)
    for ax, (title, y, col) in zip(axes, specs):
        ax.bar(x, y, 0.7, color=col)
        ax.set_ylabel(title, fontsize=10, fontfamily='Arial')
        ax.set_ylim(0, 1.0)
        ax.set_yticks(np.arange(0, 1.01, 0.2))
        if stat == "mean":
            ref = float(np.nanmean(y))
        else:
            ref = float(np.nanmedian(y))
        ax.axhline(ref, color="grey", ls="--", lw=1.0, zorder=1)
        ax.text(-0.3, ref + 0.025, f"{stat} {ref:.2f}", color="grey",
                va="bottom", ha="left", fontsize=8, fontfamily='Arial')
        for xi, yi in zip(x, y):
            if np.isnan(yi):
                continue
            ax.text(xi, yi + 0.015, f"{yi:.2f}", ha="center", va="bottom",
                    fontsize=7, fontfamily='Arial')
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        # 设置 y 轴刻度标签字体
        for label in ax.get_yticklabels():
            label.set_fontfamily('Arial')
    
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels([c.replace("chr", "") for c in chroms], fontsize=9, fontfamily='Arial')
    axes[-1].set_xlabel("Chromosome", fontfamily='Arial')
    fig.suptitle("centromere p-arm vs. q-arm phylogeny congruence", fontsize=13, fontfamily='Arial')
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    for ext in ("pdf", "png"):
        fig.savefig(f"{prefix}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[done] {prefix}.pdf / .png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", default="compare_all_summary.tsv")
    ap.add_argument("--out-grouped", default="compare_all_summary_grouped")
    ap.add_argument("--out-panels", default="compare_all_summary_panels")
    a = ap.parse_args()
    chroms, ibs, q70, mci0, mci70 = load(a.tsv)
    grouped(chroms, ibs, q70, mci70, a.out_grouped)
    panels(chroms, mci0, q70, ibs, a.out_panels, stat="median")
    panels(chroms, mci0, q70, ibs, a.out_panels + "_mean", stat="mean")


if __name__ == "__main__":
    main()