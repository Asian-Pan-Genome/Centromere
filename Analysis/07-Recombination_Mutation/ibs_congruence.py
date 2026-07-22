#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ibs_congruence.py — 着丝粒左右两侧 pairwise IBS 距离的一致性检验（不依赖树）。

对左、右两个 SNP FASTA（相同样本集），分别算每对样本的 IBS 距离
(= 两样本在双方都非缺失位点上的不一致比例)，然后:
  - Spearman / Pearson 相关 (dL vs dR over all C(n,2) pairs)
  - Mantel 置换检验 (置换 tip 标签，因 pairwise 元素不独立)
  - 散点图 dL vs dR

用法:
  python ibs_congruence.py --left L.fasta --right R.fasta --out-prefix compare/ibs \
         [--perm 9999] [--seed 1]
"""
import argparse, sys
import numpy as np
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument("--left", required=True)
ap.add_argument("--right", required=True)
ap.add_argument("--out-prefix", required=True)
ap.add_argument("--perm", type=int, default=9999)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--chrom", default="", help="标题用的染色体名(如 chr1)，留空则不带染色体名")
args = ap.parse_args()

MAP = {"A": 0, "C": 1, "G": 2, "T": 3}

def read_fa(fn):
    d, name = {}, None
    for line in open(fn):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:]; d[name] = []
        else:
            d[name].append(line)
    return {k: "".join(v) for k, v in d.items()}

def encode(seqdict, names):
    L = len(seqdict[names[0]])
    arr = np.full((len(names), L), -1, dtype=np.int8)
    for i, nm in enumerate(names):
        s = seqdict[nm]
        for j, c in enumerate(s):
            arr[i, j] = MAP.get(c, -1)   # N / other -> -1 (missing)
    return arr

def ibs_matrix(arr):
    n = arr.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        si = arr[i]
        for j in range(i + 1, n):
            sj = arr[j]
            mask = (si >= 0) & (sj >= 0)
            comp = mask.sum()
            diff = ((si != sj) & mask).sum()
            d = diff / comp if comp else np.nan
            D[i, j] = D[j, i] = d
    return D

faL, faR = read_fa(args.left), read_fa(args.right)
names = [n for n in faL if n in faR]
if set(faL) != set(faR):
    sys.stderr.write(f"[warn] tip sets differ; using {len(names)} shared\n")
print(f"samples compared: {len(names)}")

DL = ibs_matrix(encode(faL, names))
DR = ibs_matrix(encode(faR, names))

iu = np.triu_indices(len(names), k=1)
xl, xr = DL[iu], DR[iu]
ok = ~(np.isnan(xl) | np.isnan(xr))
xl, xr = xl[ok], xr[ok]
rho, p_sp = stats.spearmanr(xl, xr)
r, p_pe = stats.pearsonr(xl, xr)
print(f"pairs: {len(xl)}")
print(f"Spearman rho = {rho:.4f}  (naive p={p_sp:.2e}, NOT valid -> use Mantel)")
print(f"Pearson  r   = {r:.4f}")

# ---- Mantel test (Spearman), permute tip labels ----
rng = np.random.default_rng(args.seed)
n = len(names)
obs = rho
ge = 1
for _ in range(args.perm):
    perm = rng.permutation(n)
    DRp = DR[np.ix_(perm, perm)]
    xrp = DRp[iu][ok]
    rp, _ = stats.spearmanr(xl, xrp)
    if rp >= obs:
        ge += 1
p_mantel = ge / (args.perm + 1)
print(f"Mantel (Spearman, {args.perm} perms): p = {p_mantel:.4g}")

# ---- scatter ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, axx = plt.subplots(figsize=(6, 6))
axx.scatter(xl, xr, s=4, alpha=0.25, edgecolors="none", zorder=2)
axx.plot([0, 1], [0, 1], "r--", lw=1, zorder=3)            # y = x reference
# fixed 0-1 axes, ticks every 0.1
axx.set_xlim(0, 1); axx.set_ylim(0, 1)
ticks = np.arange(0, 1.0001, 0.1)
axx.set_xticks(ticks); axx.set_yticks(ticks)
axx.set_aspect("equal", adjustable="box")
axx.set_xlabel("p-arm pairwise IBS distance")
axx.set_ylabel("q-arm pairwise IBS distance")
_pre = (args.chrom + " ") if args.chrom else ""
axx.set_title(f"{_pre}peri-centromeric IBS congruence  ($\\it{{n}}$={n}, {len(xl)} pairs)")
# y=x label + stats in the upper-left; p in italic (mathtext)
p_disp = "< 0.0001" if p_mantel <= 1.0/(args.perm+1) else f"= {p_mantel:.4g}"
axx.text(0.04, 0.94, "y = x", color="red", transform=axx.transAxes,
         va="top", ha="left", fontsize=11)
axx.text(0.04, 0.88, f"Spearman $\\rho$ = {rho:.3f}\nMantel $\\it{{p}}$ {p_disp}",
         transform=axx.transAxes, va="top", ha="left", fontsize=11)
fig.tight_layout()
fig.savefig(args.out_prefix + "_scatter.pdf")
fig.savefig(args.out_prefix + "_scatter.png", dpi=150)
print("wrote", args.out_prefix + "_scatter.pdf / .png")

with open(args.out_prefix + "_summary.txt", "w") as o:
    o.write(f"samples\t{n}\npairs\t{len(xl)}\n")
    o.write(f"spearman_rho\t{rho:.4f}\npearson_r\t{r:.4f}\nmantel_p\t{p_mantel:.4g}\nperm\t{args.perm}\n")
print("wrote", args.out_prefix + "_summary.txt")
