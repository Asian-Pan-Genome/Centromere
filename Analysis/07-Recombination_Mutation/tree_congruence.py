#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
tree_congruence.py — 两棵树拓扑一致性（不依赖 TreeDist/Quartet R 包，纯 Python）。

指标:
  1. 归一化 Mutual Clustering Information 相似度 (Smith 2020 思路):
     每条内枝 = 一个 split(二分)；两树 split 间用互信息打分；最大权一对一匹配；
     normSim = 2*MCI / (E_A + E_B)，1 = 完全一致。给「几乎匹配」的 split 部分分。
  2. Quartet 相似度 (Monte-Carlo): 随机抽 K 个 4-样本子集，看两树诱导的解析
     (3 种二分之一 / 未解析) 是否一致。报告 resolved 一致率与 strict 一致率。

支持 bootstrap 折叠: 对每个阈值，只保留 support >= T 的 split (= 折叠低支持枝)。
显著性: 置换一侧 tip 标签得 null 分布 (pairwise/quartet 元素不独立, 故用置换)。

用法:
  python tree_congruence.py --tree-a L.treefile --tree-b R.treefile \
      --thresholds 0,50,70,95 --quartets 200000 --perm 500 --out compare/tree_congruence.txt
"""
import argparse, sys
import numpy as np
from scipy.optimize import linear_sum_assignment
from Bio import Phylo

ap = argparse.ArgumentParser()
ap.add_argument("--tree-a", required=True)
ap.add_argument("--tree-b", required=True)
ap.add_argument("--label-a", default="A")
ap.add_argument("--label-b", default="B")
ap.add_argument("--thresholds", default="0,50,70,95")
ap.add_argument("--quartets", type=int, default=200000)
ap.add_argument("--perm", type=int, default=500)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--out", default=None)
args = ap.parse_args()

def load_splits(tr, name2idx, n):
    """return list of (bool_membership[n], support, size); only 2<=size<=n-2.
    `tr` is a Bio.Phylo tree object (already pruned to the common tip set)."""
    terms = tr.get_terminals()
    out = []
    for cl in tr.get_nonterminals():
        leaves = [t.name for t in cl.get_terminals()]
        size = len(leaves)
        if size < 2 or size > n - 2:
            continue
        sup = cl.confidence
        if sup is None:
            try: sup = float(cl.name)
            except (TypeError, ValueError): sup = 0.0
        memb = np.zeros(n, dtype=bool)
        for lf in leaves:
            memb[name2idx[lf]] = True
        out.append((memb, float(sup), size))
    return out, [t.name for t in terms]

# tip universe: prune both trees to the COMMON tip set (left/right may differ
# slightly after per-side gappy-sample dropping)
trA = Phylo.read(args.tree_a, "newick"); tipsA = set(t.name for t in trA.get_terminals())
trB = Phylo.read(args.tree_b, "newick"); tipsB = set(t.name for t in trB.get_terminals())
common = tipsA & tipsB
if len(common) < 4:
    sys.exit(f"only {len(common)} common tips (A={len(tipsA)} B={len(tipsB)})")
if tipsA != tipsB:
    print(f"[info] tip sets differ (A={len(tipsA)} B={len(tipsB)}); "
          f"pruning both to {len(common)} common tips")
    for tr, ts in ((trA, tipsA), (trB, tipsB)):
        for t in [x for x in tr.get_terminals() if x.name not in common]:
            tr.prune(t)
tips = sorted(common); n = len(tips); name2idx = {t: i for i, t in enumerate(tips)}
splitsA, _ = load_splits(trA, name2idx, n)
splitsB, _ = load_splits(trB, name2idx, n)
print(f"tips: {n} | splits A={len(splitsA)} B={len(splitsB)}")

def entropy_of_size(sz, n):
    p = sz / n
    return -(p*np.log2(p) + (1-p)*np.log2(1-p)) if 0 < p < 1 else 0.0

def mi_matrix(MA, MB, n):
    """MA [SA,n] bool, MB [SB,n] bool -> mutual info matrix [SA,SB] in bits."""
    A = MA.astype(np.float64); B = MB.astype(np.float64)
    a = A @ B.T                       # |A_i ∩ B_j|
    rs = A.sum(1)[:, None]; cs = B.sum(1)[None, :]
    b = rs - a; c = cs - a; d = n - a - b - c
    out = np.zeros_like(a)
    for x, r, col in [(a, rs, cs), (b, rs, n-cs), (c, n-rs, cs), (d, n-rs, n-cs)]:
        with np.errstate(divide="ignore", invalid="ignore"):
            term = (x/n) * np.log2((x/n) / ((r/n)*(col/n)))
        out += np.where(x > 0, term, 0.0)
    return out

def mci_norm(MA, MB, eA, eB, n):
    if len(MA) == 0 or len(MB) == 0: return 0.0
    M = mi_matrix(MA, MB, n)
    ri, ci = linear_sum_assignment(-M)
    mci = M[ri, ci].sum()
    return 2*mci / (eA + eB) if (eA+eB) > 0 else 0.0

def quartet_agreement(MA, MB, n, K, rng):
    """MC: sample K quartets; per tree resolution from split membership."""
    if len(MA) == 0 or len(MB) == 0:
        return dict(agree=0, disagree=0, neutral=K, res_sim=float("nan"), strict=0.0)
    q = np.empty((K, 4), dtype=np.int64)
    for c in range(4):
        q[:, c] = rng.integers(0, n, size=K)
    # ensure 4 distinct
    bad = np.ones(K, dtype=bool)
    while bad.any():
        idx = np.where(bad)[0]
        for c in range(4):
            q[idx, c] = rng.integers(0, n, size=idx.size)
        bad = (q[:,0]==q[:,1])|(q[:,0]==q[:,2])|(q[:,0]==q[:,3])|(q[:,1]==q[:,2])|(q[:,1]==q[:,3])|(q[:,2]==q[:,3])
    def resolve(M):
        # returns array K with 0/1/2 pairing or -1 unresolved
        S = M  # [S,n] bool
        a,b,c,d = q[:,0],q[:,1],q[:,2],q[:,3]
        Sa,Sb,Sc,Sd = S[:,a],S[:,b],S[:,c],S[:,d]   # [S,K]
        # pairing0 ab|cd, pairing1 ac|bd, pairing2 ad|bc
        p0 = ((Sa==Sb)&(Sc==Sd)&(Sa!=Sc)).any(0)
        p1 = ((Sa==Sc)&(Sb==Sd)&(Sa!=Sb)).any(0)
        p2 = ((Sa==Sd)&(Sb==Sc)&(Sa!=Sb)).any(0)
        res = np.full(K, -1, dtype=np.int8)
        res[p0]=0; res[p1]=1; res[p2]=2
        return res
    rA, rB = resolve(MA), resolve(MB)
    both = (rA>=0)&(rB>=0)
    agree = int(((rA==rB)&both).sum())
    disagree = int(((rA!=rB)&both).sum())
    neutral = int((~both).sum())
    res_sim = agree/(agree+disagree) if (agree+disagree)>0 else float("nan")
    return dict(agree=agree, disagree=disagree, neutral=neutral, res_sim=res_sim, strict=agree/K)

rng = np.random.default_rng(args.seed)
thresholds = [float(x) for x in args.thresholds.split(",")]
lines = [f"# tree congruence  {args.label_a} vs {args.label_b}  (n={n} tips)",
         f"# A splits={len(splitsA)}  B splits={len(splitsB)}",
         f"{'thresh':>7} {'nA':>4} {'nB':>4} {'MCI_normSim':>11} {'MCI_p':>7} "
         f"{'Quartet_resSim':>14} {'Q_strict':>9} {'Q_neutral%':>10} {'Quartet_p':>9}"]

for T in thresholds:
    MA = np.array([m for m, s, sz in splitsA if s >= T]) if any(s>=T for m,s,sz in splitsA) else np.zeros((0,n),bool)
    MB = np.array([m for m, s, sz in splitsB if s >= T]) if any(s>=T for m,s,sz in splitsB) else np.zeros((0,n),bool)
    szA = [sz for m,s,sz in splitsA if s>=T]; szB=[sz for m,s,sz in splitsB if s>=T]
    eA = sum(entropy_of_size(sz,n) for sz in szA); eB = sum(entropy_of_size(sz,n) for sz in szB)
    obs_mci = mci_norm(MA, MB, eA, eB, n)
    q = quartet_agreement(MA, MB, n, args.quartets, rng)
    # permutation null (relabel B columns)
    ge_m = ge_q = 1
    nullm = []; nullq = []
    for _ in range(args.perm):
        perm = rng.permutation(n)
        MBp = MB[:, perm] if len(MB) else MB
        m = mci_norm(MA, MBp, eA, eB, n); nullm.append(m)
        if m >= obs_mci: ge_m += 1
    # quartet perm (fewer, cheaper K)
    Kq = min(args.quartets, 40000)
    qperm = max(1, args.perm//5)
    for _ in range(qperm):
        perm = rng.permutation(n)
        MBp = MB[:, perm] if len(MB) else MB
        qq = quartet_agreement(MA, MBp, n, Kq, rng); nullq.append(qq["res_sim"])
        if not np.isnan(qq["res_sim"]) and qq["res_sim"] >= q["res_sim"]: ge_q += 1
    p_m = ge_m/(args.perm+1); p_q = ge_q/(qperm+1)
    lines.append(f"{T:>7.0f} {len(MA):>4} {len(MB):>4} {obs_mci:>11.4f} {p_m:>7.4g} "
                 f"{q['res_sim']:>14.4f} {q['strict']:>9.4f} {100*q['neutral']/args.quartets:>10.1f} {p_q:>9.4g}")
    print(lines[-1])

txt = "\n".join(lines) + "\n"
if args.out:
    open(args.out, "w").write(txt)
    print("wrote", args.out)
