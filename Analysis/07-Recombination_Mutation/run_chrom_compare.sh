#!/usr/bin/env bash
# run_chrom_compare.sh <chr>
# 单染色体 left vs right 树比较 + tanglegram 可视化（参照 chrY）。
# 前置：已跑 prep_compare_inputs.py（生成 map / drop 清单 / full|complete FASTA）。
# 比较两版：full（去 GRCh38）与 complete（再去 filterflag!=0 的不完整样本）。
# 树保持 unrooted。可视化(tanglegram hicat+hormon)用 complete 树。
# 从 SNP_phylogeny/ 运行（../cenhaps 命中项目根 cenhaps）。
set -euo pipefail
chr=${1:?usage: bash run_chrom_compare.sh <chr>}
d=$chr
RS=${RS:-Rscript}
PY=${PY:-python}
TANG=../add_ape_phylogeny/scripts/tanglegram.R
mkdir -p "$d/compare" "$d/figure"

# 1) 改名 + 删除 -> full / complete 树（unrooted）
for side in left right; do
  tf="$d/${chr}_${side}.fasta.treefile"
  [ -f "$tf" ] || { echo "[skip] no $tf"; continue; }
  "$RS" scripts/relabel_root_tree.R --tree "$tf" --map "$d/$chr.map.tsv" \
      --drop-file "$d/$chr.drop_full.txt"     --out "$d/${chr}_${side}.full.nwk"
  "$RS" scripts/relabel_root_tree.R --tree "$tf" --map "$d/$chr.map.tsv" \
      --drop-file "$d/$chr.drop_complete.txt" --out "$d/${chr}_${side}.complete.nwk"
done

# 2) 比较：full + complete（MCI/Quartet + RF/Jaccard + IBS）
for ver in full complete; do
  echo "===== $chr $ver ====="
  "$PY" scripts/tree_congruence.py \
      --tree-a "$d/${chr}_left.$ver.nwk" --tree-b "$d/${chr}_right.$ver.nwk" \
      --label-a LEFT --label-b RIGHT --thresholds 0,50,70,95 \
      --quartets 200000 --perm 300 --out "$d/compare/${chr}_treecong_$ver.txt"
  "$RS" ../add_ape_phylogeny/scripts/compare_trees.R \
      --tree-a "$d/${chr}_left.$ver.nwk" --tree-b "$d/${chr}_right.$ver.nwk" \
      --label-a LEFT --label-b RIGHT --out-prefix "$d/compare/${chr}_cmp_$ver" \
      > "$d/compare/${chr}_cmp_$ver.log" 2>&1 || true
  "$PY" scripts/ibs_congruence.py \
      --left "$d/${chr}_left.$ver.fasta" --right "$d/${chr}_right.$ver.fasta" \
      --out-prefix "$d/compare/${chr}_ibs_$ver" --perm 9999 --seed 1
done

# 3) tanglegram（complete 树）hicat + hormon
for src in hicat hormon; do
  "$RS" "$TANG" "$chr" \
      --left-tree  "$d/${chr}_left.complete.nwk" \
      --right-tree "$d/${chr}_right.complete.nwk" \
      --hor-source "$src" --outdir "$d/figure/" || echo "[warn] tanglegram $src failed"
done
echo "[done] $chr -> $d/compare/ , $d/figure/"
