#!/usr/bin/env Rscript
# relabel_root_tree.R
# 读 IQ-TREE treefile -> 按映射表(old\tnew)重命名 tip -> 可选 root 在指定 tip
# -> 写出 newick。用于把 chrY SNP 树的 tip 改成 sample_hap / 跨树比较用的统一键。
#
# 用法:
#   Rscript relabel_root_tree.R --tree IN.treefile --map map.tsv --out OUT.treefile [--root LABEL]
#   --map  : TSV，含表头 old\tnew；--root 用「映射后」的 new 标签
suppressMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, d = NA) { i <- match(k, args); if (is.na(i)) d else args[i + 1] }
tree_p <- getarg("--tree"); map_p <- getarg("--map")
out_p  <- getarg("--out");  root_lab <- getarg("--root")
drop_s <- getarg("--drop")        # comma-separated ORIGINAL tip labels to remove
drop_f <- getarg("--drop-file")   # file of ORIGINAL tip labels (one per line)
stopifnot(!is.na(tree_p), !is.na(map_p), !is.na(out_p))

tr  <- read.tree(tree_p)
drop <- character(0)
if (!is.na(drop_s)) drop <- strsplit(drop_s, ",")[[1]]
if (!is.na(drop_f) && file.exists(drop_f))
  drop <- c(drop, readLines(drop_f))
drop <- intersect(unique(trimws(drop)), tr$tip.label)
if (length(drop) > 0) {
  tr <- drop.tip(tr, drop)
  message("dropped ", length(drop), " tips")
}
map <- read.table(map_p, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                  quote = "", comment.char = "")
lut <- setNames(map$new, map$old)

miss <- setdiff(tr$tip.label, names(lut))
if (length(miss) > 0)
  stop("tips not in map: ", paste(head(miss, 5), collapse = ", "),
       if (length(miss) > 5) " ..." else "")
tr$tip.label <- unname(lut[tr$tip.label])
if (any(duplicated(tr$tip.label)))
  stop("duplicate tip labels after relabel: ",
       paste(tr$tip.label[duplicated(tr$tip.label)], collapse = ", "))

if (!is.na(root_lab)) {
  if (!root_lab %in% tr$tip.label)
    stop("root label '", root_lab, "' not found among tips")
  tr <- root(tr, outgroup = root_lab, resolve.root = TRUE)
  message("rooted on: ", root_lab)
}

write.tree(tr, out_p)
message("wrote ", out_p, "  (", length(tr$tip.label), " tips)")
