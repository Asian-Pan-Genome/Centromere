#!/usr/bin/env Rscript
# compare_trees.R
#
# Pairwise comparison of two phylogenetic trees:
#   - Robinson-Foulds (RF) distance (raw + normalized)
#   - Common / unique tip counts
#   - Simple cophyloplot tanglegram for visual sanity check
#
# Designed for "is window X giving the same signal as window Y" sorts
# of question (e.g. chr1 LEFT W1 vs LEFT W2, or W1 vs RIGHT).
#
# Two input modes:
#   (a) raw IQ-TREE outputs       --> pass --chrom + (optional)
#                                      --{a,b}-dropped / --{a,b}-mismatch
#                                      to apply the same tip-cleaning /
#                                      renaming / outgroup-rooting that
#                                      scripts/tanglegram.R does.
#   (b) already-cleaned newicks   --> omit --chrom; tip names are used
#                                      as-is.
#
# Tip set is restricted to the intersection of the two trees before
# computing RF (otherwise the metric is undefined). RF is computed on
# the UNROOTED tree, so different root choices do not affect the result.
#
# Usage:
#   Rscript scripts/compare_trees.R \
#       --tree-a <treefile_or_nwk> --tree-b <treefile_or_nwk> \
#       [--a-dropped <f>] [--b-dropped <f>] \
#       [--a-mismatch <f>] [--b-mismatch <f>] \
#       [--chrom chr1] \
#       [--label-a "W1"] [--label-b "W2"] \
#       --out-prefix <path/to/W1_vs_W2>
#
# Outputs:
#   <out_prefix>.report.txt        text report with metrics + tip lists
#   <out_prefix>.tanglegram.pdf    cophyloplot, common tips connected

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
get_named <- function(args, key, default = NA_character_) {
  i <- match(key, args)
  if (is.na(i) || i + 1 > length(args)) default else args[i + 1]
}

tree_a_p   <- get_named(args, "--tree-a")
tree_b_p   <- get_named(args, "--tree-b")
a_drop_p   <- get_named(args, "--a-dropped")
b_drop_p   <- get_named(args, "--b-dropped")
a_mis_p    <- get_named(args, "--a-mismatch")
b_mis_p    <- get_named(args, "--b-mismatch")
chrom      <- get_named(args, "--chrom")
label_a    <- get_named(args, "--label-a", "A")
label_b    <- get_named(args, "--label-b", "B")
out_pref   <- get_named(args, "--out-prefix", "compare_trees")
boot_str   <- get_named(args, "--boot-thresholds", "30,50,70")
jacc_str   <- get_named(args, "--jacc-sizes",      "5,10,30,50")
ari_str    <- get_named(args, "--ari-k",           "2,3,4")
boot_thresholds <- as.numeric(strsplit(boot_str, ",", fixed = TRUE)[[1]])
jacc_sizes      <- as.integer(strsplit(jacc_str, ",", fixed = TRUE)[[1]])
ari_k_values    <- as.integer(strsplit(ari_str,  ",", fixed = TRUE)[[1]])

if (is.na(tree_a_p) || is.na(tree_b_p))
  stop("--tree-a and --tree-b are required")
dir.create(dirname(out_pref), recursive = TRUE, showWarnings = FALSE)

# ---- helpers (subset of tanglegram.R; only run when --chrom is given) --

normalize_iqtree <- function(s) gsub("[:#()+]", "_", s)
read_drop_list <- function(path) {
  if (is.na(path) || !file.exists(path)) return(character(0))
  lines <- readLines(path)
  if (length(lines) <= 1) return(character(0))
  vapply(lines[-1], normalize_iqtree, character(1), USE.NAMES = FALSE)
}
read_mismatch_list <- function(path) {
  if (is.na(path) || !file.exists(path)) return(character(0))
  df <- read.table(path, header = TRUE, sep = "\t", quote = "",
                   stringsAsFactors = FALSE, comment.char = "")
  if (!"status" %in% names(df) || !"aln_name" %in% names(df))
    return(character(0))
  vapply(df$aln_name[df$status == "mismatch"], normalize_iqtree,
         character(1), USE.NAMES = FALSE)
}
rename_leaf <- function(label, chrom) {
  if (grepl("at_hsa", label))        return(paste0("bonobo_", chrom))
  if (grepl("hap[0-9]+_hsa", label)) return(paste0("chimp_",  chrom))
  if (grepl("^chr[0-9XYM]+_[0-9]", label) && !grepl("_hsa", label))
    return(paste0("CHM13_", chrom))
  parts <- strsplit(label, "_", fixed = TRUE)[[1]]
  hap_set <- c("Mat", "Pat", "hap1", "hap2")
  hap_idx <- which(parts %in% hap_set)[1]
  if (!is.na(hap_idx) && hap_idx > 1) {
    return(paste0(paste(parts[seq_len(hap_idx - 1)], collapse = "_"),
                  "_", parts[hap_idx]))
  }
  sub <- strsplit(parts[1], "-", fixed = TRUE)[[1]]
  if (length(sub) >= 3) {
    sample <- sub[1]; tok <- sub[2]
    if (tok %in% c("hap1", "hap2")) return(paste0(sample, "_", tok))
    if (tok == "1") return(paste0(sample, "_Pat"))
    if (tok == "2") return(paste0(sample, "_Mat"))
  }
  NA_character_
}

clean_tree <- function(tree_p, drop_p, mis_p, chrom_arg, side) {
  tr <- read.tree(tree_p)
  if (is.na(chrom_arg)) {
    message("loaded ", side, ": ", tree_p, "  (", length(tr$tip.label),
            " tips, no rename)")
    return(tr)
  }
  to_drop <- unique(c(read_drop_list(drop_p), read_mismatch_list(mis_p)))
  to_drop_in <- intersect(to_drop, tr$tip.label)
  if (length(to_drop_in) > 0) tr <- drop.tip(tr, to_drop_in)
  new <- vapply(tr$tip.label, rename_leaf, character(1),
                chrom = chrom_arg, USE.NAMES = FALSE)
  na_idx <- which(is.na(new))
  if (length(na_idx) > 0) new[na_idx] <- tr$tip.label[na_idx]
  if (any(duplicated(new))) {
    seen <- list()
    for (j in seq_along(new)) {
      lab <- new[j]
      if (sum(new == lab) > 1) {
        seen[[lab]] <- if (is.null(seen[[lab]])) 1 else seen[[lab]] + 1
        if (seen[[lab]] > 1) new[j] <- paste0(lab, "_", seen[[lab]])
      }
    }
  }
  tr$tip.label <- new
  og <- intersect(c(paste0("bonobo_", chrom_arg),
                    paste0("chimp_",  chrom_arg)), tr$tip.label)
  if (length(og) > 0) tr <- root(tr, outgroup = og, resolve.root = TRUE)
  message("loaded ", side, ": ", tree_p, "  (", length(tr$tip.label),
          " tips, dropped ", length(to_drop_in), ", renamed, ",
          ifelse(length(og) > 0, "rooted on outgroup", "unrooted"), ")")
  tr
}

# ---- load + clean -------------------------------------------------------
ta <- clean_tree(tree_a_p, a_drop_p, a_mis_p, chrom, label_a)
tb <- clean_tree(tree_b_p, b_drop_p, b_mis_p, chrom, label_b)

n_a_orig <- length(ta$tip.label)
n_b_orig <- length(tb$tip.label)

only_a <- setdiff(ta$tip.label, tb$tip.label)
only_b <- setdiff(tb$tip.label, ta$tip.label)
if (length(only_a) > 0) ta <- drop.tip(ta, only_a)
if (length(only_b) > 0) tb <- drop.tip(tb, only_b)
n_common <- length(ta$tip.label)
stopifnot(n_common == length(tb$tip.label))
if (n_common < 4)
  stop("only ", n_common, " common tips -- need >= 4 for RF")

# ---- RF (raw + bootstrap-filtered) + Major-clade Jaccard ----------------
# Robust manual bipartition extraction via prop.part() (ape::dist.topo
# misbehaves on polytomy-heavy or residually-rooted trees). prop.part is
# rooting-invariant after canonicalization, so we feed the cleaned-and-
# rooted trees directly (the artificial root introduces no extra split
# once we dedupe canonical strings + skip trivial sides).
#
# get_bipartitions returns the SET of canonical bipartition strings whose
# corresponding internal node has bootstrap support >= `threshold`. When
# threshold = -Inf, every split is kept. When threshold = 50, only splits
# with bootstrap >= 50 are kept (the "high-support consensus" view).
get_bipartitions <- function(tr, threshold = -Inf, min_size = 2) {
  pp <- prop.part(tr)
  tips <- tr$tip.label
  n_tips <- length(tips)
  sup <- if (is.null(tr$node.label)) rep(NA_real_, tr$Nnode)
         else suppressWarnings(as.numeric(tr$node.label))
  out <- character(0)
  for (i in seq_along(pp)) {
    side1 <- sort(tips[pp[[i]]])
    side2 <- sort(setdiff(tips, side1))
    if (min(length(side1), length(side2)) < min_size) next
    sup_val <- if (i <= length(sup)) sup[i] else NA_real_
    # NA support = unknown -> keep (don't filter out unknowns).
    if (!is.na(sup_val) && sup_val < threshold) next
    canon <- if (paste(side1, collapse = ",") < paste(side2, collapse = ","))
               paste(side1, collapse = ",")
             else
               paste(side2, collapse = ",")
    out <- c(out, canon)
  }
  unique(out)
}

# Pre-compute raw bipartition sets (used for the headline RF + as the
# basis for filtered versions).
ba_raw <- get_bipartitions(ta, threshold = -Inf)
bb_raw <- get_bipartitions(tb, threshold = -Inf)
n_a_bip   <- length(ba_raw)
n_b_bip   <- length(bb_raw)
shared    <- length(intersect(ba_raw, bb_raw))
rf        <- (n_a_bip - shared) + (n_b_bip - shared)
max_rf_bin  <- 2 * (n_common - 3)
max_rf_real <- n_a_bip + n_b_bip
rf_norm     <- if (max_rf_bin > 0) rf / max_rf_bin else NA_real_

interp <- c("highly similar",
            "moderately similar",
            "substantial disagreement",
            "essentially no shared structure")[
  findInterval(rf_norm, c(0.2, 0.5, 0.9)) + 1
]

cat(sprintf("[%s vs %s] common=%d  RF raw=%d/%d  norm=%.4f  (%s)\n",
            label_a, label_b, n_common, rf, max_rf_bin, rf_norm, interp))

# Bootstrap-filtered RF at each requested threshold. Reported with two
# normalizations: binary-max (= 2(n-3), same denominator as raw RF, lets
# you see absolute drop when low-support edges are discarded) and
# real-max (= |B(A)|+|B(B)| after filter, lets you see relative agreement
# among the surviving high-support splits).
boot_rows <- list()
for (t in boot_thresholds) {
  ba_t <- get_bipartitions(ta, threshold = t)
  bb_t <- get_bipartitions(tb, threshold = t)
  sh_t <- length(intersect(ba_t, bb_t))
  rf_t <- (length(ba_t) - sh_t) + (length(bb_t) - sh_t)
  max_real_t <- length(ba_t) + length(bb_t)
  norm_bin_t  <- if (max_rf_bin > 0) rf_t / max_rf_bin else NA_real_
  norm_real_t <- if (max_real_t > 0) rf_t / max_real_t else NA_real_
  boot_rows[[as.character(t)]] <- list(
    t = t, n_a = length(ba_t), n_b = length(bb_t),
    shared = sh_t, rf = rf_t,
    norm_bin = norm_bin_t, norm_real = norm_real_t)
  cat(sprintf("  boot>=%d:  bip A/B=%d/%d shared=%d  RF=%d  norm(bin)=%.3f  norm(real)=%.3f\n",
              t, length(ba_t), length(bb_t), sh_t, rf_t,
              norm_bin_t, norm_real_t))
}

# Major-clade Jaccard: |A ∩ B| / |A ∪ B| restricted to bipartitions
# where the SMALLER side has at least `min_size` tips. This intentionally
# ignores resolution at the single-leaf / tiny-clade level (where 20kbp
# windows give noisy signals) and focuses on whether large clusters of
# samples co-occur in both trees.
jacc_rows <- list()
for (k in jacc_sizes) {
  ca <- get_bipartitions(ta, threshold = -Inf, min_size = k)
  cb <- get_bipartitions(tb, threshold = -Inf, min_size = k)
  inter <- length(intersect(ca, cb))
  un <- length(union(ca, cb))
  j <- if (un > 0) inter / un else NA_real_
  jacc_rows[[as.character(k)]] <- list(
    k = k, n_a = length(ca), n_b = length(cb), inter = inter, un = un, j = j)
  cat(sprintf("  Jaccard size>=%d:  A=%d B=%d shared=%d union=%d -> J=%.3f\n",
              k, length(ca), length(cb), inter, un, j))
}

# Adjusted Rand Index at K clusters: cut each tree into K subtrees and
# compute how much the two clusterings agree. This is a "soft" similarity
# that does not require clade memberships to match EXACTLY -- it's
# robust to leaf-level shuffling as long as the broad groupings overlap.
#
# Cutting strategy: greedy "balanced" K-cut. Start with all tips in one
# cluster; iteratively pick the bipartition whose edge sits inside some
# current cluster AND whose removal splits that cluster most evenly.
# Score = min(left_size, right_size) after the split. This is robust to
# arbitrary rooting: when an outgroup is non-monophyletic in the tree
# (e.g. chr4 W1 has chimp alone at the root, bonobo deep inside the
# ingroup), the artificial root-edge gets a low score (small_side = 1)
# and is correctly skipped in favour of biology-meaningful splits.
cut_tree_k <- function(tree, k) {
  n_tips <- length(tree$tip.label)
  tips   <- tree$tip.label
  if (k < 1) stop("k must be >= 1")
  if (k == 1) {
    out <- rep(1L, n_tips); names(out) <- tips; return(out)
  }

  # Pre-compute all bipartitions as tip-name character sets (one side per
  # bipartition; the other side = setdiff(tips, side)).
  pp <- prop.part(tree)
  bip_sides <- lapply(pp, function(p) tips[p])

  clusters <- list(tips)
  while (length(clusters) < k) {
    best_score <- -1L
    best_ci    <- NA_integer_
    best_small <- NULL
    best_big   <- NULL
    for (ci in seq_along(clusters)) {
      cl <- clusters[[ci]]
      if (length(cl) < 2) next
      for (side in bip_sides) {
        a <- intersect(cl, side)
        if (length(a) == 0 || length(a) == length(cl)) next
        b <- setdiff(cl, a)
        score <- min(length(a), length(b))
        if (score > best_score) {
          best_score <- score
          best_ci    <- ci
          best_small <- if (length(a) <= length(b)) a else b
          best_big   <- if (length(a) <= length(b)) b else a
        }
      }
    }
    if (is.na(best_ci)) break  # nothing left to split
    clusters[[best_ci]] <- best_big
    clusters[[length(clusters) + 1]] <- best_small
  }

  assignment <- rep(NA_integer_, n_tips)
  names(assignment) <- tips
  for (i in seq_along(clusters)) {
    assignment[clusters[[i]]] <- i
  }
  assignment
}

ari_score <- function(c1, c2) {
  # Adjusted Rand Index from contingency table.
  c2 <- c2[names(c1)]
  n <- length(c1)
  ct <- table(c1, c2)
  if (length(ct) <= 1) return(NA_real_)
  sum_nij_2 <- sum(choose(ct, 2))
  a_i <- rowSums(ct)
  b_j <- colSums(ct)
  sum_a_2 <- sum(choose(a_i, 2))
  sum_b_2 <- sum(choose(b_j, 2))
  n_2 <- choose(n, 2)
  expected <- (sum_a_2 * sum_b_2) / n_2
  max_idx  <- 0.5 * (sum_a_2 + sum_b_2)
  if (isTRUE(all.equal(max_idx, expected))) return(NA_real_)
  (sum_nij_2 - expected) / (max_idx - expected)
}

ari_rows <- list()
for (k in ari_k_values) {
  if (k < 2 || k > n_common) {
    cat(sprintf("  ARI K=%d: skipped (out of range)\n", k))
    next
  }
  cl_a <- cut_tree_k(ta, k)
  cl_b <- cut_tree_k(tb, k)
  ari_val <- ari_score(cl_a, cl_b)
  sz_a <- sort(as.integer(table(cl_a)), decreasing = TRUE)
  sz_b <- sort(as.integer(table(cl_b)), decreasing = TRUE)
  ari_rows[[as.character(k)]] <- list(
    k = k, ari = ari_val,
    sizes_a = paste(sz_a, collapse = ","),
    sizes_b = paste(sz_b, collapse = ","))
  cat(sprintf("  ARI K=%d:  A sizes=[%s]  B sizes=[%s]  ARI=%.3f\n",
              k,
              paste(sz_a, collapse = ","),
              paste(sz_b, collapse = ","),
              ari_val))
}

# ---- text report --------------------------------------------------------
rep_p <- paste0(out_pref, ".report.txt")
fh <- file(rep_p, open = "w")
writeLines(c(
  paste0("# Tree comparison report: ", label_a, " vs ", label_b),
  "",
  paste0("Input ", label_a, ": ", tree_a_p),
  paste0("Input ", label_b, ": ", tree_b_p),
  paste0("Chrom for cleaning: ",
         ifelse(is.na(chrom), "(none -- inputs treated as cleaned)",
                chrom)),
  "",
  paste0("Tips in ", label_a, " (after cleaning): ", n_a_orig),
  paste0("Tips in ", label_b, " (after cleaning): ", n_b_orig),
  paste0("Tips only in ", label_a, " (dropped before RF): ",
         length(only_a)),
  paste0("Tips only in ", label_b, " (dropped before RF): ",
         length(only_b)),
  paste0("Tips in common: ", n_common),
  "",
  paste0("Bipartitions in ", label_a, ":            ", n_a_bip),
  paste0("Bipartitions in ", label_b, ":            ", n_b_bip),
  paste0("Shared bipartitions:                ", shared),
  paste0("Robinson-Foulds distance (raw):     ", rf),
  paste0("Max RF, fully-binary case [2(n-3)]: ", max_rf_bin),
  paste0("Max RF, actual (= bipA + bipB):     ", max_rf_real),
  paste0("Normalized RF (RF / binary max):    ",
         sprintf("%.4f", rf_norm)),
  paste0("Interpretation:                     ", interp),
  "",
  "# Bootstrap-filtered RF (drops splits whose support < threshold).",
  "# norm(bin)  = RF / 2(n-3)   -- absolute scale, same denominator as raw RF",
  "# norm(real) = RF / (|B(A)|+|B(B)|) after filter -- relative agreement",
  "#               among the splits that survived the filter.",
  sprintf("# %-14s  %5s  %5s  %5s  %5s  %10s  %10s",
          "threshold", "|B(A)|", "|B(B)|", "share", "RF", "norm(bin)",
          "norm(real)"),
  paste(sapply(boot_rows, function(r)
    sprintf("  boot>=%-3d       %5d  %5d  %5d  %5d  %10.4f  %10.4f",
            r$t, r$n_a, r$n_b, r$shared, r$rf, r$norm_bin, r$norm_real)),
    collapse = "\n"),
  "",
  "# Major-clade Jaccard (focuses on large splits, ignores leaf noise).",
  "# Bipartitions whose smaller side has < min_size tips are excluded,",
  "# then Jaccard = |A & B| / |A | B|.",
  sprintf("# %-13s  %5s  %5s  %5s  %5s  %7s",
          "min_size", "|A|", "|B|", "share", "union", "J"),
  paste(sapply(jacc_rows, function(r)
    sprintf("  size>=%-3d      %5d  %5d  %5d  %5d  %7.3f",
            r$k, r$n_a, r$n_b, r$inter, r$un, r$j)),
    collapse = "\n"),
  "",
  "# Adjusted Rand Index at K clusters (soft cluster agreement).",
  "# Each tree is cut into K subtrees; ARI measures clustering overlap.",
  "# 1.0 = identical clusterings, 0.0 = random, negative = anti-correlated.",
  "# A/B sizes show each cluster's tip count (largest first).",
  paste(sapply(ari_rows, function(r)
    sprintf("  K=%-2d  ARI=%6.3f   A=[%s]   B=[%s]",
            r$k, r$ari, r$sizes_a, r$sizes_b)),
    collapse = "\n"),
  "",
  "Reference scale (normalized RF):",
  "  0.0       topologically identical",
  "  < 0.2     highly similar",
  "  0.2-0.5   moderately similar (major clades match)",
  "  0.5-0.9   substantial disagreement",
  "  ~ 1.0     essentially no shared structure"
), fh)
if (length(only_a) > 0)
  writeLines(c("", paste0("Tips only in ", label_a, ":"),
               paste0("  ", only_a)), fh)
if (length(only_b) > 0)
  writeLines(c("", paste0("Tips only in ", label_b, ":"),
               paste0("  ", only_b)), fh)
close(fh)
cat("wrote: ", rep_p, "\n", sep = "")

# ---- tanglegram (ggtree-based 3-panel) ----------------------------------
# Left tree | cross-tree links | right tree (mirrored). Same layout
# approach as scripts/tanglegram.R but stripped to the bare essentials
# (no annotation strips, no cenhap). Lays tips out with cladogram
# spacing so polytomy-heavy IQ-TREE outputs don't collapse to a line.

ta_p <- ta  # already pruned to intersection above
tb_p <- tb

gt_a <- ggtree(ta_p, ladderize = TRUE, branch.length = "none")
gt_b <- ggtree(tb_p, ladderize = TRUE, branch.length = "none")
dl <- gt_a$data
dr <- gt_b$data

# Mirror right tree: root ends up at right edge of its panel.
max_xr_orig <- max(dr$x, na.rm = TRUE)
dr$x <- max_xr_orig - dr$x

build_edges <- function(d) {
  e <- d[d$parent != d$node, ]
  par_row <- d[match(e$parent, d$node), c("x", "y")]
  names(par_row) <- c("p_x", "p_y")
  cbind(e[, c("node", "x", "y")], par_row)
}
el <- build_edges(dl)
er <- build_edges(dr)

tips_l <- dl[dl$isTip, c("label", "x", "y")]
tips_r <- dr[dr$isTip, c("label", "x", "y")]
max_xl <- max(dl$x, na.rm = TRUE)
max_xr <- max(dr$x, na.rm = TRUE)
max_y  <- n_common

# tunables
tip_size     <- 2.5
tip_pad      <- max_xl * 0.02
label_zone_l <- max_xl * 0.55
label_zone_r <- max_xr * 0.55
y_lims  <- c(0.5, max_y + max_y * 0.04)
y_scale <- scale_y_continuous(limits = y_lims, expand = c(0, 0))

p_left <- ggplot() +
  geom_segment(data = el, aes(x = p_x, xend = x, y = y, yend = y),
               linewidth = 0.2) +
  geom_segment(data = el, aes(x = p_x, xend = p_x, y = p_y, yend = y),
               linewidth = 0.2) +
  geom_text(data = tips_l, aes(x = x + tip_pad, y = y, label = label),
            color = "black", hjust = 0, size = tip_size) +
  annotate("text", x = max_xl / 2, y = max_y + max_y * 0.025,
           label = label_a, size = 4, fontface = "bold") +
  y_scale +
  scale_x_continuous(limits = c(-max_xl * 0.02, max_xl + label_zone_l),
                     expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() + theme(plot.margin = margin(8, 2, 8, 8))

p_right <- ggplot() +
  geom_segment(data = er, aes(x = p_x, xend = x, y = y, yend = y),
               linewidth = 0.2) +
  geom_segment(data = er, aes(x = p_x, xend = p_x, y = p_y, yend = y),
               linewidth = 0.2) +
  geom_text(data = tips_r, aes(x = x - tip_pad, y = y, label = label),
            color = "black", hjust = 1, size = tip_size) +
  annotate("text", x = max_xr / 2, y = max_y + max_y * 0.025,
           label = label_b, size = 4, fontface = "bold") +
  y_scale +
  scale_x_continuous(limits = c(-label_zone_r, max_xr + max_xr * 0.02),
                     expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() + theme(plot.margin = margin(8, 8, 8, 2))

# Links between matching tips (label-keyed). Both trees have the same
# tip set after intersection.
links <- data.frame(
  label = ta_p$tip.label,
  y_l   = tips_l$y[match(ta_p$tip.label, tips_l$label)],
  y_r   = tips_r$y[match(ta_p$tip.label, tips_r$label)]
)
p_link <- ggplot(links) +
  geom_segment(aes(x = 0, xend = 1, y = y_l, yend = y_r),
               color = "black", alpha = 0.4, linewidth = 0.18) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  y_scale +
  coord_cartesian(clip = "off") +
  theme_void() + theme(plot.margin = margin(8, 0, 8, 0))

combined <- (p_left | p_link | p_right) +
  plot_layout(widths = c(1.0, 0.5, 1.0))
combined <- combined +
  plot_annotation(title = sprintf(
    "%s vs %s  --  RF = %d / %d  (norm %.3f, %s)",
    label_a, label_b, rf, max_rf_bin, rf_norm, interp),
    theme = theme(plot.title = element_text(hjust = 0.5,
                                            face = "bold", size = 11)))

pdf_p <- paste0(out_pref, ".tanglegram.pdf")
W <- 18
H <- max(10, 0.15 * n_common)
ggsave(pdf_p, combined, width = W, height = H, limitsize = FALSE)
cat("wrote: ", pdf_p, "  (", W, " x ", round(H, 1), " in)\n", sep = "")
