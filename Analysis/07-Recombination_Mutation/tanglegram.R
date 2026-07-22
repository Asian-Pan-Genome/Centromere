#!/usr/bin/env Rscript
# tanglegram.R
#
# Build a multi-panel tanglegram comparing the two MSA-window phylogenies
# (typically the 20kbp windows left/right of the centromere) of one
# chromosome, with sample-level annotation strips and TWO cenhap tracks
# (one ordered by each tree).
#
# Panel layout (left to right):
#   1. LEFT  tree           (tip labels + bootstrap text >= threshold)
#   2. superpopulation strip
#   3. project strip
#   4. cenhap  LEFT-aligned  (each sample anchored at its min_hor_start,
#                             ordered by LEFT tree's tip y)
#   5. cross-tree link lines
#   6. cenhap  RIGHT end-aligned (each sample anchored at its max_hor_end,
#                             axis reversed so 0 Mb is on the right,
#   7. RIGHT tree           (tip labels + bootstrap text >= threshold)
#
# Tree-cleaning, renaming, rooting, intersection -- same as before.
#
# Cenhap source: --hor-source {hicat,hormon}
#   hicat : <hicat-dir>/<chrom>.HiCAT.horstv.bed
#           + <color-map> (all.horstv.hicat.update.color.txt)
#           legend uses horstv_index column directly.
#   hormon: <hormon-dir>/<chrom>.graph.hordecomposition.final.xls
#           + <hor-index> (HOR_HORindex.txt) for hor_mn -> HOR_index
#           + HiCAT bed (for the per-hor_mn color override the v2
#             reference script applies on top of horstv_color).
#
# Outgroups (bonobo_<chrom>, chimp_<chrom>) get NO annotation strip
# nor cenhap row (their positions stay blank).
#
# Outputs (per chromosome and source):
#   <chrom>_left.cleaned.nwk
#   <chrom>_right.cleaned.nwk
#   <chrom>_tanglegram.<source>.pdf
#   <chrom>_tanglegram.<source>.preview.png
#   <chrom>_cenhap_missing.<source>.tsv  (samples in tree but no cenhap rows)
#
# Usage:
#   Rscript scripts/tanglegram.R <chrom> \
#       --left-tree     <left.treefile> \
#       --right-tree    <right.treefile> \
#       [--left-dropped <f>]  [--right-dropped <f>] \
#       [--left-mismatch <f>] [--right-mismatch <f>] \
#       [--cent-chrom <f>]  [--pop <f>] \
#       [--hicat-dir <dir>] [--hormon-dir <dir>] \
#       [--color-map <f>]   [--hor-index <f>] \
#       [--hor-source hicat|hormon]  (default hicat) \
#       [--top-n-hor 10]    [--boot-thresh 60] \
#       [--outdir <dir>]

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ---- argument parsing ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_named <- function(args, key, default = NA_character_) {
  i <- match(key, args)
  if (is.na(i) || i + 1 > length(args)) default else args[i + 1]
}
get_named_num <- function(args, key, default) {
  v <- get_named(args, key, NA_character_)
  if (is.na(v)) default else as.numeric(v)
}
if (length(args) < 1) {
  stop("Usage: Rscript tanglegram.R <chrom> --left-tree <f> --right-tree <f> ...")
}
chrom         <- args[1]
left_tree_p   <- get_named(args, "--left-tree")
right_tree_p  <- get_named(args, "--right-tree")
left_drop_p   <- get_named(args, "--left-dropped")
right_drop_p  <- get_named(args, "--right-dropped")
left_mis_p    <- get_named(args, "--left-mismatch")
right_mis_p   <- get_named(args, "--right-mismatch")
cent_chrom_p  <- get_named(args, "--cent-chrom",
                           "../cenhaps/cent_chrom.txt")
pop_p         <- get_named(args, "--pop",
                           "../cenhaps/populaion.xls")
hicat_dir     <- get_named(args, "--hicat-dir",
                           "../cenhaps/HiCAT_HOR")
hormon_dir    <- get_named(args, "--hormon-dir",
                           "../cenhaps/Hormon_HOR")
color_map_p   <- get_named(args, "--color-map",
                           "../cenhaps/all.horstv.hicat.update.color.txt")
hor_index_p   <- get_named(args, "--hor-index",
                           "../cenhaps/Hormon_HOR/HOR_HORindex.txt")
hor_source    <- get_named(args, "--hor-source", "hicat")
top_n_hor     <- get_named_num(args, "--top-n-hor", 10)
boot_thresh   <- get_named_num(args, "--boot-thresh", 60)
outdir        <- get_named(args, "--outdir", ".")

if (is.na(left_tree_p) || is.na(right_tree_p)) {
  stop("--left-tree and --right-tree are required")
}
if (!hor_source %in% c("hicat", "hormon")) {
  stop("--hor-source must be 'hicat' or 'hormon'")
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Fixed color palettes
suppop_cols  <- c(AFR = "#319b62", AMR = "#939393", EAS = "#D22E77",
                  EUR = "#0070c0", SAS = "#893f8b")
project_cols <- c(APG = "#E85827", HGSVC = "#285943",
                  HPRC = "#008C8C", Ref = "#969696")

# ---- helpers (tree cleaning) ---------------------------------------------
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
  if (!"status" %in% names(df) || !"aln_name" %in% names(df)) {
    warning("contig_check tsv missing aln_name/status columns: ", path)
    return(character(0))
  }
  miss <- df$aln_name[df$status == "mismatch"]
  vapply(miss, normalize_iqtree, character(1), USE.NAMES = FALSE)
}

rename_leaf <- function(label, chrom) {
  if (grepl("at_hsa", label))         return(paste0("bonobo_", chrom))
  if (grepl("hap[0-9]+_hsa", label))  return(paste0("chimp_",  chrom))
  if (grepl("^chr[0-9XYM]+_[0-9]", label) && !grepl("_hsa", label)) {
    return(paste0("CHM13_", chrom))
  }
  parts <- strsplit(label, "_", fixed = TRUE)[[1]]
  hap_set <- c("Mat", "Pat", "hap1", "hap2")
  hap_idx <- which(parts %in% hap_set)[1]
  if (!is.na(hap_idx) && hap_idx > 1) {
    sample <- paste(parts[seq_len(hap_idx - 1)], collapse = "_")
    return(paste0(sample, "_", parts[hap_idx]))
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

rename_tree_tips <- function(tr, chrom, side) {
  new <- vapply(tr$tip.label, rename_leaf, character(1),
                chrom = chrom, USE.NAMES = FALSE)
  na_idx <- which(is.na(new))
  if (length(na_idx) > 0) {
    message("WARN [", side, "]: ", length(na_idx),
            " tips unclassified, keeping original label. First few:")
    for (i in head(na_idx, 5)) message("    ", tr$tip.label[i])
    new[na_idx] <- tr$tip.label[na_idx]
  }
  if (any(duplicated(new))) {
    seen <- list()
    for (j in seq_along(new)) {
      lab <- new[j]
      if (sum(new == lab) > 1) {
        seen[[lab]] <- if (is.null(seen[[lab]])) 1 else seen[[lab]] + 1
        if (seen[[lab]] > 1) new[j] <- paste0(lab, "_", seen[[lab]])
      }
    }
    message("WARN [", side, "]: duplicate renamed tips disambiguated")
  }
  tr$tip.label <- new
  tr
}

root_on_outgroups <- function(tr, chrom, side) {
  og <- intersect(c(paste0("bonobo_", chrom), paste0("chimp_", chrom)),
                  tr$tip.label)
  if (length(og) == 0) {
    message("WARN [", side, "]: no outgroup tip found, leaving unrooted")
    return(tr)
  }
  tryCatch(
    root(tr, outgroup = og, resolve.root = TRUE),
    error = function(e) {
      message("WARN [", side, "]: outgroup not monophyletic (",
              conditionMessage(e), "), leaving unrooted")
      tr
    }
  )
}

reorder_og_top <- function(tr, chrom_arg) {
  og <- intersect(c(paste0("bonobo_", chrom_arg),
                    paste0("chimp_", chrom_arg)),
                  tr$tip.label)
  if (length(og) == 0) return(tr)
  rotateConstr(tr, c(setdiff(tr$tip.label, og), og))
}

# ---- helpers (annotation) ------------------------------------------------

load_annot <- function(cent_p, pop_p, chrom_arg) {
  if (!file.exists(cent_p)) stop("missing cent_chrom: ", cent_p)
  if (!file.exists(pop_p))  stop("missing population: ",  pop_p)
  df_cen <- read.csv(cent_p, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)
  df_cen <- df_cen[df_cen$chrom == chrom_arg, ]
  if (nrow(df_cen) == 0) stop("no cent_chrom rows for chrom=", chrom_arg)
  df_pop <- read.csv(pop_p, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)
  df_pop$superpopulation <- ifelse(df_pop$superpopulation == "EAS-APG",
                                    "EAS", df_pop$superpopulation)
  df_cen %>%
    left_join(df_pop %>% select(sample, superpopulation), by = "sample") %>%
    mutate(
      superpopulation = ifelse(superpopulation == "EAS-APG", "EAS",
                                superpopulation),
      tip_label = ifelse(sample == "CHM13",
                         paste0("CHM13_", chrom_arg), sample_hap)
    ) %>%
    select(tip_label, sample, hap, project, superpopulation,
           cen_start = start, cen_end = end)
}

# Common back-end: given a long-format `stv` table with columns
#   tip_label, horstart, horend, legend_class, plot_color
# join in cen_start/cen_end, clip to centromere, compute per-sample
# min_hor_start anchor, and pick top-N legend_class by mean per-sample
# total length (excluding `exclude_set`).
finalize_cenhap <- function(stv, annot_df, top_n, exclude_set) {
  ann_sel <- annot_df %>% select(tip_label, cen_start, cen_end)
  rects <- stv %>%
    inner_join(ann_sel, by = "tip_label") %>%
    filter(horstart >= cen_start, horend <= cen_end) %>%
    group_by(tip_label) %>%
    mutate(min_hor_start = min(horstart)) %>%
    ungroup() %>%
    mutate(modstart   = horstart - min_hor_start + 1,
           modend     = horend   - min_hor_start + 1,
           hor_length = horend   - horstart)
  if (nrow(rects) == 0) {
    message("WARN: no cenhap rows fall inside cen_start/cen_end")
    return(NULL)
  }
  per_sample_len <- rects %>%
    filter(!legend_class %in% exclude_set) %>%
    group_by(tip_label, legend_class) %>%
    summarise(sample_total = sum(hor_length), .groups = "drop")
  avg_by_class <- per_sample_len %>%
    group_by(legend_class) %>%
    summarise(mean_total = mean(sample_total), .groups = "drop") %>%
    arrange(desc(mean_total))
  top_set <- head(avg_by_class$legend_class, top_n)
  # Per-class representative color, restricted to top_set. The bars
  # themselves are drawn from `plot_color` directly (identity scale),
  # so non-top10 HORstv classes KEEP their original colors -- only the
  # legend is collapsed to just the top N.
  color_lookup <- rects %>%
    filter(legend_class %in% top_set) %>%
    distinct(legend_class, plot_color) %>%
    arrange(match(legend_class, top_set))
  list(rects = rects, top_set = top_set, color_lookup = color_lookup)
}

# HiCAT source: <hicat-dir>/<chrom>.HiCAT.horstv.bed  (8 cols)
load_cenhap_hicat <- function(hicat_dir, color_map_p, chrom_arg,
                              annot_df, valid_tips, top_n) {
  bed_p <- file.path(hicat_dir, paste0(chrom_arg, ".HiCAT.horstv.bed"))
  if (!file.exists(bed_p)) {
    message("WARN: no HiCAT bed at ", bed_p, " -- cenhap disabled")
    return(NULL)
  }
  data <- read.csv(bed_p, header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE)
  colnames(data) <- c("sample_hap_chrom", "horstart", "horend", "hor_mn",
                      "horclass", "horclass_color", "horstv_index",
                      "horstv_color")
  recolordf <- read.csv(color_map_p, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
  colnames(recolordf) <- c("horstv_index", "horstv_newcolor")

  stv <- data %>%
    mutate(tip_label = case_when(
      sample_hap_chrom == chrom_arg ~ paste0("CHM13_", chrom_arg),
      str_detect(sample_hap_chrom, "#") ~ sample_hap_chrom %>%
        str_replace("#", "_") %>%
        str_replace(paste0("#", chrom_arg), ""),
      str_detect(sample_hap_chrom, "_") ~ str_replace(sample_hap_chrom,
                                                       paste0("_", chrom_arg), ""),
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(tip_label), tip_label %in% valid_tips) %>%
    left_join(recolordf, by = "horstv_index") %>%
    mutate(plot_color = ifelse(horstv_color == "-", "#525252",
                                ifelse(!is.na(horstv_newcolor),
                                        horstv_newcolor, horstv_color)),
           legend_class = horstv_index) %>%
    select(tip_label, horstart, horend, legend_class, plot_color)

  if (nrow(stv) == 0) {
    message("WARN: no HiCAT cenhap rows after tip filter")
    return(NULL)
  }
  finalize_cenhap(stv, annot_df, top_n, c("Mon", "Others"))
}

# HORmon source: <hormon-dir>/<chrom>.graph.hordecomposition.final.xls
#                + <hor-index>  (HOR_HORindex.txt for hor_mn -> HOR_index)
#                + HiCAT bed    (for the per-hor_mn color override
#                                described by the v2 reference script)
load_cenhap_hormon <- function(hormon_dir, hicat_dir, color_map_p,
                               hor_index_p, chrom_arg,
                               annot_df, valid_tips, top_n) {
  hormon_p <- file.path(hormon_dir,
                        paste0(chrom_arg, ".graph.hordecomposition.final.xls"))
  if (!file.exists(hormon_p)) {
    message("WARN: no HORmon xls at ", hormon_p, " -- cenhap disabled")
    return(NULL)
  }
  if (!file.exists(hor_index_p)) {
    message("WARN: no HOR_HORindex.txt at ", hor_index_p, " -- cenhap disabled")
    return(NULL)
  }
  gdata <- read.csv(hormon_p, header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE)
  colnames(gdata) <- c("sample_hap_chrom", "horstart", "horend", "hor_mn",
                       "horclass", "horclass_color", "horstv_color")
  hor_idx <- read.csv(hor_index_p, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
  # accept either HOR/HOR_index or hor_mn/horstv_index
  if (all(c("HOR", "HOR_index") %in% names(hor_idx))) {
    hor_idx <- hor_idx %>% rename(hor_mn = HOR, legend_class = HOR_index)
  } else if (all(c("hor_mn", "horstv_index") %in% names(hor_idx))) {
    hor_idx <- hor_idx %>% rename(legend_class = horstv_index)
  } else {
    stop("HOR_HORindex.txt columns must be HOR/HOR_index or hor_mn/horstv_index")
  }

  # Build hor_mn -> color override from HiCAT (recolored), matching the
  # v2 reference script's `newcolor` logic.
  hicat_bed_p <- file.path(hicat_dir,
                           paste0(chrom_arg, ".HiCAT.horstv.bed"))
  override <- data.frame(hor_mn = character(0), override_color = character(0))
  if (file.exists(hicat_bed_p) && file.exists(color_map_p)) {
    h <- read.csv(hicat_bed_p, header = FALSE, sep = "\t",
                  stringsAsFactors = FALSE)
    colnames(h) <- c("sample_hap_chrom", "horstart", "horend", "hor_mn",
                     "horclass", "horclass_color", "horstv_index",
                     "horstv_color")
    rc <- read.csv(color_map_p, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)
    colnames(rc) <- c("horstv_index", "horstv_newcolor")
    override <- h %>%
      left_join(rc, by = "horstv_index") %>%
      mutate(override_color = ifelse(!is.na(horstv_newcolor),
                                      horstv_newcolor, horstv_color)) %>%
      filter(grepl("_", hor_mn)) %>%
      select(hor_mn, override_color) %>%
      distinct(hor_mn, .keep_all = TRUE)
  } else {
    message("NOTE: HiCAT bed or color map missing; HORmon colors won't be ",
            "overridden (using horstv_color from the xls directly).")
  }

  stv <- gdata %>%
    mutate(tip_label = case_when(
      sample_hap_chrom == chrom_arg ~ paste0("CHM13_", chrom_arg),
      str_detect(sample_hap_chrom, "#") ~ sample_hap_chrom %>%
        str_replace("#", "_") %>%
        str_replace(paste0("#", chrom_arg), ""),
      str_detect(sample_hap_chrom, "_") ~ str_replace(sample_hap_chrom,
                                                       paste0("_", chrom_arg), ""),
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(tip_label), tip_label %in% valid_tips) %>%
    left_join(hor_idx, by = "hor_mn") %>%
    left_join(override, by = "hor_mn") %>%
    mutate(plot_color = ifelse(!is.na(override_color), override_color,
                                ifelse(horstv_color == "-", "#525252",
                                        horstv_color)),
           # final sanitize: any remaining "-" or NA (e.g., from a HiCAT
           # override that itself was "-") becomes the neutral gray.
           plot_color = ifelse(is.na(plot_color) | plot_color == "-",
                                "#525252", plot_color),
           legend_class = ifelse(is.na(legend_class) | legend_class == "",
                                  "Others", legend_class)) %>%
    select(tip_label, horstart, horend, legend_class, plot_color)

  if (nrow(stv) == 0) {
    message("WARN: no HORmon cenhap rows after tip filter")
    return(NULL)
  }
  finalize_cenhap(stv, annot_df, top_n, c("Mon", "Others"))
}

# ---- load + clean trees --------------------------------------------------
t_left  <- read.tree(left_tree_p)
t_right <- read.tree(right_tree_p)
message("loaded LEFT : ", left_tree_p,  "  (", length(t_left$tip.label),  " tips)")
message("loaded RIGHT: ", right_tree_p, "  (", length(t_right$tip.label), " tips)")

drop_left  <- unique(c(read_drop_list(left_drop_p),
                       read_mismatch_list(left_mis_p)))
drop_right <- unique(c(read_drop_list(right_drop_p),
                       read_mismatch_list(right_mis_p)))
drop_left_in  <- intersect(drop_left,  t_left$tip.label)
drop_right_in <- intersect(drop_right, t_right$tip.label)
message("LEFT  drop set: requested ", length(drop_left),
        ", matched in tree ", length(drop_left_in))
message("RIGHT drop set: requested ", length(drop_right),
        ", matched in tree ", length(drop_right_in))
if (length(drop_left_in)  > 0) t_left  <- drop.tip(t_left,  drop_left_in)
if (length(drop_right_in) > 0) t_right <- drop.tip(t_right, drop_right_in)

t_left  <- rename_tree_tips(t_left,  chrom, "LEFT")
t_right <- rename_tree_tips(t_right, chrom, "RIGHT")

t_left  <- root_on_outgroups(t_left,  chrom, "LEFT")
t_right <- root_on_outgroups(t_right, chrom, "RIGHT")

common0 <- intersect(t_left$tip.label, t_right$tip.label)
left_only0  <- setdiff(t_left$tip.label,  t_right$tip.label)
right_only0 <- setdiff(t_right$tip.label, t_left$tip.label)
message("intersection: ", length(common0),
        " common | dropping ", length(left_only0), " left-only, ",
        length(right_only0), " right-only")
if (length(left_only0)  > 0) message("  left-only:  ",
                                     paste(left_only0,  collapse = ", "))
if (length(right_only0) > 0) message("  right-only: ",
                                     paste(right_only0, collapse = ", "))
if (length(left_only0)  > 0) t_left  <- drop.tip(t_left,  left_only0)
if (length(right_only0) > 0) t_right <- drop.tip(t_right, right_only0)

out_left  <- file.path(outdir, paste0(chrom, "_left.cleaned.nwk"))
out_right <- file.path(outdir, paste0(chrom, "_right.cleaned.nwk"))
write.tree(t_left,  out_left)
write.tree(t_right, out_right)
message("wrote: ", out_left,  " (", length(t_left$tip.label),  " tips)")
message("wrote: ", out_right, " (", length(t_right$tip.label), " tips)")

# ---- load annotation + cenhap (source-dependent) ------------------------
all_tips <- t_left$tip.label
annot_df <- load_annot(cent_chrom_p, pop_p, chrom)
message("annotation rows for ", chrom, ": ", nrow(annot_df),
        "  (matched tree tips: ",
        sum(all_tips %in% annot_df$tip_label), ")")

cenhap <- if (hor_source == "hicat") {
  load_cenhap_hicat(hicat_dir, color_map_p, chrom,
                    annot_df, all_tips, top_n_hor)
} else {
  load_cenhap_hormon(hormon_dir, hicat_dir, color_map_p, hor_index_p,
                     chrom, annot_df, all_tips, top_n_hor)
}
if (!is.null(cenhap)) {
  message("cenhap[", hor_source, "] rects: ", nrow(cenhap$rects),
          "  top-", top_n_hor, ": ",
          paste(head(cenhap$top_set, 5), collapse = ", "), " ...")
}

# Missing-sample report: tips that exist in the (non-outgroup) tree but
# have no cenhap rows after the cen-window filter. Both trees have the
# same tip set after intersection, so report once.
og_names <- c(paste0("bonobo_", chrom), paste0("chimp_", chrom))
expected_tips <- setdiff(all_tips, og_names)
present_tips  <- if (is.null(cenhap)) {
  character(0)
} else {
  unique(cenhap$rects$tip_label)
}
missing_tips  <- setdiff(expected_tips, present_tips)
if (length(missing_tips) > 0) {
  miss_df <- data.frame(tip_label = missing_tips) %>%
    left_join(annot_df %>% select(tip_label, sample, hap),
              by = "tip_label")
  miss_p <- file.path(outdir, paste0(chrom, "_cenhap_missing.",
                                     hor_source, ".tsv"))
  write.table(miss_df, miss_p, sep = "\t", quote = FALSE,
              row.names = FALSE)
  preview <- paste(head(missing_tips, 5), collapse = ", ")
  message("cenhap missing: ", length(missing_tips),
          " samples (first 5: ", preview, ")")
  message("wrote: ", miss_p)
} else {
  message("cenhap missing: 0 samples")
}

# ---- per-tree layout (cladogram) ----------------------------------------
t_left_t  <- reorder_og_top(t_left,  chrom)
t_right_t <- reorder_og_top(t_right, chrom)

gt_left  <- ggtree(t_left_t,  ladderize = FALSE, branch.length = "none")
gt_right <- ggtree(t_right_t, ladderize = FALSE, branch.length = "none")
dl <- gt_left$data
dr <- gt_right$data

n_tips_l <- length(t_left_t$tip.label)
n_tips_r <- length(t_right_t$tip.label)
max_y <- max(n_tips_l, n_tips_r)
if (n_tips_l < max_y) dl$y <- dl$y + (max_y - n_tips_l)
if (n_tips_r < max_y) dr$y <- dr$y + (max_y - n_tips_r)

# Mirror right tree so root ends at the right edge of its panel
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
tips_l$is_og <- tips_l$label %in% og_names
tips_r$is_og <- tips_r$label %in% og_names

max_xl <- max(dl$x, na.rm = TRUE)
max_xr <- max(dr$x, na.rm = TRUE)

get_support_pts <- function(d, thresh) {
  d2 <- d[!d$isTip, ]
  sup <- suppressWarnings(as.numeric(d2$label))
  ok  <- !is.na(sup) & sup > thresh
  data.frame(x = d2$x[ok], y = d2$y[ok], support = sup[ok])
}
sup_l <- get_support_pts(dl, boot_thresh)
sup_r <- get_support_pts(dr, boot_thresh)
message("bootstrap >", boot_thresh, " nodes -- LEFT: ", nrow(sup_l),
        ", RIGHT: ", nrow(sup_r))

# ---- shared y-axis -------------------------------------------------------
y_lims  <- c(0.5, max_y + max_y * 0.04)
y_scale <- scale_y_continuous(limits = y_lims, expand = c(0, 0))

# ---- panel: LEFT tree ---------------------------------------------------
tip_size      <- 2.8
boot_text_sz  <- 1.8
tip_pad_l     <- max_xl * 0.02
label_zone_l  <- max_xl * 0.55

p_left <- ggplot() +
  geom_segment(data = el, aes(x = p_x, xend = x, y = y, yend = y),
               linewidth = 0.2) +
  geom_segment(data = el, aes(x = p_x, xend = p_x, y = p_y, yend = y),
               linewidth = 0.2) +
  geom_text(data = tips_l,
            aes(x = x + tip_pad_l, y = y, label = label),
            color = "black", hjust = 0, size = tip_size) +
  {if (nrow(sup_l) > 0)
    geom_text(data = sup_l,
              aes(x = x, y = y, label = round(support)),
              size = boot_text_sz, color = "black",
              hjust = 1.15, vjust = -0.25)} +
  annotate("text", x = max_xl / 2, y = max_y + max_y * 0.025,
           label = paste0(chrom, " p-arm"),
           size = 4, fontface = "bold") +
  y_scale +
  scale_x_continuous(limits = c(-max_xl * 0.02, max_xl + label_zone_l),
                     expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(8, 2, 8, 8))

# ---- panel: superpop strip (ordered by LEFT tree y) ---------------------
tipy_l <- tips_l %>% select(label, y)
tipy_r <- tips_r %>% select(label, y)

ann_join_l <- annot_df %>%
  filter(tip_label %in% tipy_l$label, !tip_label %in% og_names) %>%
  inner_join(tipy_l %>% rename(tip_label = label), by = "tip_label")

p_suppop <- ggplot(ann_join_l,
                   aes(x = 1, y = y, fill = superpopulation)) +
  geom_tile(width = 1, height = 0.95) +
  scale_fill_manual(values = suppop_cols, name = "Superpopulation",
                    na.value = "transparent") +
  scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  y_scale +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(8, 1, 8, 1))

p_project <- ggplot(ann_join_l, aes(x = 1, y = y, fill = project)) +
  geom_tile(width = 1, height = 0.95) +
  scale_fill_manual(values = project_cols, name = "Project",
                    na.value = "transparent") +
  scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  y_scale +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(8, 1, 8, 1))

# ---- cenhap panels (left-aligned + right-aligned/end-anchored) -------------
# Both panels render the SAME data; what differs is (a) the tip y values
# they're joined to (LEFT tree y vs RIGHT tree y), and (b) the x-axis
# anchor: left panel aligns each sample's min_hor_start at 0 (normal x);
# right panel aligns each sample's max_hor_end at 0 (scale_x_reverse).
build_cenhap_panel <- function(cenhap, tipy, max_y, y_lims, og_names,
                                reverse_x = FALSE) {
  if (is.null(cenhap)) {
    return(ggplot() + theme_void())
  }
  rects <- cenhap$rects %>%
    filter(tip_label %in% tipy$label, !tip_label %in% og_names) %>%
    inner_join(tipy %>% rename(tip_label = label), by = "tip_label")
  if (nrow(rects) == 0) return(ggplot() + theme_void())

  # Right panel: re-anchor to max_hor_end so every sample's HOR-array
  # rightmost end aligns at the panel right edge (scale_x_reverse puts
  # 0 Mb on the right).
  if (reverse_x) {
    rects <- rects %>%
      group_by(tip_label) %>%
      mutate(max_hor_end = max(horend)) %>%
      ungroup() %>%
      mutate(modstart = max_hor_end - horend + 1,
             modend   = max_hor_end - horstart + 1)
  }

  # Bars are colored by the row's own plot_color (identity scale), so
  # every HORstv class keeps its original color. The legend is then
  # manually restricted to just the top N via scale_fill_identity's
  # breaks/labels (non-top entries are drawn but absent from the key).
  top_colors <- cenhap$color_lookup$plot_color[
    match(cenhap$top_set, cenhap$color_lookup$legend_class)]
  top_labels <- cenhap$top_set

  max_x <- max(rects$modend, na.rm = TRUE)
  x_scale <- if (reverse_x) {
    scale_x_reverse(breaks = seq(0, max_x, by = 1e6),
                    labels = function(x) paste0(x / 1e6, " Mb"),
                    expand = c(0.005, 0))
  } else {
    scale_x_continuous(breaks = seq(0, max_x, by = 1e6),
                       labels = function(x) paste0(x / 1e6, " Mb"),
                       expand = c(0.005, 0))
  }

  ggplot(rects) +
    geom_rect(aes(xmin = modstart, xmax = modend,
                  ymin = y - 0.45, ymax = y + 0.45,
                  fill = plot_color), color = NA) +
    scale_fill_identity(name = "HORs", guide = "legend",
                        breaks = top_colors,
                        labels = top_labels) +
    x_scale +
    scale_y_continuous(limits = y_lims, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(axis.title    = element_blank(),
          axis.text.y   = element_blank(),
          axis.ticks.y  = element_blank(),
          axis.line.y   = element_blank(),
          axis.text.x   = element_text(size = 7),
          plot.margin   = margin(8, 2, 8, 2))
}
p_cenhap_l <- build_cenhap_panel(cenhap, tipy_l, max_y, y_lims, og_names,
                                  reverse_x = FALSE)
p_cenhap_r <- build_cenhap_panel(cenhap, tipy_r, max_y, y_lims, og_names,
                                  reverse_x = TRUE)

# ---- panel: cross-tree links --------------------------------------------
common <- intersect(t_left_t$tip.label, t_right_t$tip.label)
links <- data.frame(
  label = common,
  y_l   = tips_l$y[match(common, tips_l$label)],
  y_r   = tips_r$y[match(common, tips_r$label)]
)
p_link <- ggplot(links) +
  geom_segment(aes(x = 0, xend = 1, y = y_l, yend = y_r),
               color = "black", alpha = 0.5, linewidth = 0.15) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  y_scale +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(8, 0, 8, 0))

# ---- panel: RIGHT tree --------------------------------------------------
tip_pad_r    <- max_xr * 0.02
label_zone_r <- max_xr * 0.55

p_right <- ggplot() +
  geom_segment(data = er, aes(x = p_x, xend = x, y = y, yend = y),
               linewidth = 0.2) +
  geom_segment(data = er, aes(x = p_x, xend = p_x, y = p_y, yend = y),
               linewidth = 0.2) +
  geom_text(data = tips_r,
            aes(x = x - tip_pad_r, y = y, label = label),
            color = "black", hjust = 1, size = tip_size) +
  {if (nrow(sup_r) > 0)
    geom_text(data = sup_r,
              aes(x = x, y = y, label = round(support)),
              size = boot_text_sz, color = "black",
              hjust = -0.15, vjust = -0.25)} +
  annotate("text", x = max_xr / 2, y = max_y + max_y * 0.025,
           label = paste0(chrom, " q-arm"),
           size = 4, fontface = "bold") +
  y_scale +
  scale_x_continuous(limits = c(-label_zone_r, max_xr + max_xr * 0.02),
                     expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(8, 8, 8, 2))

# ---- assemble -----------------------------------------------------------
w_tree   <- 1.0
w_strip  <- 0.07
w_gap    <- 0.05
w_cenhap <- 1.4
w_link   <- 0.4
combined <- (p_left | plot_spacer() |
             p_suppop | plot_spacer() |
             p_project | plot_spacer() |
             p_cenhap_l | p_link | p_cenhap_r |
             plot_spacer() | p_right) +
  plot_layout(widths = c(w_tree, w_gap,
                         w_strip, w_gap,
                         w_strip, w_gap,
                         w_cenhap, w_link, w_cenhap,
                         w_gap, w_tree),
              guides = "collect") &
  theme(legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        legend.text     = element_text(size = 7),
        legend.title    = element_text(size = 8))

# ---- save ---------------------------------------------------------------
out_pdf <- file.path(outdir, paste0(chrom, "_tanglegram.",
                                    hor_source, ".pdf"))
H <- max(10, 0.1 * max_y)
W <- 32
ggsave(out_pdf, combined, width = W, height = H, limitsize = FALSE)
out_png <- file.path(outdir, paste0(chrom, "_tanglegram.",
                                    hor_source, ".preview.png"))
ggsave(out_png, combined, width = W, height = H, limitsize = FALSE,
       dpi = 60)
message("wrote: ", out_pdf, "  (", W, " x ", round(H, 1), " in)")
message("wrote: ", out_png, " (preview)")
