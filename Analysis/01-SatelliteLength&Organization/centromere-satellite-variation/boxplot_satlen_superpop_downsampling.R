rm(list = ls())

library(tidyverse)
library(ggplot2)
library(patchwork)

# NOTE: do NOT use showtext here.
# showtext::showtext_auto() draws every glyph as vector outlines (paths),
# making text non-editable in Adobe Illustrator. Instead we let cairo_pdf
# embed the real Arial font, so all text stays editable.
# Arial is a system font on Windows; cairo_pdf resolves "Arial" by family name.
base_family <- "Arial"

# ── Data loading ──────────────────────────────────────────────────────────────
satlen <- read.csv("data/all.sat.length.xls", header = TRUE, sep = "\t")
df_cen <- read.csv("data/cent_chrom.txt",     header = TRUE, sep = "\t")
popdf  <- read.csv("data/populaion.xls",      header = TRUE, sep = "\t")

popdf <- popdf %>%
  mutate(superpopulation = ifelse(superpopulation == "EAS-APG", "EAS", superpopulation))

superpopulation_colors <- c(
  "AFR" = "#319b62", "AMR" = "#939393", "EAS" = "#d22e77",
  "EUR" = "#0070c0", "SAS" = "#893f8b"
)
pop_order   <- c("AFR", "AMR", "EAS", "EUR", "SAS")
chrom_order <- c(paste0("chr", 1:22), "chrX", "chrY")

# EAS vs others — bracket rank
pair_bracket_rank <- c("AMR" = 1L, "EUR" = 2L, "AFR" = 3L, "SAS" = 4L)

free_y_sats     <- c("HSat1", "HSat2", "HSat3")
MEAN_FILTER_BP  <- 500000   # 500 Kbp (HSat1/2/3 only)

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  # Plain text (NOT plotmath) so "p = xxx" stays one editable string in Illustrator
  if (p < 0.001) {
    paste0("p = ", formatC(p, digits = 1, format = "e"))
  } else {
    paste0("p = ", formatC(p, digits = 3, format = "f"))
  }
}

set.seed(42)
K     <- 10L
MIN_N <- 5L

# ── Build per-satellite data + Wilcoxon results ──────────────────────────────
target_dfs  <- list()
results_dfs <- list()

for (target_sat in c("ASat", "HSat1", "HSat2", "HSat3", "BSat")) {
  message("\n=== ", target_sat, " ===")

  if (target_sat == "HSat1") {
    target_df <- satlen %>%
      left_join(df_cen %>% select(sample_hap_chrom, filterflag),
                by = "sample_hap_chrom", relationship = "many-to-many") %>%
      filter(grepl("HSat1", satellite), filterflag == 0) %>%
      group_by(sample_hap_chrom, sample_hap, sample, hap, chrom, project, filterflag) %>%
      summarise(length = sum(length), .groups = "drop") %>%
      mutate(len_mbp = length / 1e6, satellite = "HSat1") %>%
      left_join(popdf %>% select(sample, superpopulation), by = "sample")
  } else {
    target_df <- satlen %>%
      left_join(df_cen %>% select(sample_hap_chrom, filterflag),
                by = "sample_hap_chrom", relationship = "many-to-many") %>%
      filter(satellite == target_sat, filterflag == 0) %>%
      mutate(len_mbp = length / 1e6) %>%
      left_join(popdf %>% select(sample, superpopulation), by = "sample")
  }

  target_df <- target_df %>%
    filter(superpopulation %in% pop_order) %>%
    mutate(
      superpopulation = factor(superpopulation, levels = pop_order),
      chrom           = factor(chrom, levels = chrom_order),
      x_pos           = as.integer(superpopulation)
    )

  if (nrow(target_df) == 0) { message("No data, skipping."); next }

  # Pairwise downsampling Wilcoxon for all chrom × pop pairs
  results_list <- list()
  for (tgt_chrom in levels(droplevels(target_df$chrom))) {
    subdf <- target_df %>% filter(chrom == tgt_chrom)
    pops <- subdf %>% count(superpopulation, .drop = TRUE) %>%
      filter(n > 0) %>% pull(superpopulation) %>% as.character()
    if (length(pops) < 2) next

    for (pair in combn(pops, 2, simplify = FALSE)) {
      g1 <- pair[1]; g2 <- pair[2]
      n1 <- sum(subdf$superpopulation == g1)
      n2 <- sum(subdf$superpopulation == g2)
      min_n <- min(n1, n2)

      if (min_n < MIN_N) {
        results_list[[length(results_list) + 1]] <- tibble(
          chrom = tgt_chrom, group1 = g1, group2 = g2,
          n1 = n1, n2 = n2, mean_p = NA_real_,
          status = paste0("skipped (min_n=", min_n, ")")
        )
        next
      }
      p_vals <- vapply(seq_len(K), function(.) {
        s1 <- subdf %>% filter(superpopulation == g1) %>% slice_sample(n = min_n)
        s2 <- subdf %>% filter(superpopulation == g2) %>% slice_sample(n = min_n)
        wilcox.test(s1$len_mbp, s2$len_mbp, exact = FALSE)$p.value
      }, numeric(1))
      results_list[[length(results_list) + 1]] <- tibble(
        chrom = tgt_chrom, group1 = g1, group2 = g2,
        n1 = n1, n2 = n2, mean_p = mean(p_vals),
        status = paste0("tested (downsample_n=", min_n, ", K=", K, ")")
      )
    }
  }

  results_df <- bind_rows(results_list)
  out_tsv <- paste0("downsampling_wilcox_", target_sat, ".tsv")
  write.table(results_df, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  message("TSV saved: ", out_tsv)

  target_dfs[[target_sat]]  <- target_df
  results_dfs[[target_sat]] <- results_df
}

# ── Plot builder ─────────────────────────────────────────────────────────────
# Builds a ggplot for one satellite, including:
#   • boxplots per superpopulation
#   • sample-size labels above each box
#   • EAS-vs-other significance brackets (mean_p)
make_plot <- function(target_sat,
                      target_df,
                      results_df,
                      facet_arg   = NULL,
                      show_title  = TRUE,
                      show_legend = TRUE,
                      show_y_lab  = TRUE,
                      n_label_size = 4.8,
                      p_label_size = 5.5,
                      axis_text_y_size = 13,
                      axis_text_x_size = 13,
                      axis_title_size  = 15,
                      strip_text_size  = 15) {

  # 1. Apply mean array size filter to all satellites (drop panels with mean < threshold)
  keep_chroms <- target_df %>%
    group_by(chrom) %>%
    summarise(mean_len = mean(length, na.rm = TRUE), .groups = "drop") %>%
    filter(mean_len >= MEAN_FILTER_BP) %>%
    pull(chrom) %>% as.character()
  plot_df <- target_df %>%
    filter(chrom %in% keep_chroms) %>%
    mutate(chrom = factor(chrom, levels = keep_chroms[keep_chroms %in% chrom_order]))
  message(target_sat,
          " chromosomes kept (mean >= ", MEAN_FILTER_BP/1e3, " Kbp): ",
          paste(levels(plot_df$chrom), collapse = ", "))
  if (nrow(plot_df) == 0) return(NULL)

  # 2. Sample-size labels per (chrom, pop)
  n_per_pop <- plot_df %>%
    group_by(chrom, superpopulation) %>%
    summarise(n_pop = n(), .groups = "drop") %>%
    filter(n_pop > 0) %>%
    mutate(x_pos = as.integer(superpopulation))

  # 3. Per-chromosome max for label positioning
  chrom_max_df <- plot_df %>%
    group_by(chrom) %>%
    summarise(y_data_max = max(len_mbp, na.rm = TRUE), .groups = "drop")

  # 4. EAS-vs-other annotations
  afr_annot <- results_df %>%
    filter(!is.na(mean_p)) %>%
    mutate(
      pop_other = case_when(
        group1 == "EAS" & group2 %in% names(pair_bracket_rank) ~ group2,
        group2 == "EAS" & group1 %in% names(pair_bracket_rank) ~ group1,
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(pop_other),
           chrom %in% as.character(levels(plot_df$chrom))) %>%
    mutate(
      bracket_rank = pair_bracket_rank[pop_other],
      annotations  = vapply(mean_p, fmt_p, character(1)),
      chrom        = factor(chrom, levels = levels(plot_df$chrom))
    )

  # 5. y positions
  if (target_sat %in% free_y_sats) {
    afr_annot <- afr_annot %>% left_join(chrom_max_df, by = "chrom") %>%
      mutate(y_position = y_data_max * 1.12 + (bracket_rank - 1L) * y_data_max * 0.10)
    n_per_pop <- n_per_pop %>% left_join(chrom_max_df, by = "chrom") %>%
      mutate(n_y = y_data_max * 1.03)
    bracket_df <- afr_annot %>% mutate(
      x1    = 3L,
      x2    = as.integer(factor(pop_other, levels = pop_order)),
      x_mid = (3L + as.integer(factor(pop_other, levels = pop_order))) / 2,
      y_lo  = y_position - y_data_max * 0.018
    )
    y_scale <- scale_y_continuous(limits = c(0, NA),
                                  expand = expansion(mult = c(0.02, 0.18)))
  } else {
    global_max     <- max(plot_df$len_mbp, na.rm = TRUE)
    y_bracket_base <- global_max * 1.12
    bracket_step   <- global_max * 0.09
    y_tick_h       <- global_max * 0.015
    afr_annot <- afr_annot %>%
      mutate(y_position = y_bracket_base + (bracket_rank - 1L) * bracket_step)
    y_lim_max <- y_bracket_base + 4L * bracket_step + bracket_step * 1.2
    n_per_pop <- n_per_pop %>% mutate(n_y = global_max * 1.03)
    bracket_df <- afr_annot %>% mutate(
      x1    = 3L,
      x2    = as.integer(factor(pop_other, levels = pop_order)),
      x_mid = (3L + as.integer(factor(pop_other, levels = pop_order))) / 2,
      y_lo  = y_position - y_tick_h
    )
    y_scale <- scale_y_continuous(limits = c(0, y_lim_max))
  }

  # 6. Facet
  if (is.null(facet_arg)) {
    facet_arg <- if (target_sat == "ASat") {
      facet_wrap(~chrom, nrow = 4, ncol = 6)
    } else if (target_sat %in% free_y_sats) {
      facet_wrap(~chrom, nrow = 1, scales = "free_y", drop = TRUE)
    } else {
      facet_wrap(~chrom, nrow = 1, drop = TRUE)
    }
  }

  # 7. Build plot
  p <- ggplot(plot_df,
              aes(x = x_pos, y = len_mbp, fill = superpopulation, group = x_pos)) +
    geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.5,
                 linewidth = 0.35, width = 0.7) +
    geom_text(
      data        = n_per_pop,
      aes(x = x_pos, y = n_y, label = n_pop),
      inherit.aes = FALSE,
      size = n_label_size, vjust = 0, color = "black",
      family = base_family
    ) +
    facet_arg +
    scale_fill_manual(values = superpopulation_colors, breaks = pop_order) +
    scale_x_continuous(
      breaks = seq_along(pop_order), labels = pop_order,
      expand = expansion(add = 0.6)
    ) +
    y_scale +
    labs(
      title = if (show_title) paste0(target_sat, " satellite length across superpopulations") else NULL,
      x     = NULL,
      y     = if (show_y_lab) "Length (Mbp)" else NULL,
      fill  = "Superpopulation"
    ) +
    theme_bw(base_size = 14, base_family = base_family) +
    theme(
      text              = element_text(family = base_family, color = "black"),
      axis.text.x       = element_text(angle = 0, hjust = 0.5, size = axis_text_x_size, color = "black", family = base_family),
      axis.text.y       = element_text(size = axis_text_y_size, color = "black", family = base_family),
      axis.title        = element_text(size = axis_title_size, color = "black", family = base_family),
      axis.ticks        = element_line(color = "black"),
      strip.text        = element_text(face = "bold", size = strip_text_size, color = "black", family = base_family),
      strip.background  = element_rect(fill = "grey92", color = NA),
      panel.grid        = element_blank(),
      legend.position   = if (show_legend) "bottom" else "none",
      legend.text       = element_text(size = 13, family = base_family),
      legend.title      = element_text(size = 14, family = base_family),
      plot.title        = element_text(size = 15, face = "bold", family = base_family)
    )

  if (nrow(bracket_df) > 0) {
    p <- p +
      geom_segment(data = bracket_df,
                   aes(x = x1, xend = x2, y = y_position, yend = y_position),
                   inherit.aes = FALSE, color = "black", linewidth = 0.3) +
      geom_segment(data = bracket_df,
                   aes(x = x1, xend = x1, y = y_lo, yend = y_position),
                   inherit.aes = FALSE, color = "black", linewidth = 0.3) +
      geom_segment(data = bracket_df,
                   aes(x = x2, xend = x2, y = y_lo, yend = y_position),
                   inherit.aes = FALSE, color = "black", linewidth = 0.3) +
      geom_text(data = bracket_df,
                aes(x = x_mid, y = y_position, label = annotations),
                inherit.aes = FALSE,
                size = p_label_size, vjust = -0.3,
                family = base_family)
  }
  p
}

# ── ASat — standalone PDF (6 cols × 4 rows, square aspect, no title) ────────
if (!is.null(target_dfs$ASat)) {
  p_asat <- make_plot("ASat", target_dfs$ASat, results_dfs$ASat,
                      show_title       = FALSE,
                      axis_text_y_size = 16,
                      axis_text_x_size = 14,
                      axis_title_size  = 17,
                      strip_text_size  = 16)
  out_pdf <- "boxplot_satlen_superpop_ASat.pdf"
  ggsave(out_pdf, plot = p_asat, width = 20, height = 20, device = cairo_pdf)
  message("PDF saved: ", out_pdf)
}

# ── Combined PDF: BSat / HSat1 / HSat2 / HSat3 ───────────────────────────────
# Max 7 panels per row; satellites with more chroms span multiple rows.
# Each row has a vertical satellite-name tag on the right.
get_chroms <- function(sat_name) {
  td <- target_dfs[[sat_name]]
  if (is.null(td)) return(character(0))
  keep <- td %>%
    group_by(chrom) %>%
    summarise(mean_len = mean(length, na.rm = TRUE), .groups = "drop") %>%
    filter(mean_len >= MEAN_FILTER_BP) %>%
    pull(chrom) %>% as.character()
  chrom_order[chrom_order %in% keep]
}

chroms_per_sat <- list(
  BSat  = get_chroms("BSat"),
  HSat1 = get_chroms("HSat1"),
  HSat2 = get_chroms("HSat2"),
  HSat3 = get_chroms("HSat3")
)

MAX_COLS  <- 7L
TAG_WIDTH <- 0.35
sat_labels <- c(BSat = "βSat",
                HSat1 = "HSat1", HSat2 = "HSat2", HSat3 = "HSat3")

split_chunks <- function(v, n) {
  if (length(v) == 0) return(list())
  unname(split(v, ceiling(seq_along(v) / n)))
}

build_sat_row <- function(sat_name, chroms_use, show_legend) {
  if (length(chroms_use) == 0) return(NULL)
  td <- target_dfs[[sat_name]] %>%
    filter(chrom %in% chroms_use) %>%
    mutate(chrom = factor(chrom, levels = chroms_use))
  rd <- results_dfs[[sat_name]] %>% filter(chrom %in% chroms_use)
  facet_arg <- if (sat_name %in% free_y_sats) {
    facet_wrap(~chrom, nrow = 1, scales = "free_y", drop = TRUE)
  } else {
    facet_wrap(~chrom, nrow = 1, drop = TRUE)
  }
  make_plot(sat_name, td, rd,
            facet_arg    = facet_arg,
            show_title   = FALSE,
            show_legend  = show_legend,
            n_label_size = 4.5, p_label_size = 5.2)
}

sat_tag_plot <- function(sat_name) {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = "grey92", color = NA) +
    annotate("text", x = 0.5, y = 0.5,
             label = sat_labels[[sat_name]],
             angle = 270, fontface = "bold", size = 7,
             family = base_family) +
    theme_void(base_family = base_family) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)
}

sats_order <- c("BSat", "HSat1", "HSat2", "HSat3")
row_specs <- list()
for (sat in sats_order) {
  chunks <- split_chunks(chroms_per_sat[[sat]], MAX_COLS)
  for (i in seq_along(chunks)) {
    row_specs[[length(row_specs) + 1]] <- list(
      sat    = sat,
      chroms = chunks[[i]],
      n_cols = length(chunks[[i]])
    )
  }
}

n_rows_total <- length(row_specs)
ref_n        <- MAX_COLS

build_row_with_tag <- function(spec, is_last_row) {
  main_p <- build_sat_row(spec$sat, spec$chroms, show_legend = is_last_row)
  if (is.null(main_p)) return(NULL)
  tag_p <- sat_tag_plot(spec$sat)
  if (spec$n_cols >= ref_n) {
    wrap_plots(main_p, tag_p,
               nrow   = 1,
               widths = c(ref_n, TAG_WIDTH))
  } else {
    wrap_plots(main_p, plot_spacer(), tag_p,
               nrow   = 1,
               widths = c(spec$n_cols, ref_n - spec$n_cols, TAG_WIDTH))
  }
}

padded_rows <- lapply(seq_along(row_specs), function(i) {
  build_row_with_tag(row_specs[[i]], is_last_row = (i == n_rows_total))
})

combined <- wrap_plots(padded_rows, ncol = 1)

out_combo <- "boxplot_satlen_superpop_BSat_HSat.pdf"
ggsave(out_combo, plot = combined, width = 22, height = 24, device = cairo_pdf)
message("PDF saved: ", out_combo)

message("\nAll done.")
