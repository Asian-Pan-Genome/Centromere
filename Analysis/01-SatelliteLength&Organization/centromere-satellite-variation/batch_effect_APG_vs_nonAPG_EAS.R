rm(list = ls())

library(tidyverse)
library(ggplot2)
library(patchwork)

# NOTE: do NOT use showtext here.
# showtext::showtext_auto() draws every glyph as vector outlines (paths),
# which makes the text non-editable in Adobe Illustrator (this was the v5 bug).
# Instead we let cairo_pdf embed the real Arial font, so text stays editable.
# Arial is a system font on Windows; cairo_pdf resolves "Arial" by family name.

base_family <- "Arial"

# ── Load data ──────────────────────────────────────────────────────────────────
satlen <- read.csv("data/all.sat.length.xls", header = TRUE, sep = "\t")
df_cen <- read.csv("data/cent_chrom.txt",     header = TRUE, sep = "\t")
popdf  <- read.csv("data/populaion.xls",      header = TRUE, sep = "\t")

# ── Preprocess ─────────────────────────────────────────────────────────────────
popdf <- popdf %>%
  mutate(superpopulation = ifelse(superpopulation == "EAS-APG", "EAS", superpopulation))

eas_samples <- popdf %>% filter(superpopulation == "EAS") %>% pull(sample)

chrom_levels <- paste0("chr", c(1:22, "X", "Y"))

# Join filterflag, keep EAS + relevant projects, merge HSat1A/B → HSat1
sat_eas <- satlen %>%
  left_join(
    df_cen %>% select(sample_hap_chrom, filterflag),
    by = "sample_hap_chrom",
    relationship = "many-to-many"
  ) %>%
  filter(filterflag == 0,
         sample %in% eas_samples,
         project %in% c("APG", "HPRC", "HGSVC")) %>%
  mutate(satellite = ifelse(grepl("^HSat1", satellite), "HSat1", satellite)) %>%
  group_by(sample_hap_chrom, sample_hap, sample, hap, chrom, project, satellite) %>%
  summarise(len_mbp = sum(length) / 1e6, .groups = "drop") %>%
  mutate(
    group = factor(ifelse(project == "APG", "APG", "non-APG"), levels = c("APG", "non-APG")),
    chrom = factor(chrom, levels = chrom_levels)
  ) %>%
  filter(satellite %in% c("ASat", "BSat", "HSat1", "HSat2", "HSat3"))

# ── Wilcoxon test with downsampling ───────────────────────────────────────────
# Per satellite × chromosome:
#   n_APG > 10 and n_APG > n_other → downsample APG to n_other, 10 iterations, mean p-value
#   n_APG > 10 and n_APG ≤ n_other → single Wilcoxon test
#   n_APG ≤ 10 or n_other < 3      → skip
set.seed(42)
N_ITER <- 10

wilcox_list <- list()

for (sat in unique(sat_eas$satellite)) {
  for (chr in chrom_levels) {

    sub   <- sat_eas %>% filter(satellite == sat, chrom == chr)
    apg   <- sub %>% filter(group == "APG")
    other <- sub %>% filter(group == "non-APG")

    n_apg   <- n_distinct(apg$sample_hap_chrom)
    n_other <- n_distinct(other$sample_hap_chrom)

    if (n_apg <= 10 || n_other < 3) next

    if (n_apg > n_other) {
      p_vals <- sapply(seq_len(N_ITER), function(i) {
        ids     <- sample(unique(apg$sample_hap_chrom), n_other, replace = FALSE)
        apg_sub <- apg %>% filter(sample_hap_chrom %in% ids)
        tryCatch(
          wilcox.test(apg_sub$len_mbp, other$len_mbp, exact = FALSE)$p.value,
          error = function(e) NA_real_
        )
      })
      mean_p <- mean(p_vals, na.rm = TRUE)
    } else {
      mean_p <- tryCatch(
        wilcox.test(apg$len_mbp, other$len_mbp, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
    }

    wilcox_list[[paste(sat, chr, sep = "_")]] <- data.frame(
      satellite = sat,
      chrom     = chr,
      n_apg     = n_apg,
      n_other   = n_other,
      mean_p    = mean_p,
      stringsAsFactors = FALSE
    )
  }
}

wilcox_df <- bind_rows(wilcox_list) %>%
  mutate(
    chrom   = factor(chrom, levels = chrom_levels),
    p_label = case_when(
      is.na(mean_p)  ~ "",
      mean_p < 0.001 ~ paste0("p=", formatC(mean_p, format = "e", digits = 1)),
      TRUE           ~ paste0("p=", round(mean_p, 3))
    )
  )

write.table(wilcox_df, "batch_effect_APG_vs_nonAPG_EAS_wilcox.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Wilcoxon results saved.")

# ── Eligibility filter (per satellite × chrom) ────────────────────────────────
# Keep only combos with:
#   1. Downsampling Wilcoxon was performed (already in wilcox_df with non-NA mean_p)
#   2. Mean satellite length ≥ 500 kbp (= 0.5 Mbp)
#   3. NOT in chrY for BSat / HSat2 / HSat3
mean_len_df <- sat_eas %>%
  group_by(satellite, chrom) %>%
  summarise(mean_len_mbp = mean(len_mbp, na.rm = TRUE), .groups = "drop")

elig_df <- wilcox_df %>%
  filter(!is.na(mean_p)) %>%
  inner_join(mean_len_df, by = c("satellite", "chrom")) %>%
  filter(
    mean_len_mbp >= 0.5,
    !(satellite %in% c("BSat", "HSat2", "HSat3") & chrom == "chrY")
  )

sat_order <- c("ASat", "BSat", "HSat1", "HSat2", "HSat3")

plot_df <- sat_eas %>%
  semi_join(elig_df, by = c("satellite", "chrom")) %>%
  mutate(
    satellite = factor(satellite, levels = sat_order),
    chrom     = factor(chrom,     levels = chrom_levels)
  )

# ── Per (sat, chrom) bracket y positions ──────────────────────────────────────
y_range_df <- plot_df %>%
  group_by(satellite, chrom) %>%
  summarise(
    y_max = max(len_mbp, na.rm = TRUE),
    y_rng = pmax(max(len_mbp, na.rm = TRUE) - min(len_mbp, na.rm = TRUE),
                 max(len_mbp, na.rm = TRUE) * 0.05),
    .groups = "drop"
  ) %>%
  mutate(
    y_bar      = y_max + y_rng * 0.22,
    y_tick_bot = y_bar - y_rng * 0.06,
    y_text     = y_bar + y_rng * 0.05   # text bottom sits above the bar
  )

annot_df <- elig_df %>%
  mutate(
    satellite = factor(satellite, levels = sat_order),
    chrom     = factor(chrom,     levels = chrom_levels),
    # plain text (NOT plotmath) so "p = xxx" stays one editable string in Illustrator
    p_label   = ifelse(
      mean_p < 0.001,
      paste0("p = ", formatC(mean_p, format = "e", digits = 1)),
      paste0("p = ", formatC(round(mean_p, 3), format = "f", digits = 3))
    )
  ) %>%
  inner_join(y_range_df, by = c("satellite", "chrom"))

n_label_df <- elig_df %>%
  mutate(
    satellite = factor(satellite, levels = sat_order),
    chrom     = factor(chrom,     levels = chrom_levels),
    n_label   = paste0("n = ", n_other)
  )

# ── Per-satellite plots + patchwork (no empty cells, uniform panel widths) ────
group_colors <- c("APG" = "#E85827", "non-APG" = "#285943")
sat_labels   <- c(ASat = "αSat", BSat = "βSat",
                  HSat1 = "HSat1", HSat2 = "HSat2", HSat3 = "HSat3")

# Eligible chroms per satellite, preserving chrom_levels order
elig_chroms <- function(sat_name) {
  chrs <- elig_df %>% filter(satellite == sat_name) %>% pull(chrom) %>% as.character()
  chrom_levels[chrom_levels %in% chrs]
}

asat_chroms  <- elig_chroms("ASat")
n_half       <- ceiling(length(asat_chroms) / 2)
asat_top     <- asat_chroms[1:n_half]
asat_bot     <- asat_chroms[(n_half + 1):length(asat_chroms)]

bsat_chroms  <- elig_chroms("BSat")
hsat1_chroms <- elig_chroms("HSat1")
hsat2_chroms <- elig_chroms("HSat2")
hsat3_chroms <- elig_chroms("HSat3")

# Reference column count = ASat half-row (= widest row)
ref_n <- max(n_half, length(bsat_chroms), length(hsat1_chroms),
             length(hsat2_chroms), length(hsat3_chroms))

# Subplot builder: one satellite × subset of chroms
build_subplot <- function(sat_name, chroms_use,
                          show_x = TRUE,
                          show_y_title = TRUE,
                          show_legend = FALSE) {
  pdf <- plot_df    %>% filter(satellite == sat_name, chrom %in% chroms_use) %>%
    mutate(chrom = factor(chrom, levels = chroms_use))
  adf <- annot_df   %>% filter(satellite == sat_name, chrom %in% chroms_use) %>%
    mutate(chrom = factor(chrom, levels = chroms_use))
  ndf <- n_label_df %>% filter(satellite == sat_name, chrom %in% chroms_use) %>%
    mutate(chrom = factor(chrom, levels = chroms_use))

  p <- ggplot(pdf, aes(x = group, y = len_mbp, fill = group)) +
    geom_boxplot(outlier.size = 0.3, width = 0.6, linewidth = 0.4) +
    geom_segment(data = adf,
                 aes(x = 1, xend = 2, y = y_bar, yend = y_bar),
                 inherit.aes = FALSE, linewidth = 0.4, color = "black") +
    geom_segment(data = adf,
                 aes(x = 1, xend = 1, y = y_tick_bot, yend = y_bar),
                 inherit.aes = FALSE, linewidth = 0.4, color = "black") +
    geom_segment(data = adf,
                 aes(x = 2, xend = 2, y = y_tick_bot, yend = y_bar),
                 inherit.aes = FALSE, linewidth = 0.4, color = "black") +
    geom_text(data = adf,
              aes(x = 1.5, y = y_text, label = p_label),
              inherit.aes = FALSE,
              size = 3.0, vjust = 0, color = "black", family = base_family) +
    geom_text(data = ndf,
              aes(x = Inf, y = Inf, label = n_label),
              inherit.aes = FALSE,
              size = 3.0, hjust = 1.08, vjust = 1.3,
              color = "black", family = base_family) +
    facet_grid(satellite ~ chrom, scales = "free_y",
               labeller = labeller(satellite = sat_labels)) +
    scale_fill_manual(values = group_colors) +
    labs(x = NULL,
         y = if (show_y_title) "Length (Mbp)" else NULL,
         fill = NULL) +
    theme_bw(base_size = 9, base_family = base_family) +
    theme(
      text               = element_text(family = base_family, color = "black"),
      legend.position    = if (show_legend) "bottom" else "none",
      legend.text        = element_text(family = base_family, color = "black", size = 10),
      strip.background   = element_rect(fill = "grey92", color = NA),
      strip.text.x       = element_text(size = 10, family = base_family,
                                        color = "black", face = "bold"),
      strip.text.y.right = element_text(size = 11, family = base_family,
                                        color = "black", face = "bold", angle = -90),
      axis.text.x        = if (show_x) element_text(size = 9, family = base_family,
                                                    color = "black", angle = 30, hjust = 1)
                          else element_blank(),
      axis.ticks.x       = if (show_x) element_line(color = "black") else element_blank(),
      axis.text.y        = element_text(size = 9, family = base_family, color = "black"),
      axis.title.y       = element_text(size = 11, family = base_family, color = "black"),
      axis.ticks.y       = element_line(color = "black"),
      panel.grid         = element_blank(),
      panel.spacing.x    = unit(0.18, "lines"),
      panel.spacing.y    = unit(0.4, "lines")
    )
  p
}

# Build each row; x-axis labels only on the last row; legend only on last row
p_asat_top <- build_subplot("ASat",  asat_top,     show_x = FALSE)
p_asat_bot <- build_subplot("ASat",  asat_bot,     show_x = FALSE)
p_bsat     <- build_subplot("BSat",  bsat_chroms,  show_x = FALSE)
p_hsat1    <- build_subplot("HSat1", hsat1_chroms, show_x = FALSE)
p_hsat2    <- build_subplot("HSat2", hsat2_chroms, show_x = FALSE)
p_hsat3    <- build_subplot("HSat3", hsat3_chroms, show_x = TRUE,
                            show_legend = TRUE)

# Pad shorter rows with spacers so panel widths align
wrap_row <- function(p, n_cols, ref_n) {
  if (n_cols >= ref_n) {
    return(p)
  }
  wrap_plots(p, plot_spacer(),
             nrow = 1,
             widths = c(n_cols, ref_n - n_cols))
}

row_top  <- wrap_row(p_asat_top, length(asat_top),     ref_n)
row_bot  <- wrap_row(p_asat_bot, length(asat_bot),     ref_n)
row_bsat <- wrap_row(p_bsat,     length(bsat_chroms),  ref_n)
row_h1   <- wrap_row(p_hsat1,    length(hsat1_chroms), ref_n)
row_h2   <- wrap_row(p_hsat2,    length(hsat2_chroms), ref_n)
row_h3   <- wrap_row(p_hsat3,    length(hsat3_chroms), ref_n)

combined <- wrap_plots(row_top, row_bot, row_bsat, row_h1, row_h2, row_h3,
                       ncol = 1)

out_file <- "batch_effect_EAS_APG_vs_nonAPG_all_v6.pdf"
ggsave(out_file, plot = combined,
       width = 22, height = 14,
       device = cairo_pdf)
message("Saved: ", out_file)
