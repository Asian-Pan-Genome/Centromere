rm(list = ls())

library("tidyverse")
library("ggplot2")
library("patchwork")
library("ggh4x")        # facet_grid2: per-panel (independent) free scales
library("systemfonts")  # Arial lookup for cairo_pdf (editable text in Illustrator)

##Arial is already a registered system font on this machine; cairo_pdf will##
##embed it as live (editable) text. Just sanity-check it is available.##
stopifnot("Arial" %in% systemfonts::system_fonts()$family)

##chrom order & colors##
chrom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

chrom_colors <- c(
  "chr1" = "#9D5427", "chr2" = "#D67C1B", "chr3" = "#C6A15B", "chr4" = "#F2E86D",
  "chr5" = "#FFC000", "chr6" = "#FF6533", "chr7" = "#FF4B3E", "chr8" = "#CB3C23",
  "chr9" = "#972D07", "chr10" = "#a4c558", "chr11" = "#59b669", "chr12" = "#3c8340",
  "chr13" = "#6DAB30", "chr14" = "#285943", "chr15" = "#7ec0b4", "chr16" = "#85A7CC",
  "chr17" = "#3f66a1", "chr18" = "#3897c5", "chr19" = "#afa7d8", "chr20" = "#a874b5",
  "chr21" = "#893f8b", "chr22" = "#b13f73", "chrX" = "#574D68", "chrY" = "#C52B94"
)

satellite_order <- c("ASat", "HSat1", "HSat2", "HSat3", "BSat", "GSat")

##display labels for satellites (Greek letters for ASat / BSat / GSat)##
satellite_display <- c(
  "ASat"  = "\u03B1Sat",
  "HSat1" = "HSat1",
  "HSat2" = "HSat2",
  "HSat3" = "HSat3",
  "BSat"  = "\u03B2Sat",
  "GSat"  = "\u03B3Sat"
)

##load data##
satlen  <- read.csv("data/all.sat.length.xls", header = TRUE, sep = "\t")
df_cen  <- read.csv("data/cent_chrom.txt", header = TRUE, sep = "\t")

##join filterflag from cent_chrom (collapse duplicate rows: a chrom is bad if any record has filterflag != 0)##
cen_flag <- df_cen %>%
  group_by(sample_hap_chrom) %>%
  summarise(filterflag = max(filterflag, na.rm = TRUE), .groups = "drop")

df <- satlen %>% left_join(cen_flag, by = "sample_hap_chrom")

##keep APG only, drop rDNA##
df <- df %>% filter(project == "APG", satellite != "rDNA")

##merge HSat1A + HSat1B -> HSat1##
hsat1_df <- df %>%
  filter(satellite %in% c("HSat1A", "HSat1B")) %>%
  group_by(sample_hap_chrom, sample_hap, sample, hap, chrom, project, filterflag) %>%
  summarise(length = sum(length), .groups = "drop") %>%
  mutate(satellite = "HSat1")

df <- bind_rows(df %>% filter(!satellite %in% c("HSat1A", "HSat1B")), hsat1_df)

##pivot to Mat/Pat per (sample, chrom, satellite); keep only pairs where both filterflag == 0##
pair_df <- df %>%
  select(sample, chrom, satellite, hap, length, filterflag) %>%
  pivot_wider(names_from = hap, values_from = c(length, filterflag)) %>%
  filter(!is.na(length_Mat), !is.na(length_Pat),
         filterflag_Mat == 0, filterflag_Pat == 0) %>%
  mutate(Mat_size  = length_Mat,
         Pat_size  = length_Pat,
         diff_size = length_Mat - length_Pat)

##compute per chrom x satellite mean length (using both haplotypes)##
group_mean_df <- df %>%
  filter(filterflag == 0) %>%
  group_by(chrom, satellite) %>%
  summarise(group_mean = mean(length, na.rm = TRUE), .groups = "drop")

##unreliable sat x chrom combinations to exclude##
unreliable_sats   <- c("HSat1", "HSat2", "HSat3", "BSat", "GSat")
unreliable_chroms <- c("chr13", "chr14", "chr15", "chr21", "chr22")

##restrict to combos with mean > 500 kb and drop unreliable combos##
valid_combos <- group_mean_df %>%
  filter(group_mean > 500000) %>%
  filter(!(chrom %in% unreliable_chroms & satellite %in% unreliable_sats)) %>%
  select(chrom, satellite)

########################################################################
## Cross-superpopulation comparison of (Mat - Pat) satellite diffs.    ##
## Only trio-phased (Mat/Pat) assemblies carry a parental sign:        ##
## APG + HPRC + Ref. HGSVC is hap1/hap2 (no parental phase) -> dropped.##
## APG is entirely EAS-APG (large n -> density); the other            ##
## superpopulations come from HPRC/Ref and are sparse -> shown as pts. ##
########################################################################

##sample -> superpopulation map##
pop_map <- read.csv("data/populaion.xls",
                    header = TRUE, sep = "\t") %>%
  select(sample, superpopulation) %>%
  ##recode EAS-APG -> EAS (Q1.1.4 convention): APG joins the EAS superpopulation##
  mutate(superpopulation = if_else(superpopulation == "EAS-APG", "EAS", superpopulation))

##all-project diploid pairs (HSat1A+B merged), both haplotypes filterflag 0.   ##
##Trio-phased assemblies (APG/HPRC/Ref) give Mat/Pat; HGSVC3 is hap1/hap2 with  ##
##no parental phase, so hap1-hap2 substitutes Mat-Pat (sign arbitrary).         ##
##Unified labels: h1 = Mat or hap1, h2 = Pat or hap2.                           ##
df_allproj <- satlen %>%
  left_join(cen_flag, by = "sample_hap_chrom") %>%
  filter(satellite != "rDNA", hap %in% c("Mat", "Pat", "hap1", "hap2")) %>%
  mutate(hp = if_else(hap %in% c("Mat", "hap1"), "h1", "h2")) %>%
  select(sample, hp, chrom, satellite, project, filterflag, length)

hsat1_all <- df_allproj %>%
  filter(satellite %in% c("HSat1A", "HSat1B")) %>%
  group_by(sample, hp, chrom, project, filterflag) %>%
  summarise(length = sum(length), .groups = "drop") %>%
  mutate(satellite = "HSat1")
df_allproj <- bind_rows(df_allproj %>% filter(!satellite %in% c("HSat1A", "HSat1B")),
                        hsat1_all)

pair_allproj <- df_allproj %>%
  pivot_wider(names_from = hp, values_from = c(length, filterflag)) %>%
  filter(!is.na(length_h1), !is.na(length_h2),
         filterflag_h1 == 0, filterflag_h2 == 0) %>%
  mutate(diff_size = length_h1 - length_h2,
         phasing   = if_else(project == "HGSVC", "hap1-hap2", "Mat-Pat")) %>%
  left_join(pop_map, by = "sample") %>%
  filter(!is.na(superpopulation)) %>%
  inner_join(valid_combos, by = c("chrom", "satellite"))

##(1) counts list: complete Mat/Pat pairs per superpopulation x chrom x satellite##
superpop_counts <- pair_allproj %>%
  group_by(superpopulation, chrom, satellite) %>%
  summarise(n_pairs = n(), .groups = "drop") %>%
  mutate(chrom = factor(chrom, levels = chrom_order),
         satellite = factor(satellite, levels = satellite_order)) %>%
  arrange(superpopulation, satellite, chrom)

write.table(superpop_counts,
            file = "apg_mat_vs_pat_satellite_superpop_pair_counts.xls",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##wide version (superpopulation columns) for quick reading##
superpop_counts_wide <- superpop_counts %>%
  pivot_wider(names_from = superpopulation, values_from = n_pairs, values_fill = 0) %>%
  arrange(satellite, chrom)
write.table(superpop_counts_wide,
            file = "apg_mat_vs_pat_satellite_superpop_pair_counts_wide.xls",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##per-superpopulation summary: complete pairs per (chrom x satellite) cell,    ##
##over the valid autosomal combos (absent cells counted as 0).                 ##
superpop_pair_summary <- superpop_counts_wide %>%
  filter(chrom %in% chrom_order[1:22]) %>%
  pivot_longer(-c(chrom, satellite),
               names_to = "superpopulation", values_to = "n_pairs") %>%
  group_by(superpopulation) %>%
  summarise(
    n_combos         = n(),
    combos_with_data = sum(n_pairs > 0),
    total_pairs      = sum(n_pairs),
    mean_pairs       = round(mean(n_pairs), 2),
    median_pairs     = median(n_pairs),
    min_pairs        = min(n_pairs),
    max_pairs        = max(n_pairs),
    .groups = "drop"
  ) %>%
  arrange(desc(median_pairs), desc(mean_pairs))
write.table(superpop_pair_summary,
            file = "apg_mat_vs_pat_satellite_superpop_pair_summary.xls",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("\nComplete diploid pairs (Mat-Pat, or hap1-hap2 for HGSVC) per individual:\n")
print(pair_allproj %>% distinct(superpopulation, phasing, sample) %>%
        count(superpopulation, phasing))
cat("\nPer-superpopulation pairs per (chrom x satellite) cell:\n")
print(superpop_pair_summary)

########################################################################
## All-project per-sample summary table (HPRC, HGSVC, Ref in addition  ##
## to APG). Columns match the APG summary table, with project and      ##
## phasing metadata added. HGSVC uses hap1-hap2 (no parental phase).   ##
########################################################################
allproj_summary <- pair_allproj %>%
  mutate(chrom = factor(chrom, levels = chrom_order),
         satellite = factor(satellite, levels = satellite_order),
         diff_mb = diff_size / 1e6) %>%
  select(project, superpopulation, phasing, sample, chrom, satellite,
         h1_size = length_h1, h2_size = length_h2, diff_size, diff_mb) %>%
  arrange(project, superpopulation, phasing, chrom, satellite, sample)

write.table(allproj_summary,
            file = "allproj_mat_vs_pat_satellite_length_summary.xls",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("\nAll-project summary table written:", nrow(allproj_summary), "rows\n")
cat("Projects:\n")
print(allproj_summary %>% count(project, superpopulation, phasing))

##HPRCy1 and HGSVC3 pairs with |diff| > 10 Mbp##
large_diff <- pair_allproj %>%
  filter(project %in% c("HPRC", "HGSVC"), abs(diff_size) > 10e6) %>%
  mutate(diff_mb = round(diff_size / 1e6, 2),
         chrom   = factor(chrom, levels = chrom_order)) %>%
  arrange(project, chrom, satellite, desc(abs(diff_mb))) %>%
  select(project, superpopulation, phasing, sample, chrom, satellite,
         h1_size = length_h1, h2_size = length_h2, diff_mb)

cat("\nHPRCy1 / HGSVC3 pairs with |Mat-Pat| (or |hap1-hap2|) > 10 Mbp:\n")
print(large_diff, n = 100)
cat("Total:", nrow(large_diff), "pairs\n")

##superpopulation palette & order (canonical, see Q1.1.4 scripts)##
superpop_order  <- c("AFR", "AMR", "EAS", "EUR", "SAS")
superpop_colors <- c(AFR = "#319b62", AMR = "#939393", EAS = "#d22e77",
                     EUR = "#0070c0", SAS = "#893f8b")

pair_df <- pair_df %>%
  inner_join(valid_combos, by = c("chrom", "satellite"))

##per chrom x satellite stats (mean, median, sd, mean +/- 2sd, outliers) ##
combo_stats <- pair_df %>%
  group_by(chrom, satellite) %>%
  summarise(
    n              = n(),
    mean_diff      = mean(diff_size, na.rm = TRUE),
    median_diff    = median(diff_size, na.rm = TRUE),
    sd_diff        = sd(diff_size, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    lower_2sd = mean_diff - 2 * sd_diff,
    upper_2sd = mean_diff + 2 * sd_diff
  )

##add outlier flag to per-sample table##
summary_table <- pair_df %>%
  left_join(combo_stats %>% select(chrom, satellite, mean_diff, sd_diff,
                                    lower_2sd, upper_2sd),
            by = c("chrom", "satellite")) %>%
  mutate(
    is_outlier = diff_size < lower_2sd | diff_size > upper_2sd
  ) %>%
  mutate(chrom = factor(chrom, levels = chrom_order),
         satellite = factor(satellite, levels = satellite_order)) %>%
  select(satellite, chrom, sample, Mat_size, Pat_size, diff_size,
         mean_diff, sd_diff, lower_2sd, upper_2sd, is_outlier) %>%
  arrange(satellite, chrom, sample)

##add outlier count per combo##
combo_stats <- combo_stats %>%
  left_join(summary_table %>%
              group_by(chrom, satellite) %>%
              summarise(n_outlier = sum(is_outlier, na.rm = TRUE),
                        .groups = "drop"),
            by = c("chrom", "satellite")) %>%
  mutate(pct_outlier = round(n_outlier / n, 4)) %>%
  mutate(chrom = factor(chrom, levels = chrom_order),
         satellite = factor(satellite, levels = satellite_order)) %>%
  arrange(satellite, chrom)

write.table(summary_table,
            file = "apg_mat_vs_pat_satellite_length_summary.xls",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(combo_stats,
            file = "apg_mat_vs_pat_satellite_combo_stats.xls",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("Total pairs kept:", nrow(summary_table), "\n")
cat("Per satellite counts:\n")
print(summary_table %>% count(satellite))

##test for Mat/Pat bias per chrom x satellite (one-sample t-test against 0)##
bias_test <- pair_df %>%
  group_by(chrom, satellite) %>%
  summarise(
    n         = n(),
    mean_diff = mean(diff_size, na.rm = TRUE),
    sd_diff   = sd(diff_size, na.rm = TRUE),
    t_stat    = ifelse(n > 2 & sd_diff > 0,
                       t.test(diff_size, mu = 0)$statistic, NA_real_),
    p_value   = ifelse(n > 2 & sd_diff > 0,
                       t.test(diff_size, mu = 0)$p.value, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    bias_direction = case_when(
      is.na(p_adj)            ~ "NA",
      p_adj >= 0.05           ~ "no_bias",
      mean_diff > 0           ~ "Mat_biased",
      mean_diff < 0           ~ "Pat_biased"
    )
  ) %>%
  mutate(chrom = factor(chrom, levels = chrom_order),
         satellite = factor(satellite, levels = satellite_order)) %>%
  arrange(satellite, chrom)

write.table(bias_test,
            file = "apg_mat_vs_pat_satellite_bias_test.xls",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("\nBias direction summary (BH-adjusted p < 0.05):\n")
print(bias_test %>% count(bias_direction))
cat("\nSignificantly biased combos:\n")
print(bias_test %>% filter(bias_direction %in% c("Mat_biased", "Pat_biased")) %>%
      select(satellite, chrom, n, mean_diff, p_adj, bias_direction))

#####################################################################
## Variance-share decomposition: which satellite drives the total  ##
## (Mat - Pat) centromere-size difference on each chromosome.       ##
##                                                                  ##
## Total diff per (sample, chrom): d_total = sum_s d_s (over the      ##
## reliable/valid combos only, so the total is decomposed into the    ##
## same reliable components used throughout; acrocentrics therefore    ##
## reduce to alphaSat = 100%, their only reliable array).             ##
## Because Var(d_total) = sum_s Cov(d_s, d_total), each satellite      ##
## contributes a share = Cov(d_s, d_total) / Var(d_total), and the     ##
## shares sum to 100%. This is the formal test of "which satellite     ##
## contributes to the parental size difference".                       ##
#####################################################################

##total (Mat - Pat) diff per sample x chrom, summed over valid combos##
total_diff <- pair_df %>%
  group_by(sample, chrom) %>%
  summarise(d_total = sum(diff_size, na.rm = TRUE), .groups = "drop")

##per chrom x satellite: covariance share of the total variance##
var_share <- pair_df %>%
  select(sample, chrom, satellite, diff_size) %>%
  inner_join(total_diff, by = c("sample", "chrom")) %>%
  group_by(chrom, satellite) %>%
  summarise(
    n         = n(),
    sd_ds_mb  = sd(diff_size, na.rm = TRUE) / 1e6,
    cov_ds    = cov(diff_size, d_total, use = "complete.obs"),
    var_total = var(d_total,  na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  group_by(chrom) %>%
  mutate(
    sd_total_mb = sqrt(first(var_total)) / 1e6,
    var_share   = cov_ds / var_total              # fraction of total variance
  ) %>%
  ungroup() %>%
  ##restrict reporting to the same valid combos used elsewhere##
  inner_join(valid_combos, by = c("chrom", "satellite")) %>%
  mutate(var_share_pct = round(100 * var_share, 1)) %>%
  mutate(chrom = factor(chrom, levels = chrom_order),
         satellite = factor(satellite, levels = satellite_order)) %>%
  arrange(chrom, desc(var_share)) %>%
  select(chrom, satellite, n, sd_ds_mb, sd_total_mb, var_share, var_share_pct)

write.table(var_share,
            file = "apg_mat_vs_pat_satellite_variance_share.xls",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("\nDominant satellite per chromosome (variance share of total Mat-Pat diff):\n")
print(var_share %>% group_by(chrom) %>% slice_max(var_share, n = 1) %>%
      ungroup() %>% arrange(desc(var_share)) %>%
      select(chrom, satellite, n, sd_ds_mb, sd_total_mb, var_share_pct), n = 22)

##stacked-bar supplementary figure: % of total (Mat - Pat) variance by satellite##
bar_df <- var_share %>%
  mutate(var_share_pct_clip = pmax(var_share_pct, 0),   # clip tiny negatives for plotting
         satellite = factor(satellite_display[as.character(satellite)],
                            levels = satellite_display[satellite_order]))

satellite_fill <- c(
  "αSat" = "#891640", "HSat1" = "#18988B", "HSat2" = "#323366",
  "HSat3" = "#77A0BA", "βSat" = "#DCAED0", "γSat" = "#937860"
)

p_share <- ggplot(bar_df,
                  aes(x = chrom, y = var_share_pct_clip, fill = satellite)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = satellite_fill, name = "Satellite") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL,
       y = "Share of total (Maternal - Paternal)\nsize variance (%)") +
  theme_bw(base_family = "Arial") +
  theme(
    text            = element_text(family = "Arial", color = "black"),
    axis.title      = element_text(size = 11, family = "Arial", color = "black"),
    axis.text.x     = element_text(size = 9, family = "Arial", color = "black",
                                   angle = 45, hjust = 1),
    axis.text.y     = element_text(size = 9, family = "Arial", color = "black"),
    legend.title    = element_text(size = 10, family = "Arial", color = "black"),
    legend.text     = element_text(size = 9,  family = "Arial", color = "black"),
    panel.grid      = element_blank(),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )

ggsave("apg_mat_vs_pat_satellite_variance_share.pdf",
       plot = p_share, device = cairo_pdf, width = 9, height = 4.5)
ggsave("apg_mat_vs_pat_satellite_variance_share.png",
       plot = p_share, device = "png", width = 9, height = 4.5, dpi = 300)

##density plot: one panel per chrom x satellite##
plot_df <- summary_table %>%
  mutate(diff_mb = diff_size / 1e6) %>%
  filter(chrom %in% chrom_order[1:22])

##per-facet annotation: BH-adjusted P (italic P) + n (italic n), 2 lines stacked##
anno_df <- bias_test %>%
  select(chrom, satellite, p_adj, n) %>%
  filter(chrom %in% chrom_order[1:22]) %>%
  mutate(
    label_p = sprintf("italic(P)[adj]==%.2f", p_adj),
    label_n = sprintf("italic(n)==%d", n)
  )

##helper: split into top (chr1-11) / bottom (chr12-22) for a 2-row layout##
chr_top <- chrom_order[1:11]
chr_bot <- chrom_order[12:22]

mk_density_plot <- function(chr_subset, show_x = FALSE) {
  to_disp <- function(d) d %>%
    mutate(chrom = factor(chrom, levels = chr_subset),
           satellite = factor(satellite_display[as.character(satellite)],
                              levels = satellite_display[satellite_order]))
  pd <- plot_df %>% filter(chrom %in% chr_subset) %>% to_disp()
  ad <- anno_df %>% filter(chrom %in% chr_subset) %>% to_disp()

  ggplot(pd, aes(x = diff_mb, fill = chrom, color = chrom)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_text(data = ad, aes(x = -Inf, y = Inf, label = label_p),
              hjust = -0.05, vjust = 1.3, size = 2.9,
              color = "black", parse = TRUE, inherit.aes = FALSE, family = "Arial") +
    geom_text(data = ad, aes(x = -Inf, y = Inf, label = label_n),
              hjust = -0.05, vjust = 3.0, size = 2.9,
              color = "black", parse = TRUE, inherit.aes = FALSE, family = "Arial") +
    ##facet_grid2: keep the satellite x chrom matrix but give every panel##
    ##its own free x-axis (independent = "x"), so e.g. chr16 alpha-Sat is##
    ##not squished by the wide HSat2 range sharing the column.##
    facet_grid2(satellite ~ chrom, scales = "free", independent = "x") +
    scale_fill_manual(values = chrom_colors) +
    scale_color_manual(values = chrom_colors) +
    labs(x = if (show_x) "Satellite length difference (Maternal - Paternal, Mbp)" else NULL,
         y = "Density") +
    theme_bw(base_family = "Arial") +
    theme(
      text             = element_text(family = "Arial"),
      legend.position  = "none",
      axis.title       = element_text(size = 10, color = "black", family = "Arial"),
      axis.text.x      = element_text(size = 7, color = "black", family = "Arial"),
      axis.text.y      = element_text(size = 7, color = "black", family = "Arial"),
      strip.text.x     = element_text(size = 8, color = "black", face = "bold", family = "Arial"),
      strip.text.y     = element_text(size = 8, color = "black", face = "bold", family = "Arial"),
      strip.background = element_rect(fill = "grey95", color = NA),
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.4)
    )
}

##independent x-axes make every cell (incl. empty ones) draw its own##
##bottom axis. Blank the bottom (x) axis grob of empty cells so the##
##sparse satellite rows are not cluttered with redundant 0-centred ticks.##
##(Left/y axes are shared per row and are left intact.)##
blank_empty_panel_axes <- function(p, chr_subset) {
  pd <- plot_df %>% filter(chrom %in% chr_subset) %>%
    mutate(chrom    = factor(chrom, levels = chr_subset),
           sat_disp = factor(satellite_display[as.character(satellite)],
                             levels = satellite_display[satellite_order])) %>%
    droplevels()
  rows    <- levels(pd$sat_disp)   # top -> bottom
  cols    <- levels(pd$chrom)      # left -> right
  present <- pd %>% distinct(sat_disp, chrom)
  g <- ggplot2::ggplotGrob(p)
  for (ri in seq_along(rows)) for (ci in seq_along(cols)) {
    has <- nrow(filter(present, sat_disp == rows[ri], chrom == cols[ci])) > 0
    if (!has) {
      idx <- which(g$layout$name == sprintf("axis-b-%d-%d", ci, ri))
      if (length(idx) == 1) g$grobs[[idx]] <- ggplot2::zeroGrob()
    }
  }
  g
}

p_top <- mk_density_plot(chr_top, show_x = FALSE)
p_bot <- mk_density_plot(chr_bot, show_x = TRUE)

g_top <- blank_empty_panel_axes(p_top, chr_top)
g_bot <- blank_empty_panel_axes(p_bot, chr_bot)

##stack the two blocks; heights proportional to number of satellite rows##
n_top <- length(unique(plot_df$satellite[plot_df$chrom %in% chr_top]))
n_bot <- length(unique(plot_df$satellite[plot_df$chrom %in% chr_bot]))
p_density <- wrap_elements(full = g_top) / wrap_elements(full = g_bot) +
  plot_layout(heights = c(n_top, n_bot))

ggsave("apg_mat_vs_pat_satellite_length_density.pdf",
       plot = p_density, device = cairo_pdf, width = 14, height = 12)
ggsave("apg_mat_vs_pat_satellite_length_density.png",
       plot = p_density, device = "png", width = 14, height = 12, dpi = 300)

#####################################################################
## Superpopulation comparison figure (separate file).               ##
## Per chrom x satellite: a density for each superpopulation with    ##
## > DENS_MIN pairs, individual circles otherwise; n annotated per   ##
## superpopulation. EAS = APG + HPRC + HGSVC (EAS-APG recoded to EAS;##
## hap1-hap2 substitutes Mat-Pat for HGSVC, sign arbitrary - state   ##
## in legend). Downsampling Wilcoxon (K reps) compares EAS vs each    ##
## other superpopulation with > DENS_MIN pairs; mean p annotated.    ##
#####################################################################

DENS_MIN <- 10L          # density if n > DENS_MIN, else circles
K_DOWN   <- 10L          # downsampling repeats

sp_diff <- pair_allproj %>%
  filter(chrom %in% chrom_order[1:22]) %>%
  mutate(diff_mb = diff_size / 1e6,
         superpopulation = factor(superpopulation, levels = superpop_order))

sp_n <- sp_diff %>% count(chrom, satellite, superpopulation, name = "n_pairs")

##downsampling Wilcoxon: EAS vs each other superpop with > DENS_MIN pairs##
set.seed(42)
others <- c("AFR", "AMR", "EUR", "SAS")
combos <- sp_diff %>% distinct(chrom, satellite)
dn_list <- list()
for (i in seq_len(nrow(combos))) {
  ch <- combos$chrom[i]; sa <- combos$satellite[i]
  sub <- sp_diff %>% filter(chrom == ch, satellite == sa)
  eas <- sub$diff_mb[sub$superpopulation == "EAS"]
  for (S in others) {
    sv <- sub$diff_mb[sub$superpopulation == S]
    if (length(sv) > DENS_MIN && length(eas) > DENS_MIN) {
      ps <- replicate(K_DOWN, suppressWarnings(
              wilcox.test(sample(eas, length(sv)), sv)$p.value))
      dn_list[[length(dn_list) + 1]] <- data.frame(
        chrom = ch, satellite = sa, comparison = S,
        n_eas = length(eas), n_other = length(sv),
        mean_p = mean(ps, na.rm = TRUE))
    }
  }
}
dn_results <- bind_rows(dn_list) %>%
  mutate(p_adj = p.adjust(mean_p, method = "BH")) %>%   # BH across all comparisons
  arrange(satellite, chrom, comparison)
write.table(dn_results,
            file = "apg_mat_vs_pat_satellite_superpop_downsample_wilcox.xls",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##split density-eligible (n>DENS_MIN) vs circle (n<=DENS_MIN) per superpop##
sp_diff <- sp_diff %>% left_join(sp_n, by = c("chrom", "satellite", "superpopulation"))
sp_dens <- sp_diff %>% filter(n_pairs >  DENS_MIN)
sp_pts  <- sp_diff %>% filter(n_pairs <= DENS_MIN)

##per-panel annotation: "SP n=NN [p=PP]" one colored line per superpop##
sp_anno <- sp_n %>%
  left_join(dn_results %>% select(chrom, satellite, comparison, mean_p),
            by = c("chrom", "satellite", "superpopulation" = "comparison")) %>%
  mutate(
    label = if_else(!is.na(mean_p),
                    sprintf("%s n=%d p=%.2f", superpopulation, n_pairs, mean_p),
                    sprintf("%s n=%d", superpopulation, n_pairs)),
    vpos  = match(as.character(superpopulation), superpop_order)
  )

to_disp_sp <- function(d, chr_subset) d %>%
  filter(chrom %in% chr_subset) %>%
  mutate(chrom = factor(chrom, levels = chr_subset),
         satellite = factor(satellite_display[as.character(satellite)],
                            levels = satellite_display[satellite_order]))

mk_superpop_plot <- function(chr_subset, show_x = FALSE, show_legend = FALSE) {
  dd <- to_disp_sp(sp_dens, chr_subset)
  pp <- to_disp_sp(sp_pts,  chr_subset)
  aa <- to_disp_sp(sp_anno, chr_subset)

  ggplot() +
    geom_density(data = dd,
                 aes(x = diff_mb, fill = superpopulation, color = superpopulation,
                     group = superpopulation),
                 alpha = 0.30, linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_point(data = pp,
               aes(x = diff_mb, y = 0, color = superpopulation),
               shape = 16, size = 1.4, alpha = 0.9) +
    geom_text(data = aa,
              aes(x = -Inf, y = Inf, label = label, color = superpopulation),
              hjust = -0.05, vjust = aa$vpos * 1.35, size = 2.2,
              family = "Arial", show.legend = FALSE) +
    facet_grid2(satellite ~ chrom, scales = "free", independent = "x") +
    scale_fill_manual(values = superpop_colors, name = "Superpopulation", drop = FALSE) +
    scale_color_manual(values = superpop_colors, name = "Superpopulation", drop = FALSE) +
    labs(x = if (show_x) "Satellite length difference (Maternal - Paternal / hap1 - hap2, Mbp)" else NULL,
         y = "Density") +
    theme_bw(base_family = "Arial") +
    theme(
      text             = element_text(family = "Arial"),
      legend.position  = if (show_legend) "bottom" else "none",
      legend.title     = element_text(size = 10, color = "black", family = "Arial"),
      legend.text      = element_text(size = 9,  color = "black", family = "Arial"),
      axis.title       = element_text(size = 10, color = "black", family = "Arial"),
      axis.text.x      = element_text(size = 7, color = "black", family = "Arial"),
      axis.text.y      = element_text(size = 7, color = "black", family = "Arial"),
      strip.text.x     = element_text(size = 8, color = "black", face = "bold", family = "Arial"),
      strip.text.y     = element_text(size = 8, color = "black", face = "bold", family = "Arial"),
      strip.background = element_rect(fill = "grey95", color = NA),
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.4)
    )
}

sp_top <- mk_superpop_plot(chr_top, show_x = FALSE, show_legend = FALSE)
sp_bot <- mk_superpop_plot(chr_bot, show_x = TRUE,  show_legend = TRUE)

gs_top <- blank_empty_panel_axes(sp_top, chr_top)
gs_bot <- blank_empty_panel_axes(sp_bot, chr_bot)

p_superpop <- wrap_elements(full = gs_top) / wrap_elements(full = gs_bot) +
  plot_layout(heights = c(n_top, n_bot))

ggsave("apg_mat_vs_pat_satellite_superpop_density.pdf",
       plot = p_superpop, device = cairo_pdf, width = 14, height = 13)
ggsave("apg_mat_vs_pat_satellite_superpop_density.png",
       plot = p_superpop, device = "png", width = 14, height = 13, dpi = 300)

cat("\nDownsampling Wilcoxon (EAS vs other, K =", K_DOWN, ") comparisons:",
    nrow(dn_results), "\n")
print(dn_results)
