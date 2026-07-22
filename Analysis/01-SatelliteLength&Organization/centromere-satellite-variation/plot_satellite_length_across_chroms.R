rm(list = ls())

library("tidyverse")
library("ggplot2")
library(ggnewscale)
##load data##
satlen  <- read.csv("data/all.sat.length.xls", header=TRUE, sep="\t")
df_cen <- read.csv("data/cent_chrom.txt", header=TRUE, sep="\t")
head(df_cen)
df <- satlen %>% left_join(df_cen %>% select(sample_hap_chrom, len, gapless, CL, filterflag), by = "sample_hap_chrom", relationship = "many-to-many")
head(df)

##filtered chromosomes##
filter_chrs = c(
  'C020-CHA-E20#Mat#chr6',
  'C045-CHA-N05#Mat#chr17',
  'C019-CHA-E19#Pat#chr9',
  'C076-CHA-NE16#Pat#chr14',
  'HG02666_hap2_chr15',
  'HG01114_hap1_chr16',
  'HG02769_hap2_chr20',
  'HG03452_hap1_chr4',
  'NA19036_hap2_chr14',
  'NA19434_hap2_chr20',
  'NA20847_hap2_chr17'
)
filtereddf <- df %>% filter(filterflag == 0) %>% filter( ! sample_hap_chrom %in% filter_chrs) %>% filter(satellite != "rDNA")
head(filtereddf)

##chrom order##
chrom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 'chrX', 'chrY')
filtereddf$chrom <- factor(filtereddf$chrom, levels = chrom_order)


##chrom color##
chrom_colors <- c(
  "chr1" = "#9D5427", "chr2" = "#D67C1B", "chr3" = "#C6A15B", "chr4" = "#F2E86D",
  "chr5" = "#FFC000", "chr6" = "#FF6533", "chr7" = "#FF4B3E", "chr8" = "#CB3C23",
  "chr9" = "#972D07", "chr10" = "#a4c558", "chr11" = "#59b669", "chr12" = "#3c8340",
  "chr13" = "#6DAB30", "chr14" = "#285943", "chr15" = "#7ec0b4", "chr16" = "#85A7CC",
  "chr17" = "#3f66a1", "chr18" = "#3897c5", "chr19" = "#afa7d8", "chr20" = "#a874b5",
  "chr21" = "#893f8b", "chr22" = "#b13f73", "chrX" = "#574D68", "chrY" = "#C52B94"
)

##satellite color##
satellite_colors <- c(
  "ASat" = "#891640", "HSat1A" = "#18988B", "HSat1B" = "#18988B", "HSat1" = "#18988B", "HSat2" = "#323366",
  "HSat3" = "#77A0BA", "BSat" = "#DCAED0", "GSat" = "#758660"
)

head(filtereddf)

##add HSat1 length##
hsat1_df <- filtereddf %>%
  filter(satellite %in% c("HSat1A", "HSat1B")) %>%
  group_by(sample, sample_hap_chrom, sample_hap, hap, chrom, project, gapless, CL, filterflag) %>%
  summarize(length = sum(length)) %>%
  mutate(satellite = "HSat1")
head(hsat1_df)

##remove HSat1A and HSat1B
filtereddf <- bind_rows(filtereddf, hsat1_df) %>% filter(satellite != "HSat1A" & satellite != "HSat1B")
apg_filtereddf <- filtereddf %>% filter(project == "APG")
head(filtereddf)

##satellite order##
satellite_order <- c("ASat", "HSat1", "HSat2", "HSat3", "BSat", "GSat")
filtereddf$satellite <- factor(filtereddf$satellite, levels = satellite_order)
apg_filtereddf$satellite <- factor(apg_filtereddf$satellite, levels = satellite_order)


##references
ref_sat <- filtereddf %>% filter(project == "Ref") %>%
mutate(sample_hap = factor(sample_hap, levels = c("YAO_Mat", "YAO_Pat", "HG002_Mat", "HG002_Pat", "CN1_Mat", "CN1_Pat", "CHM13_-")),
color = case_when(
      sample_hap == "CHM13_-" ~ "#008C8C",
      sample_hap == "CN1_Mat" ~ "#fdb462",
      sample_hap == "CN1_Pat" ~ "#fdb462",
      sample_hap == "HG002_Mat" ~ "#b3de69",
      sample_hap == "HG002_Pat" ~ "#b3de69",
      sample_hap == "YAO_Mat" ~ "#bc80bd",
      sample_hap == "YAO_Pat" ~ "#bc80bd"
    ),
    shape = case_when(
      sample_hap == "CHM13_-" ~ 23,
      sample_hap == "CN1_Mat" ~ 23,
      sample_hap == "CN1_Pat" ~ 22,
      sample_hap == "HG002_Mat" ~ 23,
      sample_hap == "HG002_Pat" ~ 22,
      sample_hap == "YAO_Mat" ~ 23,
      sample_hap == "YAO_Pat" ~ 22,
    )) %>% arrange(sample_hap)
head(ref_sat)

##output stat file##
apg_stat <- apg_filtereddf %>%
  group_by(chrom, satellite) %>%
  summarise(mean_length = mean(length),
            sd_length = sd(length),
            min_len = min(length),
            max_length = max(length),
            median_length = median(length),
            q1            = quantile(length, 0.25),
            q3            = quantile(length, 0.75),
            .groups = "drop")
print(apg_stat)
write.table(apg_stat, file = "apg_satellite_length_stat.xls", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


##violin plot 
violin_plot <- ggplot(apg_filtereddf, aes(x = chrom, y = length / 1000000)) +
  geom_violin(aes(color = chrom), scale = "width", alpha = 0.8, draw_quantiles = c(0.25, 0.5, 0.75)) +
  facet_wrap(~ satellite, ncol = 1, scales = "free_y") +
  scale_color_manual(values = chrom_colors) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "plain", color = "black"),
    axis.text.x = element_text(size = 12, face = "plain", color = "black"),
    axis.text.y = element_text(size = 12, face = "plain", color = "black"),
    strip.text.x = element_text(colour = "black", face = "bold")
  ) +
  labs(x = "", y = "Length (Mb)") +
  new_scale_color() +
  geom_point(
    data = ref_sat,
    aes(x = chrom, y = length / 1000000, shape = factor(shape), fill = color, color = color),
    size = 1.5,
    position = position_nudge(x = 0.4)
  ) +
  scale_shape_manual(values = c(22, 23)) +
  scale_fill_identity() +
  scale_color_identity()

print(violin_plot)
outpdf <- paste0("apg_satellite_length_variation.pdf")
ggsave(outpdf, plot = violin_plot, device = "pdf", width = 16, height = 16)


##outlier plot
p <- ggplot(apg_filtereddf, aes(x = chrom, y = length / 1000000, color = chrom)) +
  # Error bar: mean ± SD
  stat_summary(fun.data = function(y) {
    m <- mean(y, na.rm = TRUE)
    s <- sd(y, na.rm = TRUE)
    data.frame(y = m, ymin = m - s, ymax = m + s)
  }, geom = "errorbar", width = 0.1, size = 0.3) +
  # Mean point
  stat_summary(fun = mean, geom = "point", shape = 21, size = 4,
               fill = "white", stroke = 0.3) +
  # Outliers based on mean ± 2×SD
  stat_summary(fun.data = function(y) {
    m <- mean(y, na.rm = TRUE)
    s <- sd(y, na.rm = TRUE)
    threshold <- 2
    outliers <- y[abs(y - m) > threshold * s]
    data.frame(y = outliers)
  }, geom = "point", shape = 4, size = 2) +
  facet_wrap(~ satellite, ncol = 1, scales = "free_y") +
  scale_color_manual(values = chrom_colors) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "plain", color = "black"),
    axis.text.x = element_text(size = 12, face = "plain", color = "black"),
    axis.text.y = element_text(size = 12, face = "plain", color = "black"),
    strip.text.x = element_text(colour = "black", face = "bold")
  ) +
  labs(x = "", y = "Length (Mb)") +
  new_scale_color() +
  geom_point(
    data = ref_sat,
    aes(x = chrom, y = length / 1000000, shape = factor(shape), fill = color, color = color),
    size = 1.5,
    position = position_nudge(x = 0.4)
  ) +
  scale_shape_manual(values = c(22, 23)) +
  scale_fill_identity() +
  scale_color_identity()

print(p)
outpdf <- "apg_satellite_length_variation_errorbar.pdf"
ggsave(outpdf, plot = p, device = "pdf", width = 16, height = 16)

##outlier table##
outlier_table <- apg_filtereddf %>%
  group_by(satellite, chrom) %>%
  mutate(
    group_mean   = mean(length, na.rm = TRUE),
    group_sd     = sd(length, na.rm = TRUE),
    z_score      = (length - group_mean) / group_sd,
    is_outlier   = abs(length - group_mean) > 2 * group_sd
  ) %>%
  ungroup() %>%
  select(satellite, chrom, sample, hap, sample_hap, sample_hap_chrom,
         length, group_mean, group_sd, z_score, is_outlier) %>%
  arrange(satellite, chrom, desc(abs(z_score)))

write.table(outlier_table, file = "apg_satellite_outlier_detail.xls",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

outlier_only <- outlier_table %>% filter(is_outlier)
print(outlier_only)
write.table(outlier_only, file = "apg_satellite_outlier_only.xls",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##define reliable subsets##
unreliable_sats   <- c("HSat1", "HSat2", "HSat3", "BSat", "GSat")
unreliable_chroms <- c("chr13", "chr14", "chr15", "chr21", "chr22")

##keep only large satellite groups (group_mean > 500 kb), then remove unreliable sat×chrom pairs##
large_outlier_table <- outlier_table %>%
  filter(group_mean > 500000) %>%
  filter(!(chrom %in% unreliable_chroms & satellite %in% unreliable_sats))

large_outlier_only <- large_outlier_table %>% filter(is_outlier)

##per satellite summary##
per_sat_summary <- large_outlier_table %>%
  group_by(satellite) %>%
  summarise(
    total_obs      = n(),
    n_outlier      = sum(is_outlier),
    outlier_rate   = round(n_outlier / total_obs, 4),
    n_high_outlier = sum(is_outlier & z_score > 0),
    n_low_outlier  = sum(is_outlier & z_score < 0),
    .groups = "drop"
  ) %>%
  arrange(desc(outlier_rate))
print(per_sat_summary)
write.table(per_sat_summary, file = "outlier_summary_per_satellite.xls",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##per sample summary##
per_sample_summary <- large_outlier_only %>%
  group_by(sample_hap) %>%
  summarise(
    n_outlier       = n(),
    n_high_outlier  = sum(z_score > 0),
    n_low_outlier   = sum(z_score < 0),
    chroms_affected = paste(sort(unique(as.character(chrom))), collapse = ","),
    sats_affected   = paste(sort(unique(as.character(satellite))), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(desc(n_outlier))
print(per_sample_summary)
write.table(per_sample_summary, file = "outlier_summary_per_sample.xls",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##per chrom summary##
per_chrom_summary <- large_outlier_table %>%
  filter(!chrom %in% unreliable_chroms) %>%
  group_by(chrom) %>%
  summarise(
    total_obs      = n(),
    n_outlier      = sum(is_outlier),
    outlier_rate   = round(n_outlier / total_obs, 4),
    n_high_outlier = sum(is_outlier & z_score > 0),
    n_low_outlier  = sum(is_outlier & z_score < 0),
    .groups = "drop"
  ) %>%
  mutate(chrom = factor(chrom, levels = chrom_order)) %>%
  arrange(chrom)
print(per_chrom_summary)
write.table(per_chrom_summary, file = "outlier_summary_per_chrom.xls",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
