rm(list = ls())

library("tidyverse")

## load data ##
satlen  <- read.csv("data/all.sat.length.xls", header = TRUE, sep = "\t")
df_cen  <- read.csv("data/cent_chrom.txt", header = TRUE, sep = "\t")

## 1. HSat1 = HSat1A + HSat1B (aggregate before filtering) ##
hsat1_df <- satlen %>%
  filter(satellite %in% c("HSat1A", "HSat1B")) %>%
  group_by(sample_hap_chrom, sample_hap, sample, hap, chrom, project) %>%
  summarise(length = sum(length), .groups = "drop") %>%
  mutate(satellite = "HSat1")

satlen_df <- bind_rows(satlen, hsat1_df) %>%
  filter(!satellite %in% c("HSat1A", "HSat1B"))

## 2. Keep only complete centromeres (filterflag == 0) ##
satlen_complete <- satlen_df %>%
  left_join(
    df_cen %>% select(sample_hap_chrom, filterflag),
    by = "sample_hap_chrom",
    relationship = "many-to-many"
  ) %>%
  filter(filterflag == 0)

cat("Rows before filter:", nrow(satlen_df), "\n")
cat("Rows after complete filter:", nrow(satlen_complete), "\n")

## 3. Summary function: mean, sd, Q1, median, Q3 ##
summarize_sat <- function(df) {
  df %>%
    group_by(chrom, satellite, project) %>%
    summarise(
      n     = n(),
      mean  = mean(length, na.rm = TRUE),
      sd    = sd(length, na.rm = TRUE),
      Q1    = quantile(length, 0.25, na.rm = TRUE),
      median = quantile(length, 0.50, na.rm = TRUE),
      Q3    = quantile(length, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

## 4a. Per-project summary ##
summary_by_project <- summarize_sat(satlen_complete)

## 4b. Global summary (all projects combined) ##
summary_global <- satlen_complete %>%
  group_by(chrom, satellite) %>%
  summarise(
    project = "Global",
    n     = n(),
    mean  = mean(length, na.rm = TRUE),
    sd    = sd(length, na.rm = TRUE),
    Q1    = quantile(length, 0.25, na.rm = TRUE),
    median = quantile(length, 0.50, na.rm = TRUE),
    Q3    = quantile(length, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

## 5. Combine and write output ##
final_summary <- bind_rows(summary_by_project, summary_global) %>%
  arrange(chrom, satellite, project)

head(final_summary, 20)
cat("Total rows in output:", nrow(final_summary), "\n")

write.table(final_summary, "satlen_summary_stats.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Output written to: satlen_summary_stats.tsv\n")
