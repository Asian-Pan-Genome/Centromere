rm(list = ls())

library("tidyverse")

set.seed(42)

## ============================================================================
## Load data
## ============================================================================
satlen  <- read.csv("data/all.sat.length.xls", header = TRUE, sep = "\t")
df_cen  <- read.csv("data/cent_chrom.txt", header = TRUE, sep = "\t")
pop_df  <- read.csv("data/populaion.xls",
                    header = TRUE, sep = "\t") %>%
  select(sample, superpopulation) %>%
  mutate(superpopulation = if_else(superpopulation == "EAS-APG", "EAS", superpopulation))

cat("satlen rows:", nrow(satlen), "\n")
cat("cen rows:", nrow(df_cen), "\n")
cat("pop rows:", nrow(pop_df), "\n")

## ============================================================================
## 1. HSat1 = HSat1A + HSat1B
## ============================================================================
hsat1_df <- satlen %>%
  filter(satellite %in% c("HSat1A", "HSat1B")) %>%
  group_by(sample_hap_chrom, sample_hap, sample, hap, chrom, project) %>%
  summarise(length = sum(length, na.rm = TRUE), .groups = "drop") %>%
  mutate(satellite = "HSat1")

satlen_df <- bind_rows(
  satlen %>% filter(!satellite %in% c("HSat1A", "HSat1B")),
  hsat1_df
)

cat("After HSat1 merge, satlen rows:", nrow(satlen_df), "\n")
cat("Satellites:", unique(satlen_df$satellite), "\n")

## ============================================================================
## 2. Centromere data: dedup and pair haplotypes
##    Mat vs Pat (APG/HPRC/Ref), hap1 vs hap2 (HGSVC)
##    Both filterflag == 0
## ============================================================================

## Dedup: group by key columns, take max len and max filterflag
cen_dedup <- df_cen %>%
  group_by(sample_hap_chrom, sample, hap, chrom) %>%
  summarise(
    len        = first(len, order_by = desc(len)),
    filterflag = max(filterflag, na.rm = TRUE),
    .groups = "drop"
  )

## --- Trio-phased: Mat vs Pat ---
cen_trio_pair <- cen_dedup %>%
  filter(hap %in% c("Mat", "Pat")) %>%
  select(sample, chrom, hap, len, filterflag) %>%
  pivot_wider(names_from = hap, values_from = c(len, filterflag)) %>%
  filter(!is.na(len_Mat), !is.na(len_Pat),
         filterflag_Mat == 0, filterflag_Pat == 0) %>%
  mutate(
    Maternal_censize   = len_Mat,
    Paternal_censize   = len_Pat,
    censize_difference = len_Mat - len_Pat,
    haplotype_pair     = "Mat-Pat"
  )

cat("Trio-phased centromere pairs (Mat-Pat):", nrow(cen_trio_pair), "\n")

## --- HGSVC: hap1 vs hap2 ---
cen_hgsvc_pair <- cen_dedup %>%
  filter(hap %in% c("hap1", "hap2")) %>%
  select(sample, chrom, hap, len, filterflag) %>%
  pivot_wider(names_from = hap, values_from = c(len, filterflag)) %>%
  filter(!is.na(len_hap1), !is.na(len_hap2),
         filterflag_hap1 == 0, filterflag_hap2 == 0) %>%
  mutate(
    Maternal_censize   = len_hap1,
    Paternal_censize   = len_hap2,
    censize_difference = len_hap1 - len_hap2,
    haplotype_pair     = "Hap1-Hap2"
  )

cat("HGSVC centromere pairs (hap1-hap2):", nrow(cen_hgsvc_pair), "\n")

## Combine
cen_pairs <- bind_rows(cen_trio_pair, cen_hgsvc_pair) %>%
  select(sample, chrom, Maternal_censize, Paternal_censize,
         censize_difference, haplotype_pair)

cat("Total centromere pairs:", nrow(cen_pairs), "\n")

## ============================================================================
## 3. Satellite pairing: same logic, using filterflag from centromere
## ============================================================================

## Attach filterflag from cen_dedup
sat_with_flag <- satlen_df %>%
  left_join(
    cen_dedup %>% select(sample_hap_chrom, filterflag),
    by = "sample_hap_chrom"
  )

## Sum satellite lengths per (sample, chrom, satellite, hap) to handle
## multiple segments of the same satellite type on one chromosome
sat_sum <- sat_with_flag %>%
  group_by(sample, chrom, satellite, hap, filterflag) %>%
  summarise(length = sum(length, na.rm = TRUE), .groups = "drop")

cat("Sat sum rows:", nrow(sat_sum), "\n")

## --- Trio-phased satellites (Mat vs Pat) ---
sat_trio_pair <- sat_sum %>%
  filter(hap %in% c("Mat", "Pat")) %>%
  select(sample, chrom, satellite, hap, length, filterflag) %>%
  pivot_wider(names_from = hap, values_from = c(length, filterflag)) %>%
  filter(!is.na(length_Mat), !is.na(length_Pat),
         filterflag_Mat == 0, filterflag_Pat == 0) %>%
  mutate(sat_difference = length_Mat - length_Pat)

cat("Trio-phased satellite pairs:", nrow(sat_trio_pair), "\n")

## --- HGSVC satellites (hap1 vs hap2) ---
sat_hgsvc_pair <- sat_sum %>%
  filter(hap %in% c("hap1", "hap2")) %>%
  select(sample, chrom, satellite, hap, length, filterflag) %>%
  pivot_wider(names_from = hap, values_from = c(length, filterflag)) %>%
  filter(!is.na(length_hap1), !is.na(length_hap2),
         filterflag_hap1 == 0, filterflag_hap2 == 0) %>%
  mutate(sat_difference = length_hap1 - length_hap2)

cat("HGSVC satellite pairs:", nrow(sat_hgsvc_pair), "\n")

sat_pairs <- bind_rows(sat_trio_pair, sat_hgsvc_pair) %>%
  select(sample, chrom, satellite, sat_difference)

cat("Total satellite pairs:", nrow(sat_pairs), "\n")

## ============================================================================
## 4. Pivot satellite differences to wide (one column per satellite)
## ============================================================================
sat_wide <- sat_pairs %>%
  pivot_wider(
    names_from  = satellite,
    values_from = sat_difference,
    names_glue  = "{satellite}_difference"
  )

## Ensure expected satellite columns exist
expected_sats <- c("ASat_difference", "BSat_difference", "GSat_difference",
                   "HSat1_difference", "HSat2_difference", "HSat3_difference")
for (s in expected_sats) {
  if (!(s %in% names(sat_wide))) {
    sat_wide[[s]] <- NA_real_
  }
}

cat("Sat wide rows:", nrow(sat_wide), "\n")

## ============================================================================
## 5. Join centromere pairs + satellite wide + project + superpopulation
## ============================================================================

## Get project per sample
sample_project <- df_cen %>%
  select(sample, project) %>%
  distinct()

## Join
final_df <- cen_pairs %>%
  left_join(sat_wide, by = c("sample", "chrom")) %>%
  left_join(sample_project, by = "sample") %>%
  left_join(pop_df, by = "sample")

cat("Final joined rows:", nrow(final_df), "\n")
cat("Samples with missing superpop:\n")
final_df %>% filter(is.na(superpopulation)) %>%
  distinct(sample, project) %>% print()

## ============================================================================
## 6. Outlier detection: per chromosome, mean +/- 2*sd of censize_difference
## ============================================================================
outlier_stats <- final_df %>%
  group_by(chrom) %>%
  summarise(
    chrom_mean = mean(censize_difference, na.rm = TRUE),
    chrom_sd   = sd(censize_difference, na.rm = TRUE),
    lower_2sd  = chrom_mean - 2 * chrom_sd,
    upper_2sd  = chrom_mean + 2 * chrom_sd,
    .groups = "drop"
  )

final_df <- final_df %>%
  left_join(outlier_stats, by = "chrom") %>%
  mutate(
    Whether_is_outlier = if_else(
      censize_difference < lower_2sd | censize_difference > upper_2sd,
      TRUE, FALSE
    )
  )

cat("\nOutlier counts:\n")
print(final_df %>% count(Whether_is_outlier))
cat("\nOutliers per chrom:\n")
print(final_df %>% filter(Whether_is_outlier) %>% count(chrom) %>% arrange(desc(n)))

## ============================================================================
## 7. Format final output
## ============================================================================
output_df <- final_df %>%
  mutate(
    Sample_hap_chromosome = paste0(sample, "#", chrom)
  ) %>%
  select(
    Sample_hap_chromosome,
    sample,
    chrom,
    haplotype = haplotype_pair,
    project,
    superpopulation,
    `Maternal(haplotype1)_censize` = Maternal_censize,
    `Paternal(haplotype2)_censize` = Paternal_censize,
    Difference                     = censize_difference,
    chrom_mean,
    chrom_sd,
    lower_2sd,
    upper_2sd,
    Whether_is_outlier,
    ASat_difference,
    BSat_difference,
    HSat1_difference,
    HSat2_difference,
    HSat3_difference,
    GSat_difference
  ) %>%
  arrange(project, chrom, sample)

cat("\nFinal output rows:", nrow(output_df), "\n")
cat("Columns:", ncol(output_df), "\n")
cat("\nPer project:\n")
print(output_df %>% count(project))
cat("\nPer haplotype pair type:\n")
print(output_df %>% count(haplotype))
cat("\nPer superpopulation:\n")
print(output_df %>% count(superpopulation))
cat("\nOutliers total:", sum(output_df$Whether_is_outlier, na.rm = TRUE), "\n")

## Preview
cat("\n--- First 15 rows ---\n")
print(head(output_df, 15), width = 200)

## ============================================================================
## 8. Write output
## ============================================================================
write.table(output_df,
            file = "allsamples_mat_vs_pat_comprehensive.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
            na = "", eol = "\n")

cat("\nOutput written to: allsamples_mat_vs_pat_comprehensive.tsv\n")
cat("Done.\n")
