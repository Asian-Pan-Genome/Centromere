setwd("/share/home/zhanglab/user/sunyanqing/vscode_scripts/HORmining/pan_core_type")
rm(list = ls())

library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("patchwork")

###chromosome distribution of pan core type
hicat_pancore_chrom_df <-read.csv("20250911_hicat_hor_pan_core_type_across_chromosome_stat.xls", sep="\t", header=TRUE)
hicat_pancore_chrom_df <- read.csv("20250911_hormon_hor_pan_core_type_across_chromosome_stat.xls", sep="\t", header=TRUE)
head(hicat_pancore_chrom_df,10)

hicat_pancore_chrom_df <- hicat_pancore_chrom_df %>%
  mutate(pan_hor_type = case_when(
    pan_hor_type == "Low" ~ "Low frequency",
    pan_hor_type == "Low_frequency" ~ "Low frequency",
    TRUE ~ pan_hor_type
  )) %>%
  mutate(chrom_manual_check = factor(chrom_manual_check, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                                        "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                                        "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
                                        "chr22", "chrX", "chrY"))) %>%
  mutate(pan_hor_type = factor(pan_hor_type, levels=rev(c("Core", "Dispensable", "Low frequency", "Rare", "Singleton")))) %>%
  arrange(chrom_manual_check, pan_hor_type)


head(hicat_pancore_chrom_df,10)

# pantype_color <- c(
#     "Core" = "#67001f",
#     "Dispensable" = "#b2182b",
#     "Low frequency" = "#d6604d",
#     "Rare" = "#f4a582",
#     "Singleton" = "#fddbc7"
# )

pantype_color <- c(
    "Core" = "#313695",
    "Dispensable" = "#4575b4",
    "Low frequency" = "#74add1",
    "Rare" = "#abd9e9",
    "Singleton" = "#e0f3f8"
)

# Create the stacked bar plot
p <- ggplot(hicat_pancore_chrom_df, aes(x = chrom_manual_check, y = count, fill = pan_hor_type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = pantype_color) +
    labs(x = "Chromosome", y = "Number of HOR StVs", fill = "Pan HOR Type") +
    theme_bw() +
    theme(axis.text.x = element_text())
print(p)
ggsave(p, filename="20250911_hicat_hor_pan_core_type_across_chromosome_stat.pdf", width=10, height=6)
ggsave(p, filename="20250911_hormon_hor_pan_core_type_across_chromosome_stat.pdf", width=10, height=6)


###repeat number of pan core type
repeatnum_df <- read.csv("20250911_hicat_hor_pan_core_type_repeatnum_across_stvs_stat.xls", sep="\t", header=TRUE)
repeatnum_df <- read.csv("20250911_hormon_hor_pan_core_type_repeatnum_across_stvs_stat.xls", sep="\t", header=TRUE)
head(repeatnum_df,10)

library(ggridges)

# Define the custom colors
# values <- c(
#     "Core" = "#7F7340",
#     "Dispensable" = "#C6B28D",
#     "Low frequency" = "#BB8642",
#     "Rare" = "#914532",
#     "Singleton" = "#511F20"
# )
pantype_color <- c(
    "Core" = "#313695",
    "Dispensable" = "#4575b4",
    "Low frequency" = "#74add1",
    "Rare" = "#abd9e9",
    "Singleton" = "#e0f3f8"
)


# Create the ridge plot with log10 transformation and customized settings
p2 <- ggplot(repeatnum_df, aes(x = log10(mean), y = pan_hor_type, fill = pan_hor_type)) +
    geom_density_ridges(
        scale = 1, 
        rel_min_height = 0.01, 
        quantile_lines = TRUE, 
        quantiles = c(0.25, 0.5, 0.75),
        alpha = 0.8
    ) +
    scale_fill_manual(values = pantype_color) +
    scale_x_continuous(
        breaks = log10(c(1, 10, 100, 1000, 10000)), 
        labels = c("1", "10", "100", "1000", "10000")
    ) +
    labs(x = "Average copy number of HOR StVs", y = "", fill = "Pan HOR Type") +
    theme_classic()
print(p2)
ggsave(p2, filename="20250911_hicat_hor_pan_core_type_mean_repeatnum_across_stvs_stat.pdf", width=8, height=4)
ggsave(p2, filename="20250911_hormon_hor_pan_core_type_mean_repeatnum_across_stvs_stat.pdf", width=8, height=4)

p3 <- ggplot(repeatnum_df, aes(x = log10(max), y = pan_hor_type, fill = pan_hor_type)) +
    geom_density_ridges(
        scale = 1, 
        rel_min_height = 0.01, 
        quantile_lines = TRUE, 
        quantiles = c(0.25, 0.5, 0.75),
        alpha = 0.8
    ) +
    scale_fill_manual(values = values) +
    scale_x_continuous(
        breaks = log10(c(1, 10, 100, 1000, 10000)), 
        labels = c("1", "10", "100", "1000", "10000")
    ) +
    labs(x = "Average copy number of HOR StVs", y = "", fill = "Pan HOR Type") +
    theme_bw()
print(p3)
ggsave(p3, filename="20250911_hicat_hor_pan_core_type_max_repeatnum_across_stvs_stat.pdf", width=8, height=4)
ggsave(p3, filename="20250911_hormon_hor_pan_core_type_max_repeatnum_across_stvs_stat.pdf", width=8, height=4)


###significant test
pairwise_tests <- list(
    "Core vs Dispensable" = wilcox.test(mean ~ pan_hor_type, data = repeatnum_df %>% filter(pan_hor_type %in% c("Core", "Dispensable"))),
    "Core vs Rare" = wilcox.test(mean ~ pan_hor_type, data = repeatnum_df %>% filter(pan_hor_type %in% c("Core", "Rare"))),
    "Core vs Low frequency" = wilcox.test(mean ~ pan_hor_type, data = repeatnum_df %>% filter(pan_hor_type %in% c("Core", "Low frequency"))),
    "Core vs Singleton" = wilcox.test(mean ~ pan_hor_type, data = repeatnum_df %>% filter(pan_hor_type %in% c("Core", "Singleton")))
)

# Print the results
for (test_name in names(pairwise_tests)) {
    cat("\n", test_name, "\n")
    print(pairwise_tests[[test_name]])
}
