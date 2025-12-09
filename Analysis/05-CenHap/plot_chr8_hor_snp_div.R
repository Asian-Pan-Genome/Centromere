setwd("/share/home/zhanglab/user/sunyanqing/human/alpha/part2/chr8_mnphy/shared_sequence")
rm(list = ls())

library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("patchwork")
library("ggpubr") 

df <- read.csv("all_compare.div", sep="\t", header=FALSE)
colnames(df) <- c("type", "seq1", "seq2", "value")
head(df)


p <- ggplot(df, aes(x = type, y = value)) +
    geom_boxplot(fill = "grey", color = "black") +  # Hide outliers
    # geom_violin() +
    labs(x = "", y = "SqVs", title = "") +
    theme_classic() 

print(p)


stat_df <- df %>% 
    group_by(type) %>% 
    summarize(
        mean_value = mean(value),
        sd_value = sd(value)
    )
head(stat_df)

type_order <- c("h1_11mer_intra", "h2_8mer_intra", "h3_7mer_intra", "h1_11mer_vs_h2_8mer", "h1_11mer_vs_h3_7mer", "h2_8mer_vs_h3_7mer")
stat_df$type <- factor(stat_df$type, levels = type_order)


test1 <- wilcox.test(
    value ~ type,
    data = df %>% filter(type %in% c("h1_11mer_intra", "h2_8mer_intra"))
)
cat("P-value for h1_chr13_vs_h3_chr21 vs h1_chr21_vs_h3_chr21:", test1$p.value, "\n")

test2 <- wilcox.test(
    value ~ type,
    data = df %>% filter(type %in% c("h1_11mer_intra", "h2_8mer_intra"))
)
cat("P-value for h2_chr13_vs_h3_chr21 vs h1_chr21_vs_h3_chr21:", test2$p.value, "\n")


p <- ggplot(stat_df, aes(x = type, y = mean_value)) +
    geom_point(shape = 21, size = 4, fill = "white", color = "black") +  # Circle for mean
    geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
                  width = 0.2, color = "black") +  # Error bars for SD
    labs(x = "", y = "Mean ± SD", title = "Mean and SD by Type") +
    theme_classic()
print(p)
ggsave(p, filename="chr8_horstv_sqv.pdf", height = 6, width=6)
