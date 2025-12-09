setwd("/share/home/zhanglab/user/sunyanqing/human/alpha/part2/chr1_5_16_19_mnphy/merged_hor1-3")
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

type_order <- c("chr1_intra", "chr5_intra", "chr19_intra", "chr19_vs_chr1", "chr19_vs_chr5", "chr5_vs_chr1")
stat_df$type <- factor(stat_df$type, levels = type_order)


test1 <- wilcox.test(
    value ~ type,
    data = df %>% filter(type %in% c("chr1_intra", "chr19_intra"))
)
cat("P-value for chr1_intra vs chr19_intra:", test1$p.value, "\n")

test2 <- wilcox.test(
    value ~ type,
    data = df %>% filter(type %in% c("chr5_intra", "chr19_intra"))
)
cat("P-value for chr5_intra vs chr19_intra:", test2$p.value, "\n")


p <- ggplot(stat_df, aes(x = type, y = mean_value)) +
    geom_point(shape = 21, size = 4, fill = "white", color = "black") +  # Circle for mean
    geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
                  width = 0.2, color = "black") +  # Error bars for SD
    labs(x = "", y = "Mean ± SD", title = "Mean and SD by Type") +
    theme_classic()
print(p)
ggsave(p, filename="hor1-3_horstv_sqv.pdf", height = 6, width=6)
