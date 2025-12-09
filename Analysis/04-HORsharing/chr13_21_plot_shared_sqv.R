setwd("/share/home/zhanglab/user/sunyanqing/human/alpha/part2/chr13_21_mnphy/merged_hors")
rm(list = ls())

library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("patchwork")

h1_chr13_vs_h3_chr21 <- read.csv("h1_chr13_vs_h3_chr21_inter.div", sep="\t", header=FALSE)
colnames(h1_chr13_vs_h3_chr21) <- c("type", "seq1", "seq2", "value")
head(h1_chr13_vs_h3_chr21)
print(h1_chr13_vs_h3_chr21$value %>% median())

h2_chr13_vs_h3_chr21 <- read.csv("h2_chr13_vs_h3_chr21_inter.div", sep="\t", header=FALSE)
colnames(h2_chr13_vs_h3_chr21) <- c("type", "seq1", "seq2", "value")
head(h2_chr13_vs_h3_chr21)
print(h2_chr13_vs_h3_chr21$value %>% median())

h1_chr21_vs_chr3_chr21 <- read.csv("h1_chr21_vs_h3_chr21_inter.div", sep="\t", header=FALSE)
colnames(h1_chr21_vs_chr3_chr21) <- c("type", "seq1", "seq2", "value")
head(h1_chr21_vs_chr3_chr21)
print(h1_chr21_vs_chr3_chr21$value %>% median())

df <- rbind(h1_chr13_vs_h3_chr21, h2_chr13_vs_h3_chr21, h1_chr21_vs_chr3_chr21)
print(unique(df$type))
print(df$value %>% median())

stat_df <- df %>% group_by(type) %>% summarize(median_value = median(value), iqr_value = IQR(value))
head(stat_df)

filterdf <- df %>% filter(value <100)


filterdf$type <- factor(filterdf$type, levels = c("h1_chr13_vs_h3_chr21", "h2_chr13_vs_h3_chr21", "h1_chr21_vs_h3_chr21"))

library("ggpubr") 


test1 <- wilcox.test(
    value ~ type,
    data = filterdf %>% filter(type %in% c("h1_chr13_vs_h3_chr21", "h1_chr21_vs_h3_chr21"))
)
cat("P-value for h1_chr13_vs_h3_chr21 vs h1_chr21_vs_h3_chr21:", test1$p.value, "\n")


library(exactRankTests)

# Run the exact Wilcoxon rank-sum test
test_exact <- wilcox.exact(value ~ type, data = filterdf %>%
                             filter(type %in% c("h1_chr13_vs_h3_chr21", "h1_chr21_vs_h3_chr21")))

# Get the p-value
test_exact

test2 <- wilcox.test(
    value ~ type,
    data = filterdf %>% filter(type %in% c("h2_chr13_vs_h3_chr21", "h1_chr21_vs_h3_chr21"))
)
cat("P-value for h2_chr13_vs_h3_chr21 vs h1_chr21_vs_h3_chr21:", test2$p.value, "\n")


p <- ggplot(filterdf, aes(x = type, y = value)) +
    geom_boxplot(fill = "grey", color = "black") +  # Hide outliers
    # geom_violin() +
    labs(x = "", y = "SqVs", title = "") +
    theme_classic() +
    stat_compare_means(comparisons = list(
        c("h1_chr13_vs_h3_chr21", "h1_chr21_vs_h3_chr21"),
        c("h2_chr13_vs_h3_chr21", "h1_chr21_vs_h3_chr21")
    ), method = "wilcox.test", label = "p.format")  + stat_summary(fun = mean, geom = "point", color = "red", size = 2)

print(p)
ggsave(p, filename = "h3_snp_divergence.pdf", width=5, height = 5)


stat_df <- filterdf %>% 
    group_by(type) %>% 
    summarize(
        mean_value = mean(value),
        sd_value = sd(value)
    )
head(stat_df)

# Plot mean with SD as error bars
p <- ggplot(stat_df, aes(x = type, y = mean_value)) +
    geom_point(shape = 21, size = 4, fill = "white", color = "black") +  # Circle for mean
    geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
                  width = 0.2, color = "black") +  # Error bars for SD
    labs(x = "", y = "Mean ± SD", title = "Mean and SD by Type") +
    theme_classic()

print(p)
ggsave(p, filename = "h3_snp_divergence.pdf", width=5, height = 5)


#### 7-mer ####
h1_chr13_vs_h4_chr13 <- read.csv("h1_chr13_vs_h4_chr13_inter.div", sep="\t", header=FALSE)
colnames(h1_chr13_vs_h4_chr13) <- c("type", "seq1", "seq2", "value")
head(h1_chr13_vs_h4_chr13)
print(h1_chr13_vs_h4_chr13$value %>% median())

h2_chr13_vs_h4_chr13 <- read.csv("h2_chr13_vs_h4_chr13_inter.div", sep="\t", header=FALSE)
colnames(h2_chr13_vs_h4_chr13) <- c("type", "seq1", "seq2", "value")
head(h2_chr13_vs_h4_chr13)
print(h2_chr13_vs_h4_chr13$value %>% median())

h1_chr21_vs_h4_chr13 <- read.csv("h1_chr21_vs_h4_chr13_inter.div", sep="\t", header=FALSE)
colnames(h1_chr21_vs_h4_chr13) <- c("type", "seq1", "seq2", "value")
head(h1_chr21_vs_h4_chr13)
print(h1_chr21_vs_h4_chr13$value %>% median())

df <- rbind(h1_chr13_vs_h4_chr13, h2_chr13_vs_h4_chr13, h1_chr21_vs_h4_chr13)
print(unique(df$type))
print(df$value %>% median())

filterdf <- df %>% filter(value <100)
filterdf$type <- factor(filterdf$type, levels = c("h1_chr13_vs_h4_chr13", "h2_chr13_vs_h4_chr13", "h1_chr21_vs_h4_chr13"))

p <- ggplot(filterdf, aes(x = type, y = value)) +
    geom_boxplot(fill = "grey", color = "black") +  # Hide outliers
    # geom_violin() +
    labs(x = "", y = "SqVs", title = "") +
    theme_classic() +
    stat_compare_means(comparisons = list(
        c("h1_chr13_vs_h4_chr13", "h1_chr21_vs_h4_chr13"),
        c("h2_chr13_vs_h4_chr13", "h1_chr21_vs_h4_chr13")
    ), method = "wilcox.test", label = "p.format") 

print(p)

stat_df <- filterdf %>% 
    group_by(type) %>% 
    summarize(
        mean_value = mean(value),
        sd_value = sd(value)
    )
head(stat_df)

# Plot mean with SD as error bars
p <- ggplot(stat_df, aes(x = type, y = mean_value)) +
    geom_point(shape = 21, size = 4, fill = "white", color = "black") +  # Circle for mean
    geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
                  width = 0.2, color = "black") +  # Error bars for SD
    labs(x = "", y = "Mean ± SD", title = "Mean and SD by Type") +
    theme_classic()

print(p)
ggsave(p, filename = "h4_snp_divergence.pdf", width=5, height = 5)

