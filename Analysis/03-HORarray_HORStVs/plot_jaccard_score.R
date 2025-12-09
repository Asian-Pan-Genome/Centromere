setwd("/share/home/zhanglab/user/sunyanqing/human/alpha/vsearch/20241022")
rm(list = ls())

library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")

df <- read.csv("step-4/overlap_jaccard_score_linked.txt", sep="\t", header=FALSE)
colnames(df) <- c("weighted_jaccard_score")


df <- df %>% mutate(identity_threshold = seq(87, 99, by=1))
head(df)


p <- ggplot(df, aes(x = identity_threshold, y = weighted_jaccard_score)) +
    geom_line(color = "black") + 
    geom_point(color = "#b2182b", size = 5) +  # Add dots with larger size
    scale_x_continuous(breaks = seq(87, 99, by = 1)) +  # Set x-axis ticks
    labs(
        x = "Identity Threshold",
        y = "Weighted Jaccard Score"
    ) +
    theme_minimal()  
print(p)
ggsave("step-4/jaccard_score_plot.pdf", plot = p, width = 8, height = 8)
