setwd("/share/home/zhanglab/user/sunyanqing/human/periphy/CHM13/IBS")
rm(list = ls())

library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("patchwork")
library(ggnewscale)

linked_data <- read.csv("chroms_linked_correlation_summary_update.xls", sep="\t", header=TRUE)
head(linked_data)

chrom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
    "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
    "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
    "chr19", "chr20", "chr21", "chr22", "chrX")

linked_data$chrom <- linked_data$chrom <- factor(linked_data$chrom, levels = chrom_order)


###plot left_linked##
left_linked_df <- linked_data %>% filter(type=="left_linked")
print(left_linked_df)

right_linked_df <- linked_data %>% filter(type == "right_linked")
print(right_linked_df)

# Plotting left and right #

p_left_size <- ggplot(left_linked_df, aes(x = chrom, y = length/1000000)) +
    geom_bar(stat = "identity", fill = "#B3D4DE", color="black") +
    scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2)) +
    ylim(0,2) +
    theme_classic() +
    labs(x = "Chromosome", y = "Length") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_left_size)

p_left_r <- ggplot(left_linked_df, aes(x = chrom, y = avg_r)) +
    geom_point(color = "#27418C", size = 2) +  # Plot dots
    geom_line(group = 1, color = "#27418C", size = 1) +  
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    ylim(0,1) +
    theme_classic() +
    labs(x = "Chromosome", y = "Average r") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_left_r)

p_right_size <- ggplot(right_linked_df, aes(x = chrom, y = length/1000000)) +
    geom_bar(stat = "identity", fill = "#85A665", color="black" ) +
    scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2)) +
    ylim(0,2) +
    theme_classic() +
    labs(x = "Chromosome", y = "Length") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_right_size)

p_right_r <- ggplot(right_linked_df, aes(x = chrom, y = avg_r)) +
    geom_point(color = "#254B4F", size = 2) +  # Plot dots
    geom_line(group = 1, color = "#254B4F", size = 1) +  
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    ylim(0,1) +
    theme_classic() +
    labs(x = "Chromosome", y = "Average r") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_right_r)

library(gridExtra)

combined_plot <- grid.arrange(
    p_left_size, 
    p_left_r, 
    p_right_size, 
    p_right_r, 
    nrow = 4, 
    ncol = 1
)
ggsave(combined_plot, file="left_right_linked_pattern_update.pdf", width=8, height=8)


### lr linkage ###
lr_linked_df  <- linked_data %>% filter(type=="lr_left_linked")
head(lr_linked_df)

print(lr_linked_df$avg_r %>% mean())
print(lr_linked_df$avg_r %>% sd())


p_lr_r <- ggplot(lr_linked_df, aes(x = chrom, y = avg_r)) +
    geom_point(color = "#D97E8E", size = 2) +  # Plot dots
    geom_line(group = 1, color = "#D97E8E", size = 1) +  
    geom_hline(yintercept = 0.3, linetype = "dashed", color = "grey") +  # Add dashed grey line
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    ylim(0, 1) +
    theme_classic() +
    labs(x = "Chromosome", y = "Average r") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_lr_r)
ggsave(p_lr_r, file="lr_linked_pattern_update.pdf", width=8, height=2)


###vs cenhap###
left_cent <- linked_data %>% filter(type=="left_cent")
print(left_cent)

left_tmp <- left_cent %>%
    select(chrom, avg_r, type)
head(left_tmp)

right_cent <- linked_data %>% filter(type == "right_cent")
print(right_cent)

right_tmp <- right_cent %>%
    select(chrom, avg_r, type)
head(right_tmp)


combined_cent <- bind_rows(left_tmp, right_tmp)
print(combined_cent)


lr_customized_color <- c("left_cent" = "#60C2C6", 
                        "right_cent" = "#67B82E")
p_lr_cent <- ggplot(combined_cent, aes(x = chrom, y = avg_r, color = type, group = type)) +
    geom_point(size = 2) +  # Plot dots
    geom_line(size = 1) +   # Plot lines 
    geom_hline(yintercept = 0.3, linetype = "dashed", color = "grey") +  
    scale_color_manual(values = lr_customized_color) +
    theme_classic() +
    labs(x = "Chromosome", y = "Average r", color = "Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none")
print(p_lr_cent)
ggsave(p_lr_cent, file="cent_vs_lr_linked_pattern_update.pdf", width=8, height=2)

