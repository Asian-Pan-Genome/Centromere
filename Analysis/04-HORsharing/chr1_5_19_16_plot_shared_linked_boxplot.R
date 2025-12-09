setwd("/share/home/zhanglab/user/sunyanqing/human/alpha/part2/chr1_5_16_19_mnphy/")
rm(list = ls())
library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("patchwork")

df <- read.csv("HiCAT_target_shared_stv_hicat_layered_count.xls", header=TRUE, sep="\t")
head(df)

target_stvs <- unique(df$reorder_hor)
print(target_stvs)

df <- df %>% 
    filter(chromosome %in% c("chr1", "chr5", "chr19", "chr16")) %>%
    mutate(chromosome = factor(chromosome, levels = c(
        "chr1", "chr5", "chr19", "chr16"))) %>%
    mutate(sample_hap = paste0(sample, "_", hap))


chromosomes <- unique(df$chromosome)
sample_haps <- unique(df$sample_hap)

df <- df %>% 
    complete(sample_hap = sample_haps, chromosome = chromosomes, fill = list(nrepeat = 0)) 


tmp_stv <-"21663_29601"

for (tmp_stv in target_stvs) {
    subdf <- df %>%
        filter(reorder_hor == tmp_stv) 
    head(subdf)

    chrom_color <- c(
        "chr1" = "#7fc97f", 
        "chr5" = "#fdc086",
        "chr19" = "#386cb0", 
        "chr16" = "#beaed4"
    )

    p <- ggplot(subdf, aes(x = chromosome, y = nrepeat)) +
        geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +  # Boxplot
        geom_point(aes(color = chromosome), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.4) +  # Jittered points
        geom_line(aes(group = sample_hap), color = "gray", alpha = 0.2) + 
        scale_color_manual(values =chrom_color )+
        labs(x = "", y = "HOR StVs repeat count", title = tmp_stv) +
        theme_classic() +
        theme(legend.position = "none",
        plot.title = element_text(size = 6)) 
    print(p)

    outpdf <- paste0("HiCAT_shared_stv_figures/", tmp_stv, ".shared_boxplot.pdf")
    print(outpdf)
    ggsave(p, filename=outpdf, width=5, height=5)
}
