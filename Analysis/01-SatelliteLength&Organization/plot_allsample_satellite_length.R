setwd("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/satellite_length")
rm(list = ls())

library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("patchwork")
library("pheatmap")
library(ggnewscale)
library(gridExtra)
library(ggthemes)
library(see)
library(scales)

df_cen <- read.csv("../centhap/cent_chrom.txt", header=TRUE, sep="\t")
head(df_cen)

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

complete_df_cen <- df_cen %>%
 filter(filterflag == 0) %>% 
 filter( ! sample_hap_chrom %in% filter_chrs) %>% 
 filter(project != "Ref")
head(complete_df_cen)

chrom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 'chrX', 'chrY')
complete_df_cen$chrom <- factor(complete_df_cen$chrom, levels = chrom_order)


popfile <- "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/populaion.xls"
popdf <- read.csv(popfile, header=TRUE, sep="\t")
popdf <- popdf %>% 
  select(c("sample", "superpopulation")) %>%
  mutate(superpopulation = str_replace(superpopulation, "EAS-APG", "EAS"))
head(popdf,10)


complete_df_cen <- complete_df_cen %>%
  left_join(popdf, by = "sample")
head(complete_df_cen, 10)

complete_df_cen$superpopulation <- factor(complete_df_cen$superpopulation, 
                                          levels = c("AFR", "AMR", "EAS", "EUR", "SAS"))


refdata <- df_cen %>% 
  filter(project == "Ref") %>% 
  mutate(
    color = case_when(
      sample == "CHM13" ~ "#008C8C",#65B1AF#008C8C
      sample == "CN1" ~ "#fdb462",#E85827
      sample == "HG002" ~ "#b3de69",#D8DB15
      sample == "YAO" ~ "#bc80bd"#AD218A
    ),
    shape = case_when(
      sample == "CHM13" ~ 23,  # 菱形
      sample == "CN1" & hap == "Mat" ~ 23,  # 菱形
      sample == "CN1" & hap == "Pat" ~ 22,  # 正方形
      sample == "HG002" & hap == "Mat" ~ 23,  # 菱形
      sample == "HG002" & hap == "Pat" ~ 22,  # 正方形
      sample == "YAO" & hap == "Mat" ~ 23,  # 菱形
      sample == "YAO" & hap == "Pat" ~ 22,  # 正方形
    )
  )
head(refdata, 10)

# Plotting
p <- ggplot(complete_df_cen, aes(chrom, len/1000000)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.3, shape = 21, fill = "white", color = "#525252") +
  geom_point(
    data = refdata,
    aes(chrom, len/1000000, shape = factor(shape), fill = color),
    size = 2.5,
    position = position_nudge(x = 0),
    color = "black"
  ) +
  scale_shape_manual(values = c(22, 23)) +
  scale_fill_identity() + 
  theme_classic()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "plain",color= "black"),
        axis.text.x = element_text(size = 12, face = "plain",color= "black"),
        axis.text.y = element_text(size = 12, face = "plain",color= "black")) + 
  labs(x = "", y = "Centromere Length (Mb)")+
  annotate("point", x = 23.3, y = 40, shape = 22, size = 3, fill = "white", color = "black") + 
  annotate("text", x = 24, y =40 , label = "Pat", size = 3.5, color = "black")+
  annotate("point", x = 23.3, y = 38.5, shape = 23, size = 3, fill = "white", color = "black") + 
  annotate("text", x = 24, y =38.5 , label = "Mat", size = 3.5, color = "black")+
  annotate("point", x = 23.3, y = 37, shape = 21, size = 2.5, fill = "#008C8C", color = "#008C8C") + 
  annotate("text", x = 24, y = 37, label = "CHM13", size = 3.5, color = "black") + 
  annotate("point", x = 23.3, y = 35.5, shape = 21, size = 2.5, fill = "#fdb462", color = "#fdb462") + 
  annotate("text", x = 24, y = 35.5, label = "CN1", size = 3.5, color = "black") +
  annotate("point", x = 23.3, y = 34, shape = 21, size = 2.5, fill = "#b3de69", color = "#b3de69") + 
  annotate("text", x = 24, y = 34, label = "HG002", size = 3.5, color = "black") +
  annotate("point", x = 23.3, y = 32.5, shape = 21, size = 2.5, fill = "#bc80bd", color = "#bc80bd") + 
  annotate("text", x = 24, y = 32.5, label = "YAO", size = 3.5, color = "black")+
  geom_violin(scale = "width",fill = NA, color = "black")

print(p)
ggsave(p, filename = "all_cent_size_violin.pdf", width = 10, height = 4)

superpopulation_colors <- c(
    "AFR" = "#319b62",
    "AMR" = "#939393",
    "EAS" = "#d22e77",
    "EAS-APG" = "#e85827",
    "EUR" = "#0070c0",
    "SAS" = "#893f8b"
)



offsets <- c(
  "AFR" = -0.2,
  "AMR" = -0.1,
  "EAS" = 0,
  "EUR" = 0.1,
  "SAS" = 0.2
)

chrom_levels <- levels(complete_df_cen$chrom)
chrom_map <- setNames(seq_along(chrom_levels), chrom_levels)

complete_df_cen$chrom_numeric <- chrom_map[complete_df_cen$chrom]
refdata$chrom_numeric <- chrom_map[refdata$chrom]
head(complete_df_cen, 10)
head(refdata, 10)

complete_df_cen$chrom_offset <- complete_df_cen$chrom_numeric + offsets[complete_df_cen$superpopulation]
head(complete_df_cen, 10)



p1 <- ggplot() +
  geom_jitter(
    data = complete_df_cen,
    aes(x = chrom_offset, y = len / 1e6, color = superpopulation),
    width = 0,
    height = 0.2,
    alpha = 0.6,
    size = 0.8,
    shape = 19
  ) +
  scale_color_manual(values = superpopulation_colors) +
  scale_x_continuous(
    breaks = seq_along(chrom_levels),
    labels = chrom_levels
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, hjust = 1),
    axis.text.y = element_text(size = 12)
  ) +
  labs(y = "Centromere Length (Mb)")
print(p1)
ggsave(p1, filename = "all_cent_size_superpopulation_violin.pdf", width = 12, height = 4)


library(ggpubr)

p2 <- ggplot(complete_df_cen, aes(x = superpopulation, y = len / 1e6)) +
  geom_boxplot(aes(fill = superpopulation), alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = superpopulation), width = 0.2, alpha = 0.4, size = 0.8) +
  facet_wrap(~ chrom, scales = "free_y") +
  scale_fill_manual(values = superpopulation_colors) +
  scale_color_manual(values = superpopulation_colors) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(y = "Centromere Length (Mb)") +
  stat_compare_means(
    method = "wilcox.test",
    # label = "p.signif",
    label = "p.format",
    comparisons = list(
      c("AFR", "AMR"),
      c("AFR", "EAS"),
      c("AFR", "EUR"),
      c("AFR", "SAS"),
      c("AMR", "EAS"),
      c("AMR", "EUR"),
      c("AMR", "SAS"),
      c("EAS", "EUR"),
      c("EAS", "SAS"),
      c("EUR", "SAS")
    )
  )

print(p2)
ggsave(p2, filename = "boxplot_superpopulation_chrom_significance.pdf", width = 10, height = 16)



###########  check batch effect ###############
head(complete_df_cen, 10)

df_eas_apg <- complete_df_cen %>% 
  filter(project == "APG") %>% 
  mutate(chrom = factor(chrom, levels = chrom_order))
head(df_eas_apg, 10)


df_eas_others <- complete_df_cen %>% 
  filter(project != "APG", superpopulation == "EAS") %>% 
  mutate(chrom = factor(chrom, levels = chrom_order))
head(df_eas_others, 10)

#check complete number of EAS centromere from other project#

df_eas_others_num <- df_eas_others %>%
  group_by(chrom) %>%
  summarise(count = n())
print(df_eas_others_num, n =24)

df_eas_apg <- df_eas_apg %>% 
  mutate(flag = "APG")
df_eas_others <- df_eas_others %>%
 mutate(flag = "HGSVC & HPRC")



df_eas <- bind_rows(df_eas_apg, df_eas_others)
head(df_eas, 10)

library(ggpubr)

flag_color= c("APG" = "#e85827", "HGSVC & HPRC" = "#285943")

p3 <- ggplot(df_eas, aes(x = flag, y = len / 1e6)) +
  geom_boxplot(aes(fill = flag), alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = flag), width = 0.2, alpha = 0.4, size = 0.8) +
  facet_wrap(~ chrom, scales = "free_y") +
  scale_fill_manual(values = flag_color) +
  scale_color_manual(values = flag_color) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(y = "Centromere Length (Mb)") +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    comparisons = list(
      c("APG", "HGSVC & HPRC")
    )
  )+
  geom_text(
    data = df_eas_others_num,  # Use the annotation data frame
    aes(
      x = 2,  # Adjust x position for annotation
      y = 2,  # Place annotation at the top of each facet
      label = paste0("n=", count)
    ),
    inherit.aes = FALSE,
    vjust = 1.5,  # Adjust vertical position
    size = 3
  )

print(p3)
ggsave(p3, filename = "boxplot_apg_vs_eas_batch_effect.pdf", width = 10, height = 16)



##check HGSVC##
hgsvc_df <- read.csv("/share/home/zhanglab/user/sunyanqing/vscode_scripts/satellite_length/hgsvc_chrX_horarray.xls", sep="\t", header=FALSE)
colnames(hgsvc_df) <- c("sample", "tool", "version", "contig", "len", "chrom", "strand")
head(hgsvc_df, 10)

hgsvc_df <- hgsvc_df %>% 
  left_join(popdf, by = "sample")
head(hgsvc_df, 10)

superpopulation_counts <- hgsvc_df %>%
  group_by(superpopulation) %>%
  summarise(count = n())
head(superpopulation_counts, 10)

p <- ggplot(hgsvc_df, aes(x = superpopulation, y = len, fill = superpopulation)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, color = "black") +
  geom_jitter(aes(color = superpopulation), width = 0.2, alpha = 0.6, size = 1) +
  geom_text(
    data = superpopulation_counts,
    aes(x = superpopulation, y = 4.5, label = count),
    size = 3
  ) +
  scale_fill_manual(values = superpopulation_colors) +
  scale_color_manual(values = superpopulation_colors) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(
    x = "Superpopulation",
    y = "Satellite Length"
  )

print(p)
ggsave(p, filename = "HGSVC_selfanno_boxplot_superpopulation_count.pdf", width = 5, height = 3)


##plot chr9##
chr9_complete_cen <- complete_df_cen %>% 
  filter(chrom == "chr9") %>%
  mutate(satellite = "Cent Size") %>%  
  select("sample_hap_chrom", "sample_hap", "sample", "hap", "chrom", "satellite", "len", "project")
head(chr9_complete_cen, 10)
print(length(unique(chr9_complete_cen$sample_hap_chrom)))

chr9_refdata <- refdata %>% 
  filter(chrom == "chr9")
head(chr9_refdata, 10)

sat_df <- read.csv("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/satellite_length/all.sat.length.xls", sep="\t", header=TRUE)
complete_sat_df <- sat_df %>% 
  filter(chrom == "chr9",
         sample_hap_chrom %in% chr9_complete_cen$sample_hap_chrom) 
colnames(complete_sat_df) <- c("sample_hap_chrom", "sample_hap", "sample", "hap", "chrom", "satellite", "len", "project")
head(complete_sat_df)

chr9_combined_df <- bind_rows(chr9_complete_cen, complete_sat_df) %>%
  left_join(popdf, by="sample") %>%
  mutate(satellite = factor(satellite, levels = c("Cent Size", "ASat", "HSat3", "HSat2", "BSat", "GSat")))
head(chr9_combined_df,10)

unique_counts <- chr9_combined_df %>%
  group_by(superpopulation) %>%
  summarise(unique_sample_hap_chrom = n_distinct(sample_hap_chrom))
print(unique_counts)

library(ggpubr)

# Plot
p <- ggplot(chr9_combined_df, aes(x = superpopulation, y = len/1e6, fill = superpopulation)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, color = "black") +
  geom_jitter(aes(color = superpopulation),width = 0.2, alpha = 0.6, size = 1) +
  facet_wrap(~ satellite, scales = "free") +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    comparisons = list(
      c("AFR", "EAS")
    )
  ) +
  scale_fill_manual(values = superpopulation_colors) +
  scale_color_manual(values = superpopulation_colors) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  ) + labs(x ="Superpopulation", y = "Length(Mb)")

print(p)
ggsave(p, filename = "chr9_size_sat_superpopulation_sig.pdf", width = 10, height = 10)




