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

##load satellite length data and cen info ##
satlen  <- read.csv("all.sat.length.xls", header=TRUE, sep="\t")
df_cen <- read.csv("../centhap/cent_chrom.txt", header=TRUE, sep="\t")
head(df_cen)
df <- satlen %>% left_join(df_cen %>% select(sample_hap_chrom, len, gapless, CL, filterflag), by = "sample_hap_chrom", relationship = "many-to-many")
head(df)

##filtered chromosomes##
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
filtereddf <- df %>% filter(filterflag == 0) %>% filter( ! sample_hap_chrom %in% filter_chrs) %>% filter(satellite != "rDNA")
head(filtereddf)


##chrom order##
chrom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 'chrX', 'chrY')
filtereddf$chrom <- factor(filtereddf$chrom, levels = chrom_order)


##chrom color##
chrom_colors <- c(
  "chr1" = "#9D5427", "chr2" = "#D67C1B", "chr3" = "#C6A15B", "chr4" = "#F2E86D",
  "chr5" = "#FFC000", "chr6" = "#FF6533", "chr7" = "#FF4B3E", "chr8" = "#CB3C23",
  "chr9" = "#972D07", "chr10" = "#a4c558", "chr11" = "#59b669", "chr12" = "#3c8340",
  "chr13" = "#6DAB30", "chr14" = "#285943", "chr15" = "#7ec0b4", "chr16" = "#85A7CC",
  "chr17" = "#3f66a1", "chr18" = "#3897c5", "chr19" = "#afa7d8", "chr20" = "#a874b5",
  "chr21" = "#893f8b", "chr22" = "#b13f73", "chrX" = "#574D68", "chrY" = "#C52B94"
)

##satellite color##
satellite_colors <- c(
  "ASat" = "#891640", "HSat1A" = "#18988B", "HSat1B" = "#18988B", "HSat1" = "#18988B", "HSat2" = "#323366",
  "HSat3" = "#77A0BA", "BSat" = "#DCAED0", "GSat" = "#758660"
)

#############################################################################################################################
##plot boxplot of percentage of satellite in each haplotype##
apg_sat_percent <- read.csv("apg_satellite_precentage.xls", header=TRUE, sep="\t")
ref_sat_percent <- read.csv("ref_satellite_precentage.xls", header=TRUE, sep="\t")
head(apg_sat_percent)
head(ref_sat_percent)

filter_sample_hap <- c(  'C020-CHA-E20_Mat',
  'C045-CHA-N05_Mat',
  'C019-CHA-E19_Pat',
  'C076-CHA-NE16_Pat')
plot_apg_sat_percent <- apg_sat_percent %>% 
  select(sample_hap, ASat, BSat, HSat1, HSat2, HSat3, GSat) %>%
  filter(! sample_hap %in% filter_sample_hap)
plot_ref_sat_precent <- ref_sat_percent %>% 
  select(sample_hap, ASat, BSat, HSat1, HSat2, HSat3, GSat) 

apg_sat_long <- gather(plot_apg_sat_percent, key="satellite", value="percentage", -1)
ref_sat_long <- gather(plot_ref_sat_precent, key="satellite", value="percentage", -1)
head(apg_sat_long)
head(ref_sat_long)
satellite_order <- c("ASat", "HSat1", "HSat2", "HSat3", "BSat", "GSat")
apg_sat_long$satellite <- factor(apg_sat_long$satellite, levels = satellite_order)
ref_sat_long$satellite <- factor(ref_sat_long$satellite, levels = satellite_order)

ref_sat_long <- ref_sat_long %>%
  mutate(
    sample_hap = factor(sample_hap, levels = c("YAO_Mat", "YAO_Pat", "HG002_Mat", "HG002_Pat", "CN1_Mat", "CN1_Pat", "CHM13_-")),
    color = case_when(
      sample_hap == "CHM13_-" ~ "#008C8C",
      sample_hap == "CN1_Mat" ~ "#fdb462",
      sample_hap == "CN1_Pat" ~ "#fdb462",
      sample_hap == "HG002_Mat" ~ "#b3de69",
      sample_hap == "HG002_Pat" ~ "#b3de69",
      sample_hap == "YAO_Mat" ~ "#bc80bd",
      sample_hap == "YAO_Pat" ~ "#bc80bd"
    ),
    shape = case_when(
      sample_hap == "CHM13_-" ~ 23,  # 菱形
      sample_hap == "CN1_Mat" ~ 23,  # 菱形
      sample_hap == "CN1_Pat" ~ 22,  # 正方形
      sample_hap == "HG002_Mat" ~ 23,  # 菱形
      sample_hap == "HG002_Pat" ~ 22,  # 正方形
      sample_hap == "YAO_Mat" ~ 23,  # 菱形
      sample_hap == "YAO_Pat" ~ 22,  # 正方形
    )
  )%>% arrange(sample_hap)
head(ref_sat_long)

##half violin plot of satellite percentage##
percent_sat_p <- ggplot(apg_sat_long, aes(satellite, percentage, fill = satellite)) +
  geom_violinhalf(scale = "width", position = position_nudge(x = 0.2, y = 0)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.3, shape = 21, fill = "white", color = "#525252") +
  theme_modern() +
  scale_fill_manual(values = satellite_colors) +
  new_scale("fill") + 
  geom_point(
    data = ref_sat_long,
    aes(satellite, percentage, shape = factor(shape), fill = color, color = color),
    size = 2.5,
    position = position_nudge(x = 0),
    color = "black",
    stroke = 0.1
  ) +
  scale_shape_manual(values = c(22, 23)) +
  scale_fill_identity() + 
  scale_color_identity() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, face = "plain", color = "black"),
    axis.text.x = element_text(size = 8, face = "plain", color = "black"),
    axis.text.y = element_text(size = 8, face = "plain", color = "black")
  ) + labs(x = "", y = "% of haploid assembly")+
  annotate("point", x = 5.5, y = 3.6, shape = 22, size = 3, fill = "white", color = "black") + 
  annotate("text", x = 6, y =3.6 , label = "Pat", size = 3, color = "black")+
  annotate("point", x = 5.5, y = 3.45, shape = 23, size = 3, fill = "white", color = "black") + 
  annotate("text", x = 6, y =3.45 , label = "Mat", size = 3, color = "black")+
  annotate("point", x = 5.5, y = 3.3, shape = 21, size = 3, fill = "#008C8C", color = "#008C8C") + 
  annotate("text", x = 6, y = 3.3, label = "CHM13", size = 3, color = "black") + 
  annotate("point", x = 5.5, y = 3.15, shape = 21, size = 3, fill = "#fdb462", color = "#fdb462") + 
  annotate("text", x = 6, y = 3.15, label = "CN1", size = 3, color = "black") +
  annotate("point", x = 5.5, y = 3.0, shape = 21, size = 3, fill = "#b3de69", color = "#b3de69") + 
  annotate("text", x = 6, y = 3.0, label = "HG002", size = 3, color = "black") +
  annotate("point", x = 5.5, y = 2.85, shape = 21, size = 3, fill = "#bc80bd", color = "#bc80bd") + 
  annotate("text", x = 6, y = 2.85, label = "YAO", size = 3, color = "black")

print(percent_sat_p)
outpdf <- paste0("apg_satellite_percentage_halfviolin_marker_nolegend.pdf")
ggsave(outpdf, plot = percent_sat_p, device = "pdf", width = 4.5, height = 2)

#############################################################################################################################
##each satellite length for chromsomes##
head(filtereddf)
print(unique(filtereddf$satellite))
print(filtereddf %>% filter(is.na(satellite)) )

##population info##
popdf <- read.csv("../centhap/populaion.xls", sep = "\t", header = TRUE)
head(popdf)

##merge population info##
filtereddf <- filtereddf %>% merge(popdf %>% select(sample, population, superpopulation), by="sample")
head(filtereddf)
print(unique(filtereddf$satellite))
print(filtereddf %>% filter(is.na(satellite)))

##add HSat1 length##
hsat1_df <- filtereddf %>%
  filter(satellite %in% c("HSat1A", "HSat1B")) %>%
  group_by(sample, sample_hap_chrom, sample_hap, hap, chrom, project, gapless, CL, filterflag, population, superpopulation) %>%
  summarize(length = sum(length)) %>%
  mutate(satellite = "HSat1")
head(hsat1_df)

##merge HSat1 length and delete HSat1A and HSat1B##
filtereddf <- bind_rows(filtereddf, hsat1_df) %>% filter(satellite != "HSat1A" & satellite != "HSat1B")
apg_filtereddf <- filtereddf %>% filter(project == "APG")
head(filtereddf)
print(unique(filtereddf$satellite))
print(filtereddf %>% filter(is.na(satellite)))
print(apg_filtereddf %>% filter(is.na(satellite)))

##satellite order##
satellite_order <- c("ASat", "HSat1", "HSat2", "HSat3", "BSat", "GSat")
filtereddf$satellite <- factor(filtereddf$satellite, levels = satellite_order)
apg_filtereddf$satellite <- factor(apg_filtereddf$satellite, levels = satellite_order)

##superpopulation_color and project_shape##
superpopulation_color <- c("AFR" = "#319b62", "AMR" = "#939393", "EAS" = "#D22E77", "EAS-APG" = "#E85827", "EUR" = "#0070c0", "SAS" = "#893f8b")
project_shape <- c("HPRC" = 24, "HGSVC" = 22, "Ref" = 23, "APG" = 21)

ref_sat <- filtereddf %>% filter(project == "Ref") %>%
mutate(sample_hap = factor(sample_hap, levels = c("YAO_Mat", "YAO_Pat", "HG002_Mat", "HG002_Pat", "CN1_Mat", "CN1_Pat", "CHM13_-")),
color = case_when(
      sample_hap == "CHM13_-" ~ "#008C8C",
      sample_hap == "CN1_Mat" ~ "#fdb462",
      sample_hap == "CN1_Pat" ~ "#fdb462",
      sample_hap == "HG002_Mat" ~ "#b3de69",
      sample_hap == "HG002_Pat" ~ "#b3de69",
      sample_hap == "YAO_Mat" ~ "#bc80bd",
      sample_hap == "YAO_Pat" ~ "#bc80bd"
    ),
    shape = case_when(
      sample_hap == "CHM13_-" ~ 23, 
      sample_hap == "CN1_Mat" ~ 23,  
      sample_hap == "CN1_Pat" ~ 22,  
      sample_hap == "HG002_Mat" ~ 23, 
      sample_hap == "HG002_Pat" ~ 22, 
      sample_hap == "YAO_Mat" ~ 23, 
      sample_hap == "YAO_Pat" ~ 22,  
    )) %>% arrange(sample_hap)
head(ref_sat)


all_avg_asat <- filtereddf %>% 
  filter(satellite == "ASat") %>% 
  .$length %>% 
  mean()
head(all_avg_asat)

apg_avg_asat <- apg_filtereddf %>% 
  filter(satellite == "ASat") %>% 
  .$length %>% 
  mean()
head(apg_avg_asat)

head(apg_filtereddf)
print(unique(apg_filtereddf$satellite))

apg_stat <- apg_filtereddf %>% 
  group_by(chrom, satellite) %>%
  summarise(mean_length = mean(length),
            sd_length = sd(length),
            min_len = min(length),
            max_length = max(length))
print(apg_stat)
write.table(apg_stat, file = "apg_satellite_length_stat.xls", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


violin_plot <- ggplot(apg_filtereddf, aes(x = chrom, y = length / 1000000)) +
  geom_violin(aes(color = chrom), scale = "width", alpha = 0.8, draw_quantiles = c(0.25, 0.5, 0.75)) +
  facet_wrap(~ satellite, ncol = 1, scales = "free_y") +
  scale_color_manual(values = chrom_colors) +  
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "plain", color = "black"),
    axis.text.x = element_text(size = 12, face = "plain", color = "black"),
    axis.text.y = element_text(size = 12, face = "plain", color = "black"),
    strip.text.x = element_text(colour = "black", face = "bold")
  ) + 
  labs(x = "", y = "Length (Mb)") +
  geom_hline(data = filtereddf %>% filter(satellite == "ASat"), 
             aes(yintercept = apg_avg_asat / 1000000), 
             color = "grey", linetype = "dashed") +
  new_scale_color() + 
  geom_point(
    data = ref_sat,
    aes(x = chrom, y = length / 1000000, shape = factor(shape), fill = color, color = color),
    size = 1.5,
    position = position_nudge(x = 0.4)
  ) +
  scale_shape_manual(values = c(22, 23)) +
  scale_fill_identity() +
  scale_color_identity()

print(violin_plot)
outpdf <- paste0("apg_satellite_length_variation.pdf")
ggsave(outpdf, plot = violin_plot, device = "pdf", width = 16, height = 16)


############################################################################################################################
### superpopulation ##
head(filtereddf)

pop_metadf <- filtereddf %>%
  mutate(superpopulation = ifelse(superpopulation == "EAS-APG", "EAS", superpopulation))
head(pop_metadf)


head(pop_metadf)


library(ggpubr)
asat_metadf <- pop_metadf %>% filter(satellite == "ASat", !chrom %in% c("chrY"))
asat_metadf <- pop_metadf %>% filter(satellite == "ASat")
superpopulation_colors <- c(
  "AFR" = "#319b62",
  "AMR" = "#939393",
  "EAS" = "#d22e77",
  "EUR" = "#0070c0",
  "SAS" = "#893f8b"
)

stat_sample_num <- asat_metadf %>% 
  group_by(chrom, superpopulation) %>%
  summarise(count = n_distinct(sample_hap))
head(stat_sample_num)


boxplot_asat <- ggplot(asat_metadf, aes(x = superpopulation, y = length / 1000000, fill = superpopulation)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(aes(color = superpopulation), width = 0.2, size = 1, alpha = 0.3) +  # Add jitter for individual points
  facet_wrap(~chrom, nrow = 2) +  # Facet by chromosome
  scale_color_manual(values = superpopulation_colors) + 
  scale_fill_manual(values = superpopulation_colors) +  # Apply custom colors
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "plain", color = "black"),
    axis.text.x = element_text(size = 12, face = "plain", color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, face = "plain", color = "black"),
    strip.text.x = element_text(colour = "black", face = "bold")
  ) +
  labs(x = "Superpopulation", y = "Length")+
  stat_compare_means(
    comparisons = list(
      c("AFR", "AMR"),
      c("AFR", "EAS"),
      c("AFR", "EUR"),
      c("AFR", "SAS")
    ),
    method = "wilcox.test",  # Wilcoxon rank-sum test
    label = "p.format",  # Display p-values
    # label.y = c(9,10,11,12)  # Adjust y-position for labels if needed
    label.y = c(24, 26, 28, 30)
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  geom_text(
    data = stat_sample_num,  # Use the count data
    aes(x = superpopulation, y = 24, label = count),  # Adjust y position
    inherit.aes = FALSE,
    size = 4,
    color = "black")+
  coord_cartesian(ylim = c(0,32))

# Print the plot
print(boxplot_asat)
ggsave(boxplot_asat, filename="Hsat3_size_superpopulation_sig.pdf", width=12, height = 6)


##batch effect##
head(pop_metadf)
print(unique(pop_metadf$project))

asat_metadf_batch <- pop_metadf %>% 
  filter(project == "APG" | superpopulation == "EAS") %>%
  filter(satellite == "HSat3") %>%
  mutate(batchflag = ifelse(project == "APG", "APG", "non-APG"))
head(asat_metadf_batch)

stat_asat_batch <- asat_metadf_batch %>%
  group_by(chrom, batchflag) %>%
  summarise(count = n_distinct(sample_hap))

flag_color= c("APG" = "#e85827", "non-APG" = "#285943")

boxplot_batch <- ggplot(asat_metadf_batch, aes(x = batchflag, y = length / 1000000, fill = batchflag)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(aes(color = batchflag), width = 0.2, size = 1, alpha = 0.3) +  # Add jitter for individual points
  facet_wrap(~chrom, nrow = 2) +  # Facet by chromosome
  scale_color_manual(values = flag_color) + 
  scale_fill_manual(values = flag_color) +  # Apply custom colors
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "plain", color = "black"),
    axis.text.x = element_text(size = 12, face = "plain", color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, face = "plain", color = "black"),
    strip.text.x = element_text(colour = "black", face = "bold")
  ) +
  labs(x = "Superpopulation", y = "Length")+
  stat_compare_means(
    comparisons = list(
      c("APG", "non-APG")
    ),
    method = "wilcox.test",  # Wilcoxon rank-sum test
    label = "p.format",  # Display p-values
    label.y = c(13)  # Adjust y-position for labels if needed
    # label.y = c(25, 26, 27, 28)
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  geom_text(
    data = stat_asat_batch,  # Use the count data
    aes(x = batchflag, y = 12, label = count),  # Adjust y position
    inherit.aes = FALSE,
    size = 4,
    color = "black")+
  coord_cartesian(ylim = c(0, 15))
print(boxplot_batch)
ggsave(boxplot_batch, filename="asat_size_superpopulation_sig_batch.pdf", width=10, height = 5)

#############################################################################################################################
##mat_vs_pat##
diploiddf <- read.csv("mat_vs_pat_constrain/datasets/apg_diploid_mat_pat_length.xls",header = TRUE, sep = "\t")
# diploiddf <- read.csv("mat_vs_pat_constrain/datasets/apg_diploid_mat_pat_relative_length.xls",header = TRUE, sep = "\t")
# diploiddf <- read.csv("mat_vs_pat_constrain/datasets/apg_diploid_mat_pat_ratio.xls",header = TRUE, sep = "\t")
diploiddf <- diploiddf %>% filter(filterflag_Mat == 0 & filterflag_Pat == 0)
diploiddf$chrom <- factor(diploiddf$chrom, levels = chrom_order)
satellite_order <- c("ASat", "HSat1", "HSat2", "HSat3", "BSat", "GSat")
diploiddf$satellite <- factor(diploiddf$satellite, levels = satellite_order)


head(diploiddf)
print(unique(diploiddf$satellite))



#####################################only plot ASat####################################
asat_diploiddf <- diploiddf %>% filter(satellite == "ASat") 
asat_diploiddf$chrom <- factor(asat_diploiddf$chrom, levels = rev(chrom_order))
asat_diploiddf$chrom <- factor(asat_diploiddf$chrom, levels = chrom_order)
head(asat_diploiddf)

###stat asat difference between mat_vs_pat###
mean_v <- mean(abs(asat_diploiddf$value))
sd_v <- mean(abs(asat_diploiddf$value))
lower_bound <- mean_v - 1.96*sd_v
upper_bound <- mean_v + 1.96*sd_v
cat(lower_bound, upper_bound)



check_tmp <- asat_diploiddf %>% filter(value >-2507229 & value < 2507229 )
dim(check_tmp)
dim(cen_diploiddf)

density_plot <- ggplot(asat_diploiddf, aes(x = abs(value)/1000000)) +
  geom_density(fill = "grey", alpha = 0.5) +  
  geom_vline(xintercept = 2.51, linetype = "dashed", color = "black", size = 1) +
  theme_classic() +  # Use a clean theme
  labs(
    title = "",
    x = "aSat array difference (Mb)",
    y = "Density"
  ) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
print(density_plot)
ggsave(density_plot, filename = "mat_vs_pat_constrain/asat_difference_density.pdf", device = "pdf", width = 5, height = 5)

p2 <- ggplot(asat_diploiddf, aes(x = chrom, y = value/1000000, color = chrom)) +
  # stat_summary(fun = median, geom = "point", shape = 21, size =4, aes(color = chrom), stroke = 1.5) +
  # stat_summary(fun = median, geom = "point", shape = 21, size = 1, aes(fill = chrom), stroke = 1.5) +
  stat_summary(fun.data = function(y) {
    mean_y <- mean(y)
    sd_y <- sd(y)
    ymin <- mean_y - sd_y
    ymax <- mean_y + sd_y
    data.frame(y = mean_y, ymin = ymin, ymax = ymax)
  }, geom = "errorbar", width = 0.2, size = 1) +
  # geom_boxplot(outlier.shape = NA) +  # Remove default outliers
  stat_summary(fun.data = function(y) {
    outliers <- boxplot.stats(y)$out
    data.frame(y = outliers)
  }, geom = "point", shape = 4, size = 2, aes(fill = chrom)) +  # Shape 4 corresponds to "X"
  stat_summary(fun = median, geom = "point", shape = 21, size = 4, aes(fill = I("white")), stroke = 1.5) +
  scale_color_manual(values = chrom_colors) +
  scale_fill_manual(values = chrom_colors) +
  # geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 2.51, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -2.51, linetype = "dashed", color = "grey") +
  labs(x = "",
       y = "ASat length (Mat - Pat) (Mb)") +
  scale_y_continuous(breaks = seq(-14, 14, by = 2)) +
  theme_classic() +
  theme(legend.position = "none",
    # axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "plain", color = "black"),
    axis.text.x = element_text(size = 10, face = "plain", color = "black"),
    axis.text.y = element_text(size = 10, face = "plain", color = "black"),
    strip.text.x = element_text(colour = "black", face = "bold"))
    # coord_flip()
print(p2)
outpng <- paste0("mat_vs_pat_constrain/asat/", "apg_mat_vs_pat_diff_len.pdf" )
ggsave(outpng, plot = p2, device = "pdf", width = 5, height = 2)

########################whole centromere ##################################################
library(ggbreak)


cen_diploiddf <- df_cen %>% filter(filterflag ==0, project=="APG")%>%
  select(sample, hap, chrom, len, filterflag) %>%
  group_by(sample, hap, chrom) %>%
  pivot_wider(names_from = hap, values_from = c("len", "filterflag")) %>%
  filter(filterflag_Mat == 0 & filterflag_Pat == 0) %>%
  mutate(value= len_Mat - len_Pat)
head(cen_diploiddf)

cen_diploiddf$chrom <- factor(cen_diploiddf$chrom, levels = chrom_order)
print(max(cen_diploiddf$value))

#######abs value####
mean_v <- mean(abs(cen_diploiddf$value))
sd_v <- mean(abs(cen_diploiddf$value))
lower_bound <- mean_v - 1.96*sd_v
upper_bound <- mean_v + 1.96*sd_v
cat(lower_bound, upper_bound)

check_tmp <- cen_diploiddf %>% filter(value >-4330018 & value < 4330018 )
dim(check_tmp)
dim(cen_diploiddf)

stat_tmp <- cen_diploiddf %>%
  group_by(chrom) %>%
  summarise(mean_value = mean(value), sd_value = sd(value))
print(stat_tmp, n=24)


library(ggplot2)

# Density plot for the `value` column
density_plot <- ggplot(cen_diploiddf, aes(x = abs(value)/1000000)) +
  geom_density(fill = "grey", alpha = 0.5) +  
  geom_vline(xintercept = 4.33, linetype = "dashed", color = "black", size = 1) +
  theme_classic() +  # Use a clean theme
  labs(
    title = "",
    x = "centromere size difference (Mb)",
    y = "Density"
  ) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
print(density_plot)
ggsave(density_plot, filename = "mat_vs_pat_constrain/cent_size_difference_density.pdf", device = "pdf", width = 5, height = 5)

write.csv(cen_diploiddf, "mat_vs_pat_constrain/apg_mat_vs_pat_diff_centromere.xls", row.names = FALSE)

p3 <- ggplot(cen_diploiddf, aes(x = chrom, y = value/1000000, color = chrom)) +
  # stat_summary(fun = median, geom = "point", shape = 21, size =4, aes(color = chrom), stroke = 1.5) +
  # stat_summary(fun = median, geom = "point", shape = 21, size = 1, aes(fill = chrom), stroke = 1.5) +
  stat_summary(fun.data = function(y) {
    mean_y <- mean(y)
    sd_y <- sd(y)
    ymin <- mean_y - sd_y
    ymax <- mean_y + sd_y
    q1 <- quantile(y, 0.25)
    q3 <- quantile(y, 0.75)
    data.frame(y = mean_y, ymin = ymin, ymax = ymax)
  }, geom = "errorbar", width = 0.1, size = 0.3) +
  # geom_boxplot(outlier.shape = NA) +  # Remove default outliers
  stat_summary(fun.data = function(y) {
    outliers <- boxplot.stats(y)$out
    data.frame(y = outliers)
  }, geom = "point", shape = 4, size = 2, aes(fill = chrom)) +  # Shape 4 corresponds to "X"
  stat_summary(fun = median, geom = "point", shape = 21, size = 4, aes(fill = I("white")), stroke = 0.3) +
  scale_color_manual(values = chrom_colors) +
  scale_fill_manual(values = chrom_colors) +
  # geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  # geom_hline(yintercept = 2, linetype = "dashed", color = "grey") +
  # geom_hline(yintercept = -2, linetype = "dashed", color = "grey") +
  labs(x = "",
       y = "centromere length (Mat - Pat) (Mb)") +
  scale_y_continuous(breaks = seq(-15, 15, by = 5)) +
  coord_cartesian(ylim = c(-16, 16)) +
  theme_classic() +
  theme(legend.position = "none",
    # axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, face = "plain", color = "black"),
    axis.text.x = element_text(size = 8, face = "plain", color = "black"),
    axis.text.y = element_text(size = 8, face = "plain", color = "black"),
    strip.text.x = element_text(colour = "black", size = 10, face = "plain")) +
  scale_y_break(c(6, 10), space = 0.2, scales = 1) +
  scale_y_break(c(-10, -6), space = 0.2, scales = 1.5)
print(p3)
outpng <- paste0("mat_vs_pat_constrain/", "apg_mat_vs_pat_diff_centromere_break.pdf" )
ggsave(outpng, plot = p3, device = "pdf", width = 8, height = 4)

########################all satellites#####################################################
p2 <- ggplot(diploiddf, aes(x = chrom, y = value, color = chrom)) +
  # stat_summary(fun = median, geom = "point", shape = 21, size =4, aes(color = chrom), stroke = 1.5) +
  # stat_summary(fun = median, geom = "point", shape = 21, size = 1, aes(fill = chrom), stroke = 1.5) +
  stat_summary(fun.data = function(y) {
    mean_y <- mean(y)
    sd_y <- sd(y)
    ymin <- mean_y - sd_y
    ymax <- mean_y + sd_y
    data.frame(y = mean_y, ymin = ymin, ymax = ymax)
  }, geom = "errorbar", width = 0.2, size = 1) +
  # geom_boxplot(outlier.shape = NA) +  # Remove default outliers
  stat_summary(fun.data = function(y) {
    outliers <- boxplot.stats(y)$out
    data.frame(y = outliers)
  }, geom = "point", shape = 4, size = 2, aes(fill = chrom)) +  # Shape 4 corresponds to "X"
  stat_summary(fun = median, geom = "point", shape = 21, size = 4, aes(fill = I("white")), stroke = 1.5) +
  scale_color_manual(values = chrom_colors) +
  scale_fill_manual(values = chrom_colors) +
  labs(x = "",
       y = "Length difference (Mat - Pat) (Mb)") +
  theme_bw() +
  theme(legend.position = "none",
    # axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "plain", color = "black"),
    axis.text.x = element_text(size = 10, face = "plain", color = "black"),
    axis.text.y = element_text(size = 10, face = "plain", color = "black"),
    strip.text.x = element_text(colour = "black", face = "bold")) +
  facet_wrap(~satellite, scale="free_y", ncol = 1)
  # geom_hline(data = diploiddf %>% filter(satellite %in% satellite_order[1:4]), 
  #            aes(yintercept = 2), linetype = "dashed", color = "grey") +
  # geom_hline(data = diploiddf %>% filter(satellite %in% satellite_order[1:4]), 
  #            aes(yintercept = -2), linetype = "dashed", color = "grey")
    # coord_flip()

print(p2)
outpng <- paste0("mat_vs_pat_constrain/", "apg_satellite_mat_vs_pat_diff_len.pdf" )
ggsave(outpng, plot = p2, device = "pdf", width = 16, height = 16)

############compare extremum between population and mat_vs_pat #######################################################
head(filtereddf)
head(diploiddf)
population_extremum <- filtereddf %>%
  group_by(chrom, satellite) %>%
  pivot_wider(names_from = satellite, values_from = length) %>%
  mutate(ASat = replace_na(ASat, 0),
         HSat1A = replace_na(HSat1A, 0),
         HSat1B = replace_na(HSat1B, 0),
         HSat2 = replace_na(HSat2, 0),
         HSat3 = replace_na(HSat3, 0),
         BSat = replace_na(BSat, 0),
         GSat = replace_na(GSat, 0)) %>%
  mutate(HSat1 = HSat1A + HSat1B) %>%
  group_by(chrom) %>%
  summarise(ASat = max(ASat) - min(ASat),
            HSat1 = max(HSat1) - min(HSat1),
            HSat2 = max(HSat2) - min(HSat2),
            HSat3 = max(HSat3) - min(HSat3),
            BSat = max(BSat) - min(BSat),
            GSat = max(GSat) - min(GSat))%>%
  pivot_longer(cols = c(ASat, HSat1, HSat2, HSat3, BSat, GSat),
               names_to = "satellite",
               values_to = "pop_extremum")
head(population_extremum)


diploid_extremun <- diploiddf %>%
  mutate(length_Mat = replace_na(length_Mat, 0),
         length_Pat = replace_na(length_Pat, 0),
         value = length_Mat - length_Pat) %>%
  filter(value !=0) %>%
  group_by(chrom, satellite) %>%
  summarise(diploid_extremun = max(abs(value)))
head(diploid_extremun)

vs_extremum <- merge(population_extremum, diploid_extremun, by = c("chrom", "satellite"), all=FALSE)
vs_extremum <- vs_extremum %>% mutate(diff = diploid_extremun - pop_extremum)
head(vs_extremum)
write.csv(vs_extremum, "mat_vs_pat_constrain/apg_pop_vs_diploid_extremum.csv", row.names = FALSE)

filtered_vs_extremum <- vs_extremum %>%
  filter(diff > 0)
head(filtered_vs_extremum)

#########distribution of satellite length difference between Mat and Pat##############
simulation_df <- read.csv("mat_vs_pat_constrain/datasets/apg_random_vs_real_length_diff_iter5k.xls",header = TRUE, sep = "\t")
head(simulation_df)

observed_df <- simulation_df %>% filter(type=="observed") %>% filter(! is.na(value))
head(observed_df)
random_df <- simulation_df %>% filter(type=="random") %>% mutate(value = replace_na(value, 0))
head(random_df)

observed_df_count <- observed_df %>% group_by(satellite, chrom) %>% summarise(count = n()) 
head(observed_df_count)

sampled_dfs <- lapply(1:nrow(observed_df_count), function(i) {
  row <- observed_df_count[i, ] 
  subdf <- random_df %>% filter(chrom == row$chrom, satellite == row$satellite) 
  sampled_subdf <- subdf %>% sample_n(200, replace = TRUE) 
  sampled_subdf <- subdf %>% sample_n(1000, replace = TRUE) 
  return(sampled_subdf)
})

random_sampled_df <- do.call(rbind, sampled_dfs)
head(random_sampled_df)

combined_df <- bind_rows(random_sampled_df, observed_df)
combined_df$chrom <- factor(combined_df$chrom, levels = chrom_order)
satellite_order <- c("ASat", "HSat1", "HSat2", "HSat3", "BSat", "GSat")
combined_df$satellite <- factor(combined_df$satellite, levels = satellite_order)

ggplot(combined_df, aes(x = chrom, y = value, fill = type)) +
  geom_violinhalf(scale="width") +
  facet_wrap(~satellite, ncol = 1, scale="free_y") +
  theme_minimal() +
  labs(title = "Half Violin Plot of Satellite Lengths",
       x = "Chromosome",
       y = "Value",
       fill = "Type")

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1, "group"]
  newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
      1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


custom_colors <- c("random" = "#FCBC29", "observed" = "#36907E")

p4 <- ggplot(combined_df, aes(x = chrom, y = value / 1000000, fill = type)) +
  geom_split_violin(scale = "width", alpha = 0.8) +
  facet_wrap(~satellite, ncol = 1, scale = "free_y") +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  labs(x = "",
       y = "Length difference (Mat - Pat) (Mb)")+
  theme(
    axis.title.y = element_text(size = 10, face = "plain", color = "black"),
    axis.text.x = element_text(size = 10, face = "plain", color = "black"),
    axis.text.y = element_text(size = 10, face = "plain", color = "black"),
    strip.text.x = element_text(colour = "black", face = "bold")) 
print(p4)
outpng <- paste0("mat_vs_pat_constrain/", "apg_satellite_observed_vs_random_iter02.pdf" )
ggsave(outpng, plot = p4, device = "pdf", width = 16, height = 16)
