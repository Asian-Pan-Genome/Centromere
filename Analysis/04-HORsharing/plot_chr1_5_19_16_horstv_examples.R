setwd("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap")
rm(list = ls())

library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("patchwork")

target_df <- read.csv("/share/home/zhanglab/user/sunyanqing/vscode_scripts/plot_scripts/examples_chr22.xls", sep="\t", header=F)
colnames(target_df) <- c("sample_hap_chrom")
target_df <- target_df %>% mutate(y=n()-row_number()+1)
head(target_df)
print(target_df$sample_hap_chrom)


chr<- "chr22"
df_cen <- read.csv("cent_chrom.txt", header=TRUE, sep="\t")
cen_info <- df_cen %>% 
    filter(sample_hap_chrom %in% target_df$sample_hap_chrom) %>% 
    left_join(target_df, by = "sample_hap_chrom") %>%
  select(sample_hap_chrom, sample_hap, start, end, len, project, filterflag, y)
colnames(cen_info) <- c("sample_hap_chrom", "sample_hap", "cen_start", "cen_end", "cen_len", "project", "filterflag", "y")
print(head(cen_info,10))


infile<-paste0(chr, ".HiCAT.horstv.bed")
data <- read.csv(infile, header=FALSE, sep="\t")
colnames(data) <- c("sample_hap_chrom", "horstart", "horend", "hor_mn", "horclass", "horclass_color", "horstv_index", "horstv_color")
print(data[20000, ])


recolor_file <- "all.horstv.hicat.update.color.txt"
recolordf  <- read.csv(recolor_file, header=TRUE, sep="\t")
colnames(recolordf) <- c("horstv_index", "horstv_newcolor")
print(head(recolordf,10))


stv <- data  %>%
  inner_join(cen_info, by = "sample_hap_chrom") %>%
  left_join(recolordf, by = "horstv_index") %>%
  mutate(horstv_newcolor = ifelse(horstv_color == '-', "#525252",
                                  ifelse(!is.na(horstv_newcolor), horstv_newcolor, horstv_color)))  
head(stv, n=50)


stv_min_hor_start <- stv %>%
  group_by(sample_hap_chrom) %>%
  filter(horclass == "94") %>%
  summarise(min_hor_start = min(horstart, na.rm = TRUE),
            max_hor_start = max(horend, na.rm = TRUE)) %>%
  ungroup()
head(stv_min_hor_start,10)

plot_stv_df <- stv %>%
  filter(horstart >= cen_start & horend <= cen_end) %>%
  inner_join(stv_min_hor_start, by = c("sample_hap_chrom")) %>%
  filter(horstart >= min_hor_start & horend <= max_hor_start) %>%
  mutate(
    modstart = if_else(min_hor_start < cen_start, horstart - cen_start + 1, horstart - min_hor_start + 1),
    modend = if_else(min_hor_start < cen_start, horend - cen_start + 1, horend - min_hor_start + 1)
  ) 
head(plot_stv_df,10)
print(unique(plot_stv_df$horstv_newcolor))
print(length(unique(plot_stv_df$sample_hap_chrom)))
print(plot_stv_df %>% filter(horstv_index == "C25H4(8)"))


max_len <- max(plot_stv_df$modend)
y_labels <- plot_stv_df %>%distinct(sample_hap_chrom, y)

p3 <- ggplot() +
  geom_rect(
    data = plot_stv_df,
    aes(xmin = modstart, xmax = modend, ymin = y - 0.3, ymax = y + 0.3, fill = horstv_newcolor),
    color = NA 
  ) +
  scale_fill_identity() + 
  scale_x_continuous(
      breaks = seq(0, max_len, by = 1000000), 
      labels = function(x) x / 1e6 #Mb
    ) +  
  scale_y_continuous(breaks = y_labels$y, labels = y_labels$sample_hap_chrom) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 6)
  )

print(p3)
ggsave("Figures/chr22_examples.cenhap.pdf", plot = p3 , device = "pdf", width = 10, height = 1)

