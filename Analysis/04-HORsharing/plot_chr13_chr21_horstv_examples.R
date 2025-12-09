setwd("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap")
rm(list = ls())

library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("patchwork")

##load examples file
target_df <- read.csv("examples_chr13_chr21_id.xls", sep="\t", header=F)
colnames(target_df) <- c("sample_hap_chrom")

target_df <- target_df %>% filter(grepl("chr21",sample_hap_chrom))
head(target_df)

target_df <- target_df %>% mutate(y=n()-row_number()+1)
head(target_df)
print(target_df$sample_hap_chrom)

##load cent info
target_chrs <- c("chr13", "chr21")

df_cen <- read.csv("cent_chrom.txt", header=TRUE, sep="\t")
cen_info <- df_cen %>% 
    filter(chrom %in% target_chrs, sample_hap_chrom %in% target_df$sample_hap_chrom) %>% 
    left_join(target_df, by = "sample_hap_chrom") %>%
  select(sample_hap_chrom, sample_hap, start, end, len, project, filterflag, y)
colnames(cen_info) <- c("chrom", "sample_hap", "cen_start", "cen_end", "cen_len", "project", "filterflag", "y")
print(head(cen_info,10))

##load HiCAT best layer file##
chr13_hicat_file <- "chr13.HiCAT.horstv.bed"
chr21_hicat_file <- "chr21.HiCAT.horstv.bed"
chr13_hicat_data <- read.csv(chr13_hicat_file, header=FALSE, sep="\t")
chr21_hicat_data <- read.csv(chr21_hicat_file, header=FALSE, sep="\t")
colnames(chr13_hicat_data) <- c("sample_hap_chrom", "horstart", "horend", "hor_mn", "horclass", "horclass_color", "horstv_index", "horstv_color")
colnames(chr21_hicat_data) <- c("sample_hap_chrom", "horstart", "horend", "hor_mn", "horclass", "horclass_color", "horstv_index", "horstv_color")
combined_data <- rbind(chr13_hicat_data, chr21_hicat_data) %>%
  mutate(sample_hap_chrom = if_else(sample_hap_chrom %in% target_chrs, 
                                    paste0("CHM13#", sample_hap_chrom), 
                                    sample_hap_chrom))
head(combined_data, 10)





color_map <- c(
  "#dfc27d" = "#C5C8D9",
  "#B66E64" = "#1a9850",
  "#F25CA2" = "#8C3074",  # 示例：添加更多颜色映射
  "#4a37b6" = "#F25CA2",
  "#f46d43" = "#058187",
  "#446b4c" = "#F2B33D",
  "#46c4c4" = "#0C4DA7",
  "#877DD4" = "#01325E",
  "#A50161" = "#D9D2BF"
)



recolor_file <- "/share/home/zhanglab/user/sunyanqing/vscode_scripts/plot_scripts/all.horstv.hicat.update.color.txt"
recolordf  <- read.csv(recolor_file, header=TRUE, sep="\t")
recolordf <- recolordf %>%
  mutate(new_color = recode(new_color, !!!color_map))
colnames(recolordf) <- c("horstv_index", "horstv_newcolor")
print(head(recolordf,10))

stv <- combined_data  %>%
  inner_join(cen_info, by = c("sample_hap_chrom" = "chrom")) %>%
  left_join(recolordf, by = "horstv_index") %>%
  mutate(horstv_newcolor = ifelse(horstv_color == '-', "#525252",
                                  ifelse(!is.na(horstv_newcolor), horstv_newcolor, horstv_color))) 
head(stv)

newcolor <- stv %>%
  filter(grepl("_", hor_mn), grepl("C42", horstv_index)) %>%
  distinct(hor_mn, horstv_newcolor)
head(newcolor,10)


######plot HORmon#######
##load HORmon file##
chr13_graph_file <-"chr13.graph.hordecomposition.final.xls"
chr21_graph_file <-"chr21.graph.hordecomposition.final.xls"
chr13_graph_data <- read.csv(chr13_graph_file, header=FALSE, sep="\t")
colnames(chr13_graph_data) <- c("sample_hap_chrom", "horstart", "horend", "hor_mn", "horclass", "horclass_color",  "horstv_color")
chr21_graph_data <- read.csv(chr21_graph_file, header=FALSE, sep="\t")
colnames(chr21_graph_data) <- c("sample_hap_chrom", "horstart", "horend", "hor_mn", "horclass", "horclass_color",  "horstv_color")

combined_data <- rbind(chr13_graph_data, chr21_graph_data) %>%
  mutate(sample_hap_chrom = if_else(sample_hap_chrom %in% target_chrs, 
                                    paste0("CHM13#", sample_hap_chrom), 
                                    sample_hap_chrom))

head(combined_data, 10)

target_data <- combined_data %>% 
  filter(sample_hap_chrom %in% target_df$sample_hap_chrom) %>%
  left_join(newcolor, by = "hor_mn") %>% # Join on hor_mn
  mutate(horstv_color = ifelse(!is.na(horstv_newcolor), horstv_newcolor, horstv_color)) %>% # Update horstv_color
  select(-horstv_newcolor)
head(target_data, 10)
print(unique(target_data$sample_hap_chrom))
print(target_data %>% filter(sample_hap_chrom == "CN1#Mat#chr21") %>% head(10))

head(cen_info,10)
graphstv <- target_data %>%
  left_join(cen_info, by = c("sample_hap_chrom" = "chrom")) 
print(graphstv %>% filter(sample_hap_chrom == "CN1#Mat#chr21") %>% head(10))


print(colnames(graphstv))
print(head(graphstv, 10))
print(unique(graphstv$sample_hap))
print(length(unique(graphstv$sample_hap_chrom)))

stv_min_hor_start <- graphstv %>% 
  filter(horclass == 42) %>%
  group_by(sample_hap_chrom) %>%
  summarise(min_hor_start = min(horstart, na.rm = TRUE)) %>%
  ungroup()
head(stv_min_hor_start,10)

graph_plot_stv_df <- graphstv %>%
  filter(horstart >= cen_start & horend <= cen_end) %>%
  inner_join(stv_min_hor_start, by = c("sample_hap_chrom")) %>%
  mutate(
    modstart = if_else(min_hor_start < cen_start, horstart - cen_start + 1, horstart - min_hor_start + 1),
    modend = if_else(min_hor_start < cen_start, horend - cen_start + 1, horend - min_hor_start + 1)
  ) %>%
  filter(modstart >= 0 & modend >= 0)
head(graph_plot_stv_df,10)
print(graph_plot_stv_df %>% filter(sample_hap_chrom == "CN1#Mat#chr21") %>% head(10))
print(length(unique(graph_plot_stv_df$sample_hap_chrom)))

max_len <- max(graph_plot_stv_df$modend)
y_labels <- graph_plot_stv_df %>%distinct(sample_hap_chrom, y)


checkdf <- graph_plot_stv_df %>%
  filter(horstv_color == "-")
head(checkdf,10)

p3 <- ggplot() +
  geom_rect(
    data = graph_plot_stv_df,
    aes(xmin = modstart, xmax = modend, ymin = y - 0.3, ymax = y + 0.3, fill = horstv_color),
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
ggsave("Figures/chr21_examples.cenhap.pdf", plot = p3 , device = "pdf", width = 12, height = 2)

