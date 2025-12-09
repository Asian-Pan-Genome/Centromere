setwd("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap")
rm(list = ls())

library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("patchwork")
library(grid)
library(ComplexHeatmap)
library(cowplot)
library("ggplotify")
library(circlize)
library(RColorBrewer)

chr <- "chrY"

##load cent_length file##
df_cen <- read.csv("cent_chrom.txt", header=TRUE, sep="\t")
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
df_cen <- df_cen %>% 
    filter(!sample_hap_chrom %in% filter_chrs) %>% 
  select(sample_hap_chrom,  start, end, len, project, filterflag)
colnames(df_cen) <- c("chrom", "cen_start", "cen_end", "cen_len", "project", "filterflag")
print(head(df_cen,10))

##load cenanno file##
cenanno_file <- paste0(chr, ".merged.cenanno.bed")
cenanno  <- read.csv(cenanno_file, header=FALSE, sep="\t")
colnames(cenanno) <- c("chrom", "start", "end", "satellite", "score", "strand", "s", "e", "color")
cenanno <- cenanno %>% mutate(color = ifelse(color == "", "#891640", color))
print(head(cenanno,10))

##integrate data##
plot_cenanno_df <- cenanno %>% 
  inner_join(df_cen, by = c("chrom")) %>%
  filter(start >= cen_start & end <= cen_end) %>%
  mutate(modstart = start - cen_start + 1 ) %>%
  mutate(modend = end - cen_start + 1 )
head(plot_cenanno_df,10)


##summary satellite blocks number##
cenanno_summary <- plot_cenanno_df %>% filter(satellite != "rDNA", end-start >= 10000, filterflag == 0) %>% 
  group_by(chrom, satellite) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = satellite, values_from = count, values_fill = list(count = 0))
head(cenanno_summary,10)

##heatmap for composition of continuous satellite blocks ##
cenanno_matrix <- as.matrix(cenanno_summary %>% column_to_rownames(var = "chrom"))

ht <- Heatmap(
  cenanno_matrix, 
  name = "Satellite Count",
  row_names_gp = gpar(fontsize = 3), 
  heatmap_legend_param = list(
    legend_position = "topright" 
  )
)

print(max(cenanno_matrix))
max_value <- max(cenanno_matrix)
min_value <- min(cenanno_matrix)
breaks <- seq(min_value, max_value, length.out = 9)
colors <- rev(brewer.pal(9, "Spectral"))
# colors[1] <- "#abd9e9"
# print(colors)

ht <- Heatmap(
  cenanno_matrix, 
  name = "Satellite Count",
  # rect_gp = gpar(col = "white", lwd = 1),
  col = colorRamp2(breaks, colors),
  # row_km = 7,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  heatmap_legend_param = list(
    legend_position = "topright"
  ))
##draw heatmap and get row order##
ht_drawn <- draw(ht)
# pdf("chr4_satellite_organization_heatmap.pdf", width = 6, height = 6)
# draw(ht)
# dev.off()

# row_order <- row_order(ht_drawn)
# inverted_row_order <- unlist(lapply(row_order, rev))
# cenanno_matrix_inverted <- cenanno_matrix[inverted_row_order, ]

# ht_inverted <- Heatmap(
#   cenanno_matrix_inverted, 
#   name = "Satellite Count",
#   col = colorRamp2(breaks, colors),
#   show_column_dend = FALSE,
#   show_row_names = FALSE,
#   column_names_rot = 90,  # Rotate column labels to be horizontal
#   heatmap_legend_param = list(
#     legend_position = "topright"
#   )
# )
# pdf("chr4_satellite_organization_heatmap.pdf", width = 6, height = 6)
# draw(ht_inverted)
# dev.off()
row_order <- row_order(ht_drawn)
print(row_order)
ordered_row_labels <- rownames(cenanno_matrix)[row_order]
print(ordered_row_labels)

ordered_row_df <- data.frame(
  chrom = ordered_row_labels,
  y = length(ordered_row_labels) - seq_along(ordered_row_labels) + 1
)
head(ordered_row_df,10)

##2nd cluster##
head(cenanno_summary)

##change column names##
cenanno_summary2 <- cenanno_summary
old_cols <- colnames(cenanno_summary2)
new_cols <- sapply(old_cols, function(col) {
  if (col == "chrom") {
    col
  } else {
    paste0(col, "_count")
  }
})
print(new_cols)
colnames(cenanno_summary2) <- new_cols
head(cenanno_summary2)

#merge cenanno_summary2 and plot_cenanno_df
satellite_marker <- c("ASat"=1, "HSat1A"=2, "HSat1B"=3, "HSat2"=4, "HSat3"=5, "BSat"=6, "GSat"=7)
head(plot_cenanno_df)
df2 <- plot_cenanno_df %>% inner_join(cenanno_summary2, by = c("chrom")) %>%
  filter(end-start >= 10000) %>%
  mutate(marker = satellite_marker[satellite])
head(df2,20)

subdf <- df2 %>% filter(ASat_count == 8)
head(subdf,10)
wide_subdf <- subdf %>%
  group_by(chrom) %>%
  mutate(row = row_number()) %>%
  select(chrom, row, marker) %>%
  spread(key = row, value = marker)

head(wide_subdf,10)

## merge cenanno_summary2 and df_cen ##
cenanno_summary2 <- cenanno_summary2 %>% 
  inner_join(df_cen, by = c("chrom"))
head(cenanno_summary2,10)

##plot different group lenth variation##
plot_apg_cenanno_summary2 <- cenanno_summary2 %>% filter(project == "APG")
p2 <- ggplot(plot_apg_cenanno_summary2, aes(x = factor(ASat_count), y = cen_len)) +
    geom_boxplot() +
    coord_flip()
print(p2)

##add order info##
plot_cenanno_df <- plot_cenanno_df %>% 
  inner_join(ordered_row_df, by = c("chrom")) 
head(plot_cenanno_df,10)

##plot##
p1 <- ggplot() +
    geom_rect(
      data = plot_cenanno_df,
      aes(xmin = modstart, xmax = modend, ymin = y - 0.4, ymax = y + 0.4, fill = color),
      color = NA 
    ) +
    scale_fill_identity() +
    # geom_rect( 
    #   data = combined_df %>% distinct(chrom, y, cen_len),
    #   aes(
    #     xmin = 0, 
    #     xmax = cen_len, 
    #     ymin = y - 0.4,
    #     ymax = y + 0.4
    #   ),
    #   fill = NA, 
    #   color = "black", 
    #   size = 1
    # ) +
    scale_x_continuous(
      breaks = seq(0, max(plot_cenanno_df$modend), by = 5000000), 
      labels = function(x) x / 1e6 #Mb
    ) +
    scale_y_continuous(
      breaks = unique(plot_cenanno_df$y), 
      labels = unique(plot_cenanno_df$chrom)
    ) +
    coord_cartesian(
      ylim = c(min(plot_cenanno_df$y) - 0.4, max(plot_cenanno_df$y) + 0.4)
    ) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 3),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank())
print(p1)

##merge two plots##
heatmap_ggplot <- grid.grabExpr(draw(ht))
# heatmap_ggplot <- ggplotify::as.ggplot(heatmap_ggplot) + theme(aspect.ratio = 1)  # Flatten the heatmap
combined_plot <- plot_grid(heatmap_ggplot, p1, ncol = 2, rel_widths = c(1, 1))  # Combine plots vertically
print(combined_plot)
ggsave(paste0("satellite_organization/", chr, "_satellite_organization_plot.pdf"), combined_plot, width = 10, height = 16)

##save output##
reordered_matrix <- cenanno_matrix[row_order, ]
write.csv(reordered_matrix, file = paste0("satellite_organization/", chr, "_reordered_matrix.csv"), row.names = TRUE)
print(reordered_matrix)



####### after manual check and plot ########
rm(list = ls())
data <- read.csv("satellite_transition_20250308.xls", sep="\t", header = TRUE)
data$chr <- factor(data$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY"), ordered = TRUE)


stat <- data %>%
  group_by(chr, type) %>%
  summarise(count = n()) %>%
  arrange(chr, desc(count)) %>%
  mutate(index = row_number())%>%
  mutate(level = case_when(
    index %in% 1:10 ~ paste0("top", index),
    TRUE ~ "others"
  )) 
head(stat,10)


color_values <- c(
  "top1" = "#1a3951",
  "top2" = "#087D6F",
  "top3" = "#15959F",
  "top4" = "#F1E4B3",
  "top5" = "#FFC753",
  "top6" = "#FFB1A8",
  "top7" = "#EC9770",
  "top8" = "#C7402D",
  "top9" = "#D9430D",
  "top10" = "#A60A0A",
  "others" = "grey"
)

stat$level <- factor(stat$level, levels = rev(c(paste0("top", 1:10), "others")), ordered = TRUE)
stat <- stat %>% arrange(chr, level)


p <- ggplot(stat, aes(x = chr, y = count, fill = level)) +
  geom_bar(position = 'stack', stat = "identity") +
  theme_classic() +
  labs(x = "", y = "The number of centromeres") +
  scale_fill_manual(values = color_values) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black")
  )
print(p)
outpdf <- paste0("satellite_organization/satellite_transition_nolegend.pdf")
ggsave(outpdf, plot = p , device = "pdf", width = 10, height = 5)

### plot examples of chr4 ###
target_chrs <- read.csv("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/chr4_duplication/example_chr4_ids.xls", header = FALSE)
target_chrs <- read.csv("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/chr4_duplication/example_chr4_ids_v2.xls", header = FALSE)
colnames(target_chrs) <- c("chrom")
target_chrs <- target_chrs %>% mutate(y = 2*(n() - row_number())+ 1) 
head(target_chrs,20)
head(plot_cenanno_df,20)

target_cenanno_df <- plot_cenanno_df %>% filter(chrom %in% target_chrs$chr) %>% inner_join(target_chrs, by = c("chrom")) 
# target_cenanno_df <- target_cenanno_df%>% 
#   mutate(color = case_when(
#     color == "#77A0BA" ~ "#00B0F0",
#     color == "#891640" ~ "#C91D32",
#     color == "#18988B" ~ "#248E87",
#     TRUE ~ color
#   ))
head(target_cenanno_df)


##load te data##
tefile <- paste0(chr, ".te.totalform.xls")
tedf <- read.csv(tefile, header = FALSE, sep = "\t")
colnames(tedf) <- c("project", "superpopulation", "population", "sample", "hap", "chr", "start", "end", "length", "strand", "motif", "class", "family", "flag")
colnames(tedf) <- make.names(colnames(tedf), unique = TRUE)
head(tedf,10)

plot_te_file <- tedf %>% mutate(chrom = paste(sample, hap, chr, sep = "#")) %>% filter(chrom %in% target_chrs$chrom, flag != NA) 
head(plot_te_file,10)
te_num <- plot_te_file %>% group_by(chrom) %>% summarise(count = n())
print(te_num)

p1 <- ggplot() +
    geom_rect(
      data = target_cenanno_df,
      aes(xmin = modstart, xmax = modend, ymin = y - 0.4, ymax = y + 0.4, fill = color),
      color = NA 
    ) +
    scale_fill_identity() +
    # geom_rect( 
    #   data = combined_df %>% distinct(chrom, y, cen_len),
    #   aes(
    #     xmin = 0, 
    #     xmax = cen_len, 
    #     ymin = y - 0.4,
    #     ymax = y + 0.4
    #   ),
    #   fill = NA, 
    #   color = "black", 
    #   size = 1
    # ) +
    scale_x_continuous(
      breaks = seq(0, max(target_cenanno_df$modend), by = 2000000), 
      labels = function(x) x / 1e6 #Mb
    ) +
    scale_y_continuous(
      breaks = unique(target_cenanno_df$y), 
      labels = unique(target_cenanno_df$chrom)
    ) +
    coord_cartesian(
      ylim = c(min(target_cenanno_df$y) - 0.4, max(target_cenanno_df$y) + 0.4)
    ) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 3),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank())
print(p1)

outpdf <- paste0("satellite_organization/chr4_examples_barplot_v2.pdf")
ggsave(outpdf, plot = p1 , device = "pdf", width = 10, height = 8)


#for single_barplot##
for (ichr in unique(target_cenanno_df$chrom)) {
  chr_data <- target_cenanno_df %>% filter(chrom == ichr)
  print(head(chr_data,10))
  
  p1 <- ggplot() +
    geom_rect(
      data = chr_data,
      aes(xmin = modstart, xmax = modend, ymin = -0.3, ymax = 0.3, fill = color),
      color = NA 
    ) +
    scale_fill_identity() +
    # scale_x_continuous(
    #   breaks = seq(0, max(target_cenanno_df$modend), by = 2000000), 
    #   labels = function(x) x / 1e6 #Mb
    # ) +
    # scale_y_continuous(
    #   breaks = 0, 
    #   labels = ichr
    # ) +
    coord_cartesian(
      ylim = c(-1,1)
    ) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank())
  print(p1)
  outpdf <- paste0("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/chr4_duplication/barplot/", ichr, "_barplot.pdf")
  ggsave(outpdf, plot = p1, device = "pdf", width = 10, height = 1)
  print(paste(ichr, "finished"))
}



### plot examples of chr1 ###
target_chrs <- read.csv("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/chr1_duplication_inversion/example_chr1_ids.xls", header = FALSE)
colnames(target_chrs) <- c("chrom")
target_chrs <- target_chrs %>% mutate(y = 2*(n() - row_number())+ 1) 
head(target_chrs,20)
head(plot_cenanno_df,20)

target_cenanno_df <- plot_cenanno_df %>% filter(chrom %in% target_chrs$chr) %>% inner_join(target_chrs, by = c("chrom")) 
# target_cenanno_df <- target_cenanno_df%>% 
#   mutate(color = case_when(
#     color == "#77A0BA" ~ "#00B0F0",
#     color == "#891640" ~ "#C91D32",
#     color == "#18988B" ~ "#248E87",
#     TRUE ~ color
#   ))
head(target_cenanno_df)


##load te data##
tefile <- paste0(chr, ".te.totalform.xls")
tedf <- read.csv(tefile, header = FALSE, sep = "\t")
colnames(tedf) <- c("project", "superpopulation", "population", "sample", "hap", "chr", "start", "end", "length", "strand", "motif", "class", "family", "flag")
colnames(tedf) <- make.names(colnames(tedf), unique = TRUE)
head(tedf,10)

plot_te_file <- tedf %>% mutate(chrom = paste(sample, hap, chr, sep = "#")) %>% filter(chrom %in% target_chrs$chrom, flag != NA) 
head(plot_te_file,10)
te_num <- plot_te_file %>% group_by(chrom) %>% summarise(count = n())
print(te_num)

p1 <- ggplot() +
    geom_rect(
      data = target_cenanno_df,
      aes(xmin = modstart, xmax = modend, ymin = y - 0.4, ymax = y + 0.4, fill = color),
      color = NA 
    ) +
    scale_fill_identity() +
    # geom_rect( 
    #   data = combined_df %>% distinct(chrom, y, cen_len),
    #   aes(
    #     xmin = 0, 
    #     xmax = cen_len, 
    #     ymin = y - 0.4,
    #     ymax = y + 0.4
    #   ),
    #   fill = NA, 
    #   color = "black", 
    #   size = 1
    # ) +
    scale_x_continuous(
      breaks = seq(0, max(target_cenanno_df$modend), by = 5000000), 
      labels = function(x) x / 1e6 #Mb
    ) +
    scale_y_continuous(
      breaks = unique(target_cenanno_df$y), 
      labels = unique(target_cenanno_df$chrom)
    ) +
    coord_cartesian(
      ylim = c(min(target_cenanno_df$y) - 0.4, max(target_cenanno_df$y) + 0.4)
    ) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 3),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank())
print(p1)

outpdf <- paste0("satellite_organization/chr1_examples_barplot.pdf")
ggsave(outpdf, plot = p1 , device = "pdf", width = 10, height = 8)


#for single_barplot##
for (ichr in unique(target_cenanno_df$chrom)) {
  chr_data <- target_cenanno_df %>% filter(chrom == ichr)
  print(head(chr_data,10))
  
  p1 <- ggplot() +
    geom_rect(
      data = chr_data,
      aes(xmin = modstart, xmax = modend, ymin = -0.3, ymax = 0.3, fill = color),
      color = NA 
    ) +
    scale_fill_identity() +
    # scale_x_continuous(
    #   breaks = seq(0, max(target_cenanno_df$modend), by = 2000000), 
    #   labels = function(x) x / 1e6 #Mb
    # ) +
    # scale_y_continuous(
    #   breaks = 0, 
    #   labels = ichr
    # ) +
    coord_cartesian(
      ylim = c(-1,1)
    ) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank())
  print(p1)
  outpdf <- paste0("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/chr1_duplication_inversion/", ichr, "_barplot.pdf")
  ggsave(outpdf, plot = p1, device = "pdf", width = 10, height = 1)
  print(paste(ichr, "finished"))
}

### plot all chromosome ###
head(plot_cenanno_df,10)

plot_cenanno_df <- plot_cenanno_df %>% arrange(filterflag)
head(plot_cenanno_df,10)

unique_chrom_df <- data.frame(chrom = unique(plot_cenanno_df$chrom))
head(unique_chrom_df,10)
unique_chrom_df <- unique_chrom_df %>% mutate(y = (n() - row_number())+ 1)
head(unique_chrom_df,10)
plot_cenanno_df <- plot_cenanno_df %>% inner_join(unique_chrom_df, by = c("chrom"))
head(plot_cenanno_df,10)

check <- plot_cenanno_df %>% filter(filterflag == 1)
head(check,10)

p1 <- ggplot() +
    geom_rect(
      data = plot_cenanno_df,
      aes(xmin = modstart, xmax = modend, ymin = y - 0.4, ymax = y + 0.4, fill = color),
      color = NA 
    ) +
    scale_fill_identity() +
    scale_x_continuous(
      breaks = seq(0, max(plot_cenanno_df$modend), by = 5000000), 
      labels = function(x) x / 1e6 #Mb
    ) +
    scale_y_continuous(
      breaks = unique(plot_cenanno_df$y), 
      labels = unique(plot_cenanno_df$chrom)
    ) +
    coord_cartesian(
      ylim = c(min(plot_cenanno_df$y) - 0.4, max(plot_cenanno_df$y) + 0.4)
    ) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 2),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank())
print(p1)
outpdf <- paste0("satellite_organization/", chr, ".all_chromosome_barplot.pdf")
print(outpdf)
ggsave(outpdf, plot = p1, device = "pdf", width = 15, height = 15)
