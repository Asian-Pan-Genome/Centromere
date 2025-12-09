setwd("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap")
rm(list = ls())

library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("patchwork")
library("pheatmap")
library(ggnewscale)
library(RColorBrewer)
library(grid)

# args <- commandArgs(trailingOnly = TRUE)
# chr <- args[1]

chr <- "chr16"
distype <- "hor"
distype <- "merge"

save_pheatmap_pdf <- function(pheatmap, filename, width = 10, height = 10) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(pheatmap$gtable)
  dev.off()
}

get_clusters <- function(hclust_obj, n_clusters) {
  clusters <- cutree(hclust_obj, k = n_clusters)
  return(clusters)
}


# Load the distance matrix #
distance_file <- paste0("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/distance/", distype ,"_distance_",chr,".csv")
distance_matrix <- read.csv(distance_file, row.names = 1, header = TRUE)
colnames(distance_matrix) <- rownames(distance_matrix)
head(distance_matrix[1:10, 1:10])
print(rownames(distance_matrix))
dim(distance_matrix)


# Load the cent info #           
df_cen <- read.csv("cent_chrom.txt", header=TRUE, sep="\t")
df_cen <- df_cen %>% filter(chrom == chr)
head(df_cen,10)
complete_df_cen <- df_cen %>% filter(filterflag == 0)
print(dim(complete_df_cen))


df_pop <- read.csv("populaion.xls", header=TRUE, sep="\t")
df_cen <- df_cen  %>% left_join(df_pop %>% select(sample, superpopulation), by = "sample") %>%
    mutate(sample_hap = str_replace(sample_hap, "CHM13_-", "CHM13"))
head(df_cen,10)
dim(df_cen)

selected_sample_hap_chrom <- df_cen %>% filter(filterflag == 0) %>%
    mutate(sample_hap = str_replace(sample_hap, "CHM13_-", "CHM13"))
head(selected_sample_hap_chrom,10)
print(length(unique(selected_sample_hap_chrom$sample_hap_chrom)))

# select complete centromere for distance matrix #
filtered_names <- selected_sample_hap_chrom$sample_hap
filtered_names <- filtered_names[filtered_names %in% rownames(distance_matrix) & filtered_names %in% colnames(distance_matrix)]

if (chr== "chr9"){
  filter_rows <- c("HG00514_hap2", 
                   "C019-CHA-E19_Mat", 
                   "C051-CHA-N11_Mat",
                   "HG00864_hap2", 
                   "HG02011_hap2",
                   "HG01890_hap2",
                   "HG02011_hap1"
                   )
  filtered_names <-filtered_names[!filtered_names %in% filter_rows]
}
print(filtered_names)
print("HG00514_hap2" %in% filtered_names)
ff_selected_to_filtered <- setdiff(selected_sample_hap_chrom$sample_hap, filtered_names)
print(ff_selected_to_filtered)
filtered_distance_matrix <- distance_matrix[filtered_names, filtered_names]
head(filtered_distance_matrix[1:10, 1:10])

filtered_distance_matrix <- as.matrix(filtered_distance_matrix)
max_value <- max(filtered_distance_matrix, na.rm = TRUE)

##heatmap##
annotation  <- df_cen %>% filter(sample_hap %in% filtered_names) %>%
    select(sample_hap, superpopulation) %>%
    mutate(superpopulation = ifelse(superpopulation == "EAS-APG", "EAS", superpopulation)) %>%
    column_to_rownames("sample_hap")
superpopulation_colors <- c("AFR" = "#319b62", "AMR" = "#939393", "EAS" = "#d22e77", "EUR" = "#0070c0", "SAS" = "#893f8b")
head(annotation,10)
dim(annotation)

custom_colors <- colorRampPalette(brewer.pal(11, "Spectral"))(50)

# Check for NA, NaN, or Inf values
any(is.na(filtered_distance_matrix))  # Check for NA
any(is.nan(filtered_distance_matrix)) # Check for NaN
any(is.infinite(filtered_distance_matrix)) # Check for Inf


##Cenhap cluster##
cenhap_file <- paste0("cenhap_classification/", chr, "_cenhap.xls")
print(cenhap_file)


p1 <- pheatmap(filtered_distance_matrix, 
                 cluster_rows = TRUE, 
                 cluster_cols = TRUE,
                 fontsize_row = 1,
                 fontsize_col = 0.5,
                 color = custom_colors,
                 breaks = seq(0, max_value, length.out = length(custom_colors) + 1),
                 annotation_col = annotation,
                 annotation_colors = list(superpopulation = superpopulation_colors),
                 legend = TRUE,
                 annotation_names_col = FALSE,
                 clustering_method = "average")
print(p1)
out_p1 <- paste0("Figures/",chr, ".distance_", distype, "_matrix.heatmap.updated.pdf")
print(out_p1)
save_pheatmap_pdf(p1, out_p1, width = 10, height = 8)

##save clustered distance matrix##
clustered_distance_matrix <- filtered_distance_matrix[p1$tree_row$order, p1$tree_col$order]
out_clustered_matrix <- paste0("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/Clustering/", chr, ".", distype , ".distance_clustered.updated.csv")
write.csv(clustered_distance_matrix, out_clustered_matrix, row.names = TRUE, quote = FALSE)

##save clustered result##
clustered_rows <- data.frame(Row = rownames(filtered_distance_matrix))
for (n in 2:10) {
  clusters <- get_clusters(p1$tree_row, n)
  clustered_rows <- cbind(clustered_rows, clusters)
  colnames(clustered_rows)[n] <- paste0("Cluster_n", n)
}
clustered_rows <- clustered_rows[p1$tree_row$order, ] %>% select(-Row)
print(head(clustered_rows))
out_combined_clusters <- paste0("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/Clustering/", chr,".", distype , ".distance_clustered.result.updated.csv")
write.csv(clustered_rows, out_combined_clusters, row.names = TRUE, quote = FALSE)


##HORhap##
orderdf <- rownames_to_column(clustered_rows, var = "sample_hap") %>%
  mutate(y = n() - row_number()+1) 
rownames(orderdf) <- NULL
head(orderdf,10)
print(orderdf %>% filter(sample_hap =="CHM13"))



##HORstv##
infile<-paste0(chr, ".HiCAT.horstv.bed")
data <- read.csv(infile, header=FALSE, sep="\t")
colnames(data) <- c("sample_hap_chrom", "horstart", "horend", "hor_mn", "horclass", "horclass_color", "horstv_index", "horstv_color")
print(data[20000, ])

##load HORstv recolor file##
recolor_file <- "/share/home/zhanglab/user/sunyanqing/vscode_scripts/plot_scripts/all.horstv.hicat.update.color.txt"
recolordf  <- read.csv(recolor_file, header=TRUE, sep="\t")
colnames(recolordf) <- c("horstv_index", "horstv_newcolor")
print(head(recolordf,10))

stv <- data %>%
  mutate(
    sample_hap = case_when(
      sample_hap_chrom == chr ~ "CHM13",
      str_detect(sample_hap_chrom, "#") ~ sample_hap_chrom %>%
        str_replace("#", "_") %>%
        str_replace(paste0("#", chr), ""),
      str_detect(sample_hap_chrom, "_") ~ str_replace(sample_hap_chrom, paste0("_", chr), ""),
      TRUE ~ NA_character_ # 默认情况
    )
  ) %>%
  filter(sample_hap %in% filtered_names) %>%
  left_join(recolordf, by = "horstv_index") %>%
  mutate(horstv_newcolor = ifelse(horstv_color == '-', "#525252", 
                                  ifelse(!is.na(horstv_newcolor), horstv_newcolor, horstv_color))) %>%
  select(sample_hap_chrom, sample_hap, horstart, horend, hor_mn, horstv_index, horstv_newcolor)
print(head(stv, 10))
print(unique(stv$sample_hap))
print(length(unique(stv$sample_hap_chrom)))

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
cen_info <- df_cen %>% 
    filter(!sample_hap_chrom %in% filter_chrs) %>% 
  select(sample_hap_chrom,  start, end, len, project, filterflag) %>%
  mutate(sample_hap_chrom = str_replace(sample_hap_chrom, "CHM13#", "")) 

colnames(cen_info) <- c("sample_hap_chrom", "cen_start", "cen_end", "cen_len", "project", "filterflag")

print(head(cen_info,10))


stv_min_hor_start <- stv %>%
  group_by(sample_hap_chrom) %>%
  summarise(min_hor_start = min(horstart, na.rm = TRUE)) %>%
  ungroup()
head(stv_min_hor_start,10)
print(stv_min_hor_start%>%filter(sample_hap_chrom == "chr2"))

print(cen_info %>% filter(sample_hap_chrom == "chr2"))

plot_stv_df <- stv %>%
  inner_join(cen_info, by = c("sample_hap_chrom")) %>%
  filter(horstart >= cen_start & horend <= cen_end) %>%
  inner_join(stv_min_hor_start, by = c("sample_hap_chrom")) %>%
  mutate(
    modstart = if_else(min_hor_start < cen_start, horstart - cen_start + 1, horstart - min_hor_start + 1),
    modend = if_else(min_hor_start < cen_start, horend - cen_start + 1, horend - min_hor_start + 1)
  ) %>%
  inner_join(orderdf, by = c("sample_hap"))
head(plot_stv_df,10)

print(unique(plot_stv_df$horstv_newcolor))
print(length(unique(plot_stv_df$sample_hap_chrom)))
# ff_stv <- setdiff( unique(stv$chrom), unique(plot_stv_df$chrom))
# print(ff_stv)

max_len <- max(plot_stv_df$modend)
y_labels <- plot_stv_df %>%distinct(sample_hap_chrom, y)

p2 <- ggplot() +
  geom_rect(
    data = plot_stv_df,
    aes(xmin = modstart, xmax = modend, ymin = y - 0.3, ymax = y + 0.3, fill = horstv_newcolor),
    color = NA 
  ) +
  scale_fill_identity() + 
  # new_scale("fill") +  
  # geom_rect( 
  #     data = plot_stv_df %>% distinct(chrom, y, modend),
  #     aes(xmin = 0, xmax = modend,  ymin = y - 0.3, ymax = y + 0.3),
  #     fill = NA, 
  #     color = "black", 
  #     size = 0.1
  #   ) +
  # new_scale("fill") +  
  # geom_point(
  #   data = target_te,
  #   aes(
  #     x = (modstart + modend) / 2, 
  #     y = y,  
  #     fill = motif,    
  #     color = motif 
  #   ),
  #   shape = 25,  
  #   size = 0.2,
  #   stroke = 0   
  # ) +
  # scale_color_manual(values = motif_color) + 
  # scale_fill_manual(values = motif_color) +   
  scale_x_continuous(
      breaks = seq(0, max_len, by = 2000000), 
      labels = function(x) x / 1e6 #Mb
    ) +  
  scale_y_continuous(breaks = y_labels$y, labels = y_labels$sample_hap_chrom) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 1),
    axis.title.y = element_text(size = 5)
  )

print(p2)


##graph##
graph_file <-paste0(chr, ".graph.hordecomposition.final.xls")
print(graph_file)
data <- read.csv(graph_file, header=FALSE, sep="\t")
colnames(data) <- c("sample_hap_chrom", "horstart", "horend", "hor_mn", "horclass", "horclass_color",  "horstv_color")
head(data)
print(data[30000, ])

head(stv)
newcolor <- stv %>%
  filter(grepl("_", hor_mn), grepl("C42", horstv_index)) %>%
  distinct(hor_mn, horstv_newcolor)
head(newcolor,10)


data <- data %>%
  left_join(newcolor, by = "hor_mn") %>% # Join on hor_mn
  mutate(horstv_color = ifelse(!is.na(horstv_newcolor), horstv_newcolor, horstv_color)) %>% # Update horstv_color
  select(-horstv_newcolor)
head(data, 10)

graphstv <- data %>%
  mutate(
    sample_hap = case_when(
      sample_hap_chrom == chr ~ "CHM13",
      str_detect(sample_hap_chrom, "#") ~ sample_hap_chrom %>%
        str_replace("#", "_") %>%
        str_replace(paste0("#", chr), ""),
      str_detect(sample_hap_chrom, "_") ~ str_replace(sample_hap_chrom, paste0("_", chr), ""),
      TRUE ~ NA_character_ # 默认情况
    )
  ) %>%
  filter(sample_hap %in% filtered_names) 
print(head(stv, 10))
print(unique(stv$sample_hap))
print(length(unique(stv$sample_hap_chrom)))

stv_min_hor_start <- graphstv %>%
  group_by(sample_hap_chrom) %>%
  summarise(min_hor_start = min(horstart, na.rm = TRUE)) %>%
  ungroup()
head(stv_min_hor_start,10)



graph_plot_stv_df <- graphstv %>%
  inner_join(cen_info, by = c("sample_hap_chrom")) %>%
  filter(horstart >= cen_start & horend <= cen_end) %>%
  inner_join(stv_min_hor_start, by = c("sample_hap_chrom")) %>%
  mutate(
    modstart = if_else(min_hor_start < cen_start, horstart - cen_start + 1, horstart - min_hor_start + 1),
    modend = if_else(min_hor_start < cen_start, horend - cen_start + 1, horend - min_hor_start + 1)
  ) %>%
  inner_join(orderdf, by = c("sample_hap"))
head(graph_plot_stv_df,10)
print(length(unique(graph_plot_stv_df$sample_hap_chrom)))

max_len <- max(graph_plot_stv_df$modend)
y_labels <- graph_plot_stv_df %>%distinct(sample_hap_chrom, y)


p3 <- ggplot() +
  geom_rect(
    data = graph_plot_stv_df,
    aes(xmin = modstart, xmax = modend, ymin = y - 0.3, ymax = y + 0.3, fill = horstv_color),
    color = NA 
  ) +
  scale_fill_identity() + 
  scale_x_continuous(
      breaks = seq(0, max_len, by = 2000000), 
      labels = function(x) x / 1e6 #Mb
    ) +  
  scale_y_continuous(breaks = y_labels$y, labels = y_labels$sample_hap_chrom) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 1),
    axis.title.y = element_text(size = 5)
  )

print(p3)



# library(gridExtra)
# pp1 <- p1$gtable
# p  <- grid.arrange(topp_grob, p1, p3, ncol = 3)
# p <- grid.arrange(pp1, p2, p3, ncol = 3, widths = c(2,1,1))
# outpdf <- paste0("Figures/", chr, ".", distype, ".distance.updated.cenhap.png")
# ggsave(outpdf, plot = p , device = "png", width = 12, height = 6, dpi=720)

library(gridExtra)
outp <- grid.arrange(p2, p3, ncol = 2, widths = c(1,1))
outpdf <- paste0("Figures/", chr, ".", distype, ".distance.updated.cenhap.pdf")
print(outpdf)
ggsave(outpdf, plot = outp , device = "pdf", width = 6, height = 6)
