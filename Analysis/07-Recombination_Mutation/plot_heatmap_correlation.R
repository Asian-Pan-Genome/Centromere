library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("reshape2")
library("ggcorrplot")
library("cowplot")

setwd("/share/home/zhanglab/user/sunyanqing/human/periphy/CHM13/IBS")
setwd("/share/home/zhanglab/user/sunyanqing/human/periphy/CHM13/IBS_EAS")
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# args <- commandArgs(trailingOnly = TRUE)
# pearson_r_file <- args[1]
# pearson_p_file <- args[2]
# outpdf <- args[3]

chr <- "chr8"
pearson_r_file <- paste0(chr, "/out_pearson.xls")
pearson_p_file <- paste0(chr, "/out_pearson_p.xls")
outpdf <- paste0(chr, "/", chr, "_correlation_heatmap_EAS.pdf")
cat(pearson_r_file, pearson_p_file, outpdf)

pearson_r_file <- paste0(chr, "/out_pearson_hor.xls")
pearson_p_file <- paste0(chr, "/out_pearson_p_hor.xls")
outpdf <- paste0(chr, "/", chr, "_correlation_heatmap_hor_update.pdf")
cat(pearson_r_file, pearson_p_file, outpdf)

# pearson_r_file <- "chr6/out_pearson_t99.xls"
# pearson_p_file <- "chr6/out_pearson_p_t99.xls"
# outpdf <- "chr6/chr6_correlation_heatmap_t99.pdf"


#load pearson_r, p-value matrix##
# pearson_r_file <- "chr17/out_pearson.xls" 
r_matrix <- read.csv(pearson_r_file, row.names = 1, header = TRUE, sep="\t", check.names = FALSE)
head(r_matrix)
head(r_matrix[1:10, 1:10])

# pearson_p_file <- "chr17/out_pearson_p.xls"
p_matrix <- read.csv(pearson_p_file, row.names = 1, header = TRUE, sep="\t", check.names = FALSE)
head(p_matrix[1:10, 1:10])

r_long <- r_matrix %>%
  rownames_to_column(var = "row") %>%
  gather(key = "col", value = "r", -row)
head(r_long)
print(r_long$row) 
dim(r_long)
p_long <- p_matrix %>%
  rownames_to_column(var = "row") %>%
  gather(key = "col", value = "p", -row)
head(p_long)
print(p_long$row)
dim(p_long)

merged_long <- merge(r_long, p_long, by = c("row", "col"))
head(merged_long)

print(merged_long$row)
# print(merged_long %>% filter(row == 2 & col == 10))
# print(merged_long %>% filter(row == 10 & col == 2))

filtered_long <- merged_long %>%
  # filter(p < 0.05) %>%
  select(row, col, r, p) %>%
  mutate(
    r = if_else(row == col, NA_real_, r),
    r = if_else(r < 0, 0, r)
  )

head(filtered_long)
print(filtered_long$row)
print(filtered_long$col)

# filtered_long <- filtered_long %>%
#   filter(row < 14)
# selected_rows <- c('42', '43', '44', '45', 'CEN', '47', '48', '49', '50', '51','52', '54', '56', '57', '58', '59', '60')
# filtered_long <- filtered_long %>%
#   filter(row %in% selected_rows) %>%
#   filter(col %in% selected_rows)
dim(filtered_long)

filtered_long$row <- factor(filtered_long$row, levels = rownames(r_matrix), order=TRUE)
filtered_long$col <- factor(filtered_long$col, levels = colnames(r_matrix), order=TRUE)


filtered_long <- filtered_long %>%
  filter(as.numeric(row) >= as.numeric(col) | is.na(as.numeric(row)) | is.na(as.numeric(col)))

head(filtered_long,30)

print(filtered_long$row)
print(filtered_long$col)

###ggplot2##

# p1 <- ggplot(filtered_long, aes(x = col, y = row, fill = r)) +
#   geom_tile() +
#   scale_fill_gradient2(low = "#4575b4", mid = "#fee090", high = "#d73027", midpoint = 0.5, name = "R", limits = c(0, 1), na.value = "white") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=4),
#         axis.text.y = element_text(size =4),
#         legend.position = "none") +
#   labs(x = "region_index", y = "region_index")
# print(p1)
# # outpdf <- "chr17/correlation_heatmap.pdf"
# ggsave(outpdf, plot = p1 , device = "pdf", width = 5, height =5)


##change to Spectral palette##
# p2 <- ggplot(filtered_long, aes(x = col, y = row, fill = r)) +
#   geom_tile() +
#   scale_fill_distiller(palette = "Spectral", direction = -1, name = "R", limits = c(0, 1), na.value = "white") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=4),
#         axis.text.y = element_text(size =4),
#         legend.position = "none") +
#   labs(x = "region_index", y = "region_index")
# print(p2)
# ggsave(outpdf, plot = p2 , device = "pdf", width = 5, height =5)


##change to Spectral palette with cumtomized limit and breaks##
library(RColorBrewer)

p3 <- ggplot(filtered_long, aes(x = col, y = row, fill = r)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = rev(RColorBrewer::brewer.pal(11, "Spectral")), # 使用 Spectral 调色板
    values = c(0, 0.5, 1), # 指定颜色位置：0 对应低值，0.5 对应中间值，1 对应高值
    limits = c(0, 1), # 数据范围
    na.value = "white" # 缺失值颜色
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4),
        axis.text.y = element_text(size = 4),
        legend.position = "none") +
  labs(x = "region_index", y = "region_index")
print(p3)
ggsave(outpdf, plot = p3 , device = "pdf", width = 5, height =5)
