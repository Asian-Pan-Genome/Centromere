setwd("/share/home/zhanglab/user/sunyanqing/human/periphy/CHM13/graph-diploid/chr1")
rm(list = ls())

library(ggplot2)
library(reshape2)
library("dplyr")
library("tidyr")
library("tidyverse")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Please provide input file and output file as command line arguments.")
}
input_file <- args[1]
output_file <- args[2]
distance_index_file <- args[3]
cent_start <- as.numeric(args[4])

input_file <- "chr1_customized_exclude_cent.ld"
distance_index_file <- "chr1_customized_exclude_cent.csv"
cent_start <- 121619172

ld_data <- read.table(input_file,header=TRUE)

unique_bp <- unique(c(ld_data$BP_A, ld_data$BP_B))
bp_index <- data.frame(BP = unique_bp, index = seq_along(unique_bp))
write.csv(bp_index, distance_index_file, row.names = FALSE)

cen_pos <- bp_index %>% filter(BP >= cent_start) %>% pull(index) %>% min()
print(cen_pos)

ld_data <- merge(ld_data, bp_index, by.x = "BP_A", by.y = "BP", all.x = TRUE)
colnames(ld_data)[ncol(ld_data)] <- "index_A"

ld_data <- merge(ld_data, bp_index, by.x = "BP_B", by.y = "BP", all.x = TRUE)
colnames(ld_data)[ncol(ld_data)] <- "index_B"

print(ld_data[1:10, ])


# ggplot(ld_data,aes(x=index_A,y=index_B,color=R2))+
# 	geom_point(size=0.2)+
# 	#geom_tile()+
# 	scale_color_gradient(low="#3f66a1",high="red") + 
#         scale_alpha(range = c(0.4, 1)) + 
# 	theme_minimal() +
# 	labs(fill="R2")

# p1 <- ggplot(ld_data, aes(x = index_A, y = index_B, color = R2)) +
#   geom_point(size = 0.2) +
#   #geom_tile()+
#   scale_color_gradient(low = "#3f66a1", high = "red") + 
#   scale_alpha(range = c(0.4, 1)) + 
#   theme_classic() +
#   labs(fill = "R2") +
#   geom_segment(aes(x = cen_pos, xend = cen_pos + 100, y = 0, yend = 0), color = "black", linewidth=2)+
#   scale_x_continuous(breaks = seq(0, max(ld_data$index_A), by = 200)) +
#   scale_y_continuous(breaks = seq(0, max(ld_data$index_B), by = 200)) +
#   theme(legend.position = "none")

#p1 <- ggplot(ld_data, aes(x = index_A, y = index_B, color = R2)) +
#  geom_point(size = 0.2) +
#  #geom_tile()+
#  scale_color_gradient2(low = "#4575b4", mid = "#fee090", high = "#d73027", midpoint = 0.5, limits = c(0, 1), na.value = "white") +
#  theme_classic() +
#  labs(fill = "R2") +
#  geom_segment(aes(x = cen_pos, xend = cen_pos + 100, y = 0, yend = 0), color = "black", linewidth=2) +
#  scale_x_continuous(breaks = seq(0, max(ld_data$index_A), by = 200)) +
#  scale_y_continuous(breaks = seq(0, max(ld_data$index_B), by = 200)) +
#  theme(legend.position = "none")

p1 <- ggplot(ld_data, aes(x = index_A, y = index_B, color = R2)) +
  geom_point(size = 0.2) +
  scale_color_gradientn(
    colors = rev(RColorBrewer::brewer.pal(11, "Spectral")), # 使用 Spectral 调色板
    values = c(0, 0.5, 1), # 指定颜色位置：0 对应低值，0.5 对应中间值，1 对应高值
    limits = c(0, 1), # 数据范围
    na.value = "white" # 缺失值颜色
  ) +
  theme_classic() +
  labs(fill = "R2") +
  geom_segment(aes(x = cen_pos, xend = cen_pos + 100, y = 0, yend = 0), color = "black", linewidth = 2) +
  scale_x_continuous(breaks = seq(0, max(ld_data$index_A), by = 200)) +
  scale_y_continuous(breaks = seq(0, max(ld_data$index_B), by = 200)) +
  theme(legend.position = "none")

# ggsave(output_file, dpi=300)
# print(p1)
ggsave(output_file, plot = p1, device = "png", dpi=600, width = 5, height = 5)
