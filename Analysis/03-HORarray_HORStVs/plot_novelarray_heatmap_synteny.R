setwd("/share/home/zhanglab/user/sunyanqing/human/alpha/part2/ancient_HOR_array/212/")

library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("patchwork")
library(ComplexHeatmap)
library(circlize)


data <- read.csv("/share/home/zhanglab/user/sunyanqing/human/alpha/part2/ancient_HOR_array/212/C094-CHA-C14_Pat_vs_panpan_pairwise.xls", header=FALSE, sep="\t")
data <- read.csv("/share/home/zhanglab/user/sunyanqing/human/alpha/part2/ancient_HOR_array/212/C094-CHA-C14_Pat_vs_gor_pairwise.xls", header=FALSE, sep="\t")
data <- read.csv("/share/home/zhanglab/user/sunyanqing/human/alpha/part2/ancient_HOR_array/212/CHM13_vs_gor_pairwise.xls", header=FALSE, sep="\t")
data <- read.csv("/share/home/zhanglab/user/sunyanqing/human/alpha/part2/ancient_HOR_array/212/C003-CHA-E03_Pat_vs_gor_pairwise.xls", header=FALSE, sep="\t")
data <- read.csv("/share/home/zhanglab/user/sunyanqing/human/alpha/part2/ancient_HOR_array/249/C119-CKZ01#Mat_vs_self_pairwise.xls", header=FALSE, sep="\t")
data <- read.csv("/share/home/zhanglab/user/sunyanqing/human/alpha/part2/ancient_HOR_array/249/C119-CKZ01#Mat_vs_pantro_pairwise.xls", header=FALSE, sep="\t")

colnames(data) <- c("human", "ape", "sim")
head(data, 10)
print(min(data$sim))
print(max(data$sim))

# Define the array of values
values_to_remove <- c(196, 215, 234, 253, 272, 291, 310, 329, 348, 367)

# Filter the data to remove rows where 'human' or 'ape' contains any value in the array
data <- data %>%
  filter(!(human %in% values_to_remove | ape %in% values_to_remove))

# Check the resulting data
head(data, 10)




wide_data <- data %>%
  pivot_wider(names_from = human, values_from = sim)
head(wide_data, 10)
dim(wide_data)

wide_data <- data %>%
  pivot_wider(names_from = ape, values_from = sim)
head(wide_data, 10)
dim(wide_data)

heatmap_matrix <- as.matrix(wide_data[,-1])
##gorilla
heatmap_matrix <- as.matrix(wide_data[1:1000, -1])
heatmap_matrix <- as.matrix(wide_data[10:nrow(wide_data), 1:1000])
##panpan
heatmap_matrix <- as.matrix(wide_data[1:710, -1])

col_fun <- colorRamp2(c(min(heatmap_matrix, na.rm = TRUE), 0.9, max(heatmap_matrix, na.rm = TRUE)), 
                      c("#4575b4", "white", "#a50026"))

col_fun <- colorRamp2(c(0.9, 0.93, max(heatmap_matrix, na.rm = TRUE)), 
            c("white", "#fee0d2", "#a50026"))

col_fun <- colorRamp2(c(0.9, 0.93, max(heatmap_matrix, na.rm = TRUE)), 
            c("#023858", "#e0f3db", "white"))


pdf("/share/home/zhanglab/user/sunyanqing/human/alpha/part2/ancient_HOR_array/249/C119-CKZ01#Mat_vs_self_pairwise.dotplot.pdf", width = 6, height = 6) # Specify the file name and dimensions
Heatmap(heatmap_matrix, 
        col = col_fun, 
        show_column_names = FALSE, 
        show_row_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        )

dev.off()
