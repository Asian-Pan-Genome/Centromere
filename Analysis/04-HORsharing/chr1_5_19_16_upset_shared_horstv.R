setwd("/share/home/zhanglab/user/sunyanqing/vscode_scripts/HORsharing/chr1-5-19")
getwd()
rm(list = ls())

library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library(UpSetR)


##load horstv ##
horstv <- read.csv("/share/home/zhanglab/user/sunyanqing/vscode_scripts/HORmining/20250519_hicat_hor_final_summary.xls", header = T, sep = "\t")
head(horstv)


shared_stv <- horstv %>%
    filter(horclass == 25,
        grepl('_', chrom_manual_check) ) %>%
    mutate(chrom_manual_check_reordered = 
        sapply(strsplit(chrom_manual_check, "_"), 
                function(x) paste(sort(x), collapse = "_")))
head(shared_stv)

chr1_shared_stv <- shared_stv %>%
    filter(grepl('chr1_', chrom_manual_check_reordered)) %>%
    select(reorder_hor) %>%
    pull(reorder_hor)
print(length(chr1_shared_stv))
print(chr1_shared_stv)

chr5_shared_stv <- shared_stv %>%
    filter(grepl('chr5', chrom_manual_check)) %>%
    select(reorder_hor) %>%
    pull(reorder_hor)
print(length(chr5_shared_stv))

chr19_shared_stv <- shared_stv %>%
    filter(grepl('chr19', chrom_manual_check)) %>%
    select(reorder_hor) %>%
    pull(reorder_hor)
print(length(chr19_shared_stv))

chr16_shared_stv <- shared_stv %>%
    filter(grepl('chr16', chrom_manual_check)) %>%
    select(reorder_hor) %>%
    pull(reorder_hor)
print(length(chr16_shared_stv))
print(chr16_shared_stv)

listInput <- list(chr1 = chr1_shared_stv, chr5 = chr5_shared_stv, chr19 = chr19_shared_stv, chr16 = chr16_shared_stv)
print(listInput)

pdf("chr1_5_16_19_shared_stv_hicat_upset_ordered.pdf", width = 8, height = 6)
upset(fromList(listInput), 
    point.size = 4, 
    line.size = 1, 
    mainbar.y.label = "The number of shared HOR StVs", 
    sets.x.label = "Total number of shared HOR StVs",
    order.by = "freq"
    )
dev.off()    


# Convert to sets for easier operations
chr1_set <- unique(chr1_shared_stv)
chr5_set <- unique(chr5_shared_stv)
chr19_set <- unique(chr19_shared_stv)
chr16_set <- unique(chr16_shared_stv)

# Calculate the complex intersections
# Only shared by chr1 and chr19 (excluding elements in larger sets)
only_chr1_chr19 <- setdiff(
  intersect(chr1_set, chr19_set),
  union(chr5_set, chr16_set)
)
print(only_chr1_chr19)

# Only shared by chr5 and chr19 (excluding elements in larger sets)
only_chr5_chr19 <- setdiff(
  intersect(chr5_set, chr19_set),
  union(chr1_set, chr16_set)
)
print(only_chr5_chr19)


# Only shared by chr1 and chr5 (excluding elements in larger sets)
only_chr1_chr5 <- setdiff(
  intersect(chr1_set, chr5_set),
  union(chr19_set, chr16_set)
)

# Only shared by chr1, 5, 19 (excluding chr16)
only_chr1_5_19 <- setdiff(
  Reduce(intersect, list(chr1_set, chr5_set, chr19_set)),
  chr16_set
)

# Only shared by chr1, 19, 16 (excluding chr5)
only_chr1_19_16 <- setdiff(
  Reduce(intersect, list(chr1_set, chr19_set, chr16_set)),
  chr5_set
)

# Only shared by chr1, 5, 19, 16 (all four)
only_chr1_5_19_16 <- Reduce(intersect, list(chr1_set, chr5_set, chr19_set, chr16_set))

# Print results
print("Only shared by chr1 and chr19:")
print(length(only_chr1_chr19))
print(only_chr1_chr19)

print("Only shared by chr5 and chr19:")
print(length(only_chr5_chr19))
print(only_chr5_chr19)

print("Only shared by chr1 and chr5:")
print(length(only_chr1_chr5))
print(only_chr1_chr5)

print("Only shared by chr1, 5, 19:")
print(length(only_chr1_5_19))
print(only_chr1_5_19)

print("Only shared by chr1, 19, 16:")
print(length(only_chr1_19_16))
print(only_chr1_19_16)

print("Only shared by chr1, 5, 19, 16:")
print(length(only_chr1_5_19_16))
print(only_chr1_5_19_16)