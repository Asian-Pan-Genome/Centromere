library("dplyr")
library("tidyr")
library("tidyverse")

############HiCAT#################
hicat_df <- read.csv("/share/home/zhanglab/user/sunyanqing/vscode_scripts/HORmining/20250519_hicat_hor_final_summary.xls", sep="\t", header=TRUE)
head(hicat_df)

hicat_df <- hicat_df %>%
    filter(HORarray_reported_by_CHM13 != "Novel(single_mn)" & HORarray_reported_by_CHM13 != "Novel(dimeric_mns)") %>%
    separate_rows(chrom_manual_check, sep = "_")

##check##
print(hicat_df %>% filter(reorder_hor == "21663_29601_21663_29759"))

##stat_hor_num_each_chrom##
hicat_hor_num_df <- hicat_df %>%
    group_by(chrom_manual_check, horclass) %>%
    summarise(hor_num = n(), .groups = "drop") %>%
    arrange(chrom_manual_check, desc(hor_num))
head(hicat_hor_num_df)
write.table(hicat_hor_num_df, file = "/share/home/zhanglab/user/sunyanqing/vscode_scripts/HORmining/circos/hicat_hor_num_each_chrom.xls", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


############HORmon#################
hormon_df <- read.csv("/share/home/zhanglab/user/sunyanqing/vscode_scripts/HORmining/20250519_graph_hor_final_summary.xls", sep="\t", header=TRUE)
head(hormon_df)

hormon_df <- hormon_df %>%
    filter(HORarray_reported_by_CHM13 != "Novel(single_mn)" & HORarray_reported_by_CHM13 != "Novel(dimeric_mns)") %>%
    separate_rows(chrom_manual_check, sep = "_")

print(hormon_df %>% filter(hor == "6948_9927_14362_15040_15486_15282_15087_14780_15271_14732_15418_8309_14932_15433_14824_15158_15505_14867_15618_8326"))
##stat_hor_num_each_chrom##
hormon_hor_num_df <- hormon_df %>%
    group_by(chrom_manual_check, hicat_horclass) %>%
    summarise(hor_num = n(), .groups = "drop") %>%
    arrange(chrom_manual_check, desc(hor_num))
head(hormon_hor_num_df)
write.table(hormon_hor_num_df, file = "/share/home/zhanglab/user/sunyanqing/vscode_scripts/HORmining/circos/hormon_hor_num_each_chrom.xls", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
