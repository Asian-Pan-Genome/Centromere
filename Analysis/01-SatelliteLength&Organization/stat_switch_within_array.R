setwd("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/satellite_strand_switch")
rm(list = ls())

library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")


##load strand switch data of satellites##
apg_df  <- read.csv("apg_satellite_inversion_breakpoint_count.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
other_df <- read.csv("HPRC_HGSVC_satellite_inversion_breakpoint_count.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

apg_df <- apg_df %>% 
  mutate(project = "APG")

other_df <- other_df[, c("chrom", "chr", "breakpoint_preference", "count", "project")]


head(apg_df)
head(other_df)

##merge two dataframes##
df <- rbind(apg_df, other_df)
head(df)
print(length(unique(df$chrom)))

##add cent info ##
df_cen <- read.csv("../centhap/cent_chrom.txt", header=TRUE, sep="\t")
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

complete_df_cen <- df_cen %>%
 filter(filterflag == 0) %>% 
 filter( ! sample_hap_chrom %in% filter_chrs, project != "Ref")
head(complete_df_cen)
print(length(unique(complete_df_cen$sample_hap_chrom)))

##filter incomplete cent in df ##
complete_df <- df %>%
  filter(chrom %in% complete_df_cen$sample_hap_chrom)
print(length(unique(complete_df$chrom)))


##check##
unique_values <- setdiff(complete_df_cen$sample_hap_chrom, df$chrom)
print(unique_values)

##get general pattern of satellite strand switch##
gp_df <- df %>%
    group_by(chr, breakpoint_preference) %>%
    summarise(
        shared_haploids = n_distinct(chrom),
        average_switch_count = mean(count, na.rm = TRUE),
        sd_switch_count = sd(count, na.rm = TRUE),
    ) %>%
    arrange(desc(average_switch_count), desc(shared_haploids))

print(gp_df)
write.table(gp_df, file = "all_satellite_strand_switch_general_pattern.xls", row.names = FALSE, sep="\t", quote=FALSE)

##get apg general pattern##
apg_gp_df <- df %>%
    filter(project == "APG") %>%
    group_by(chr, breakpoint_preference) %>%
    summarise(
        shared_haploids = n_distinct(chrom),
        average_switch_count = mean(count, na.rm = TRUE),
        sd_switch_count = sd(count, na.rm = TRUE),
    ) %>%
    arrange(desc(average_switch_count), desc(shared_haploids))

print(apg_gp_df)
write.table(apg_gp_df, file = "apg_satellite_strand_switch_general_pattern.xls", row.names = FALSE, sep="\t", quote=FALSE)

##get bsat pattern##

##get bsat length
satlen_df <- read.csv("../satellite_length/all.sat.length.xls", header=TRUE, sep="\t", stringsAsFactors = FALSE)
head(satlen_df)
bsatlen_df <- satlen_df %>%
  filter(satellite == "BSat") %>%
  select(sample_hap_chrom, length) %>%
  rename(bsat_length = length)
head(bsatlen_df)


##global pattern##
head(df)
head(complete_df_cen)

bsat_merged_df <- df %>%
  left_join(complete_df_cen %>% select(sample_hap, sample, sample_hap_chrom, hap, filterflag),
              by = c("chrom" = "sample_hap_chrom")) %>%
  filter(filterflag == 0,
         breakpoint_preference == "within BSat array") %>%
  left_join(bsatlen_df, by = c("chrom" = "sample_hap_chrom")) %>%
  mutate(switch_freq = count / bsat_length * 1000000)
  
head(bsat_merged_df)
write.table(bsat_merged_df, file = "global_bsat_strand_switch.xls", row.names = FALSE, sep="\t", quote=FALSE)

##plot for each chrom##
chrom_colors <- c(
  "chr1" = "#9D5427", "chr2" = "#D67C1B", "chr3" = "#C6A15B", "chr4" = "#F2E86D",
  "chr5" = "#030303", "chr6" = "#FF6533", "chr7" = "#FF4B3E", "chr8" = "#CB3C23",
  "chr9" = "#972D07", "chr10" = "#a4c558", "chr11" = "#59b669", "chr12" = "#3c8340",
  "chr13" = "#6DAB30", "chr14" = "#285943", "chr15" = "#7ec0b4", "chr16" = "#85A7CC",
  "chr17" = "#3f66a1", "chr18" = "#3897c5", "chr19" = "#afa7d8", "chr20" = "#a874b5",
  "chr21" = "#893f8b", "chr22" = "#b13f73", "chrX" = "#574D68", "chrY" = "#C52B94"
)

chrom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 'chrX', 'chrY')
bsat_merged_df$chr <- factor(bsat_merged_df$chr, levels = chrom_order)

target_chrs <- c("chr1", "chr13", "chr14", "chr15", "chr21", "chr22")
target_bsat_df <- bsat_merged_df %>%
  filter(chr %in% target_chrs)

cor_test <- cor.test(target_bsat_df$count, target_bsat_df$bsat_length / 1000000, method = "pearson")
pearson_r <- round(cor_test$estimate, 3)  # 提取相关系数 (R)
pearson_p <- signif(cor_test$p.value, 10) # 提取 p-value
print(pearson_r)
print(pearson_p)


pbsat <- ggplot(bsat_merged_df %>% filter(chr %in% target_chrs), aes(x = bsat_length / 1000000, y = count, color = chr)) +
  geom_point(size = 2, alpha = 0.5) + 
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2, color = "black") +
  labs(
    x = "BSat Length (Mb)",
    y = "Switch Count",
    color = "Chromosome"
  ) +
  scale_color_manual(values = chrom_colors) +
  theme_classic() +
  annotate(
    "text",
    x = 0.1,  # Adjust x position
    y = 1500,                  # Adjust y position
    label = paste0("r = ", pearson_r, "\nR² = ", round(pearson_r^2, 2), "\nP = ", pearson_p),
    size = 4,
    hjust = 0
  )

print(pbsat)
ggsave("bsat_strand_switch_chrom.pdf", plot = pbsat, width = 6, height = 5)

##get hsat3 patter##

hsat3len_df <- read.csv("/share/home/zhanglab/user/fukezhi/yanqing_tasks/all.chr9.hsat3.breakpoint.num.stats/all.info.chr9.hsat3.breakpoint.stats.txt",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, fileEncoding = "UTF-16LE")
head(hsat3len_df)
write.table(hsat3len_df, file = "global_hsat3_strand_switch.xls", row.names = FALSE, sep="\t", quote=FALSE, fileEncoding = "UTF-8")


cor_test <- cor.test(hsat3len_df$breakpoint_counts, hsat3len_df$HSat3_length / 1000000, method = "pearson")
pearson_r <- round(cor_test$estimate, 3)  # 提取相关系数 (R)
pearson_p <- signif(cor_test$p.value, 3) # 提取 p-value
print(pearson_r)
print(pearson_p)


phsat3 <- ggplot(hsat3len_df, aes(x = HSat3_length / 1000000, y = breakpoint_counts, color = project)) +
  geom_point(size = 2, alpha = 0.5) + 
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2, color = "black") +
  labs(
    x = "HSat3 Length on chr9 (Mb)",
    y = "Switch Count",
    color = "project"
  ) +
  theme_classic() +
  annotate(
    "text",
    x = 0.1,  # Adjust x position
    y = 400,                  # Adjust y position
    label = paste0("r = ", pearson_r, "\nR² = ", round(pearson_r^2, 2), "\nP = ", pearson_p),
    size = 4,
    hjust = 0
  )

print(phsat3)
ggsave("hsat3_strand_switch_chrom.pdf", plot = phsat3, width = 6, height = 5)


#### stat asat ####
asat_horclass_switch_df <- read.csv("all.HORclass.inversion_breakpoint.bed", header=TRUE, sep="\t", stringsAsFactors = FALSE)
colnames(asat_horclass_switch_df) <- c("sample_hap_chrom", "start", "end", "interval", "switch_type", "horclass_anno")
asat_horclass_switch_df <- asat_horclass_switch_df[, -ncol(asat_horclass_switch_df)]
head(asat_horclass_switch_df)
print(length(unique(asat_horclass_switch_df$sample_hap_chrom)))


filtered_asat_horclass_switch_df <- asat_horclass_switch_df %>%
  left_join(complete_df_cen %>% select(sample_hap, sample, chrom, sample_hap_chrom, hap, project, filterflag),
            by = "sample_hap_chrom") %>%
  filter(filterflag == 0,
         switch_type != "+/+",
         switch_type != "-/-") %>%
  mutate(
    switch_anno = case_when(
      horclass_anno == "999/999" ~ "monomeric_switch",
      grepl("999", horclass_anno) ~ "near_monomeric",
      map_chr(str_split(horclass_anno, "/"), ~ .x[1]) != map_chr(str_split(horclass_anno, "/"), ~ .x[2]) ~ "array transition switch",
      TRUE ~ horclass_anno
    )
  )
print(length(unique(filtered_asat_horclass_switch_df$sample_hap_chrom)))
head(filtered_asat_horclass_switch_df, n =22)

##get general pattern of asat strand switch##
asat_switch_stat <- filtered_asat_horclass_switch_df %>%
  group_by(sample_hap_chrom, sample, sample_hap, hap, chrom, switch_anno, project) %>%
  summarise(
    switch_count = n(),
    .groups = 'drop'
  )
head(asat_switch_stat)

asat_switch_sum <- asat_switch_stat %>%
  group_by(sample_hap) %>%
  summarise(
    total_switch_count = sum(switch_count, na.rm = TRUE),
    monomeric_switch_count = sum(switch_count[switch_anno %in% c("monomeric_switch", "near_monomeric")], na.rm = TRUE),
    transition_switch_count = sum(switch_count[switch_anno %in% c("array transition switch")], na.rm = TRUE)) %>%
   mutate(monomeric_switch_percentage = monomeric_switch_count / total_switch_count * 100,
          transition_switch_percentage = transition_switch_count / total_switch_count * 100)
head(asat_switch_sum)

print(length(unique(asat_switch_sum$sample_hap)))
print(mean(asat_switch_sum$monomeric_switch_percentage))
print(mean(asat_switch_sum$transition_switch_percentage))


asat_switch_summary <- asat_switch_stat %>%
  group_by(chrom, switch_anno) %>%
  summarise(
    shared_haploids = n_distinct(sample_hap_chrom),
    average_switch_count = mean(switch_count, na.rm = TRUE),
    sd_switch_count = sd(switch_count, na.rm = TRUE)
  ) %>%
  arrange(desc(average_switch_count), desc(shared_haploids))
head(asat_switch_summary)
write.table(asat_switch_summary, file = "global_asat_HORclass_switch_general_pattern.xls", row.names = FALSE, sep="\t", quote=FALSE)

asat_switch_summary <- asat_switch_stat %>%
  filter(project == "APG") %>%
  group_by(chrom, switch_anno) %>%
  summarise(
    
    shared_haploids = n_distinct(sample_hap_chrom),
    average_switch_count = mean(switch_count, na.rm = TRUE),
    sd_switch_count = sd(switch_count, na.rm = TRUE)
  ) %>%
  arrange(desc(average_switch_count), desc(shared_haploids))
head(asat_switch_summary)
write.table(asat_switch_summary, file = "global_asat_HORclass_switch_general_pattern.xls", row.names = FALSE, sep="\t", quote=FALSE)


###specific inversion of chr1 active array##
popfile <- "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/populaion.xls"
popdf <- read.csv(popfile, header=TRUE, sep="\t")
popdf <- popdf %>% 
  select(c("sample", "superpopulation")) %>%
  mutate(superpopulation = str_replace(superpopulation, "EAS-APG", "EAS"))
head(popdf,10)


chr1_active_array_df <- read.csv("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/satellite_strand_switch/chr1_25/chr1_25_inversion_region_length.xls", sep="\t", header=FALSE, stringsAsFactors = FALSE)
colnames(chr1_active_array_df) <- c("sample_hap_chrom", "start", "end", "strand", "horclass", "rid", "length")
head(chr1_active_array_df)

complete_chr1_df <- chr1_active_array_df %>%
  left_join(df_cen %>% select(sample_hap, sample, chrom, sample_hap_chrom, hap, project, filterflag),
            by = "sample_hap_chrom") %>%
  filter( ! sample_hap_chrom %in% filter_chrs) %>%
  filter(filterflag == 0) 
head(complete_chr1_df)

print(length(unique(complete_chr1_df$sample_hap_chrom)))
print(length(unique(chr1_active_array_df$sample_hap_chrom)))

diff_values <- setdiff(unique(chr1_active_array_df$sample_hap_chrom), unique(complete_chr1_df$sample_hap_chrom))
print(diff_values)

print(mean(complete_chr1_df %>% filter(strand == "-") %>% pull(length)))
print(sd(complete_chr1_df %>% filter(strand == "-") %>% pull(length)))

###plot violin###
complete_chr1_df$rid <- as.factor(complete_chr1_df$rid)  
p1 <- ggplot(complete_chr1_df, aes(x = rid, y = length/1000000, fill = rid)) +
  geom_violin() +
  labs(
       x = "RID",
       y = "Length") +
  theme_classic()
print(p1)
ggsave(p1, filename = "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/satellite_strand_switch/chr1_25/rid_length_violin.pdf", width = 8, height = 8)

chr1_switch_df <- filtered_asat_horclass_switch_df %>%
  filter(chrom == "chr1", switch_anno != "monomeric_switch", switch_anno != "near_monomeric")
head(chr1_switch_df)
print(unique(chr1_switch_df$horclass_anno))

chr1_active_switch_num <- chr1_switch_df %>%
  left_join(popdf, by="sample") %>%
  filter(horclass_anno == "25/25") %>%
  group_by(sample_hap, sample, superpopulation) %>%
  summarise(
    total_switch_count = n()
  ) 
diff_values <- setdiff(unique(chr1_active_switch_num$sample_hap), unique(complete_chr1_df$sample_hap))
print(diff_values)
head(chr1_active_switch_num)
dim(chr1_active_switch_num)

chr1_active_switch_type <- chr1_active_switch_num %>%
  group_by(total_switch_count) %>%
  summarise(
    sample_hap_count = n(),
    .groups = 'drop'
  )
head(chr1_active_switch_type, n =30)

print(chr1_active_switch_num %>% filter(total_switch_count >= 4))

head(complete_df_cen)
chr1_complete_superpopulation_num <- complete_df_cen %>%
  filter(chrom == "chr1") %>%
  left_join(popdf, by="sample") %>%
  group_by(superpopulation) %>%
  summarise(
    sample_hap_count = n_distinct(sample_hap_chrom)
  )
print(chr1_complete_superpopulation_num)

