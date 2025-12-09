setwd("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap")
rm(list = ls())

library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("patchwork")
library(ggnewscale)

chr <- "chr19"

##load phynode-sample_hap file##
nodefile <- "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/Phy/phynode_sample_hap.xls"
node_df <- read.csv(nodefile, header=FALSE, sep="\t")
colnames(node_df) <-  c("node", "sample_hap")
head(node_df,10)

##load left phy order file##
left_phyfile <- paste0("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/Phy/", chr, "/", chr,"_left_ordered_phy.xls")
left_orderdf <- read.csv(left_phyfile, header=FALSE, sep="\t")
colnames(left_orderdf) <- c("left_node", "left_cluster")
head(left_orderdf,10)

##load right phy order file##
right_phyfile <- paste0("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/Phy/",chr, "/", chr,"_right_ordered_phy.xls")
right_orderdf <- read.csv(right_phyfile, header=FALSE, sep="\t")
colnames(right_orderdf) <- c("right_node", "right_cluster")
head(right_orderdf,10)

##load cent info file##
filtered_cent <- c(
  "C019-CHA-E19#Pat#chr9",
  "C045-CHA-N05#Mat#chr17",
  "C076-CHA-NE16#Pat#chr14",
  "C020-CHA-E20#Mat#chr6",
  "HG02666_hap2_chr15",
  "HG01114_hap1_chr16",
  "HG02769_hap2_chr20",
  "HG03452_hap1_chr4",
  "NA19036_hap2_chr14",
  "NA19434_hap2_chr20",
  "NA20847_hap2_chr17"
)
df_cen <- read.csv("cent_chrom.txt", header=TRUE, sep="\t")
df_cen <- df_cen %>% filter(chrom == chr) %>% 
  filter(!sample_hap_chrom %in% filtered_cent) %>%
  mutate(sample_hap = str_replace(sample_hap, "CHM13_-", "CHM13")) %>%
  mutate(sample_hap_chrom = str_replace(sample_hap_chrom, "CHM13#", ""))
colnames(df_cen) <-  c("sample_hap_chrom", "sample_hap", "sample", "hap", "chrom", "cen_start", "cen_end", "cen_len", "project", "gapless", "CL", "filterflag")
head(df_cen,10)

##load pop file##
popfile <- "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/populaion.xls"
popdf <- read.csv(popfile, header=TRUE, sep="\t")
popdf <- popdf %>% 
  select(c("sample", "superpopulation")) %>%
  mutate(superpopulation = str_replace(superpopulation, "EAS-APG", "EAS"))
head(popdf,10)

##load cenhap file##
cenhapfile <- paste0("/share/home/zhanglab/user/sunyanqing/vscode_scripts/cenhap/", chr, "_cenhap.xls")
cenhap_df <- read.csv(cenhapfile, header=TRUE, sep="\t")
head(cenhap_df)
print(unique(cenhap_df$Final_subcluster))


##load te file##
te_color_df <- read.csv("te_color.xls", header=T, sep="\t")
head(te_color_df)

tefile <- paste0("TE/", chr, "_te.bed")
te_df <- read.csv(tefile, header=FALSE, sep="\t")
colnames(te_df) <- c("project", "superpopulation", "population", "sample", "hap", "chrom", "te_start", "te_end", "length", "strand", "motif", "class", "family", "type")
head(te_df,10)

print(unique(te_df$sample))
asat_te_df <- te_df %>% 
  filter(type == "ASat") %>%
  mutate(sample_hap = paste0(sample, "_", hap)) %>%
  mutate(sample_hap = str_replace(sample_hap, "CHM13_-", "CHM13")) %>%
  select(c("sample_hap", "te_start", "te_end", "motif")) %>%
  mutate(te_color = ifelse(motif %in% te_color_df$motif, 
                            te_color_df$color[match(motif, te_color_df$motif)], 
                            "grey"))

head(asat_te_df,10)
print(unique(asat_te_df$sample_hap))


stat_te_df <- asat_te_df %>%
  group_by(motif, te_color) %>%
  summarise(count = n()) %>%
  arrange(desc(count))
head(stat_te_df, n=20)

##merge left and right phy order file##
left_df <- left_orderdf %>%
  inner_join(right_orderdf, by = c("left_node" = "right_node")) %>%
  inner_join(node_df, by = c("left_node" = "node")) %>%
  inner_join(df_cen, by = "sample_hap") %>%
  filter(filterflag == 0) %>%
  inner_join(popdf, by = "sample")
head(left_df,10)
dim(left_df)
write.table(left_df, file = paste0("Phy/", chr, "/", chr, "_left_ordered_out.xls"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

left_trimed <- left_orderdf %>% 
  filter(! left_node %in% left_df$left_node) %>%
  select(left_node) %>%
  mutate(left_node = str_replace(left_node, " ", "_"))
head(left_trimed,10)
dim(left_trimed)
write.table(left_trimed, file = paste0("Phy/", chr, "/", chr, "_left_trimed.xls"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

left_superpop_df <- left_df %>%
  select(left_node, superpopulation) %>%
  mutate(left_node = str_replace(left_node, " ", "_")) %>%
  group_by(left_node, superpopulation) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = superpopulation, values_from = count, values_fill = 0)
head(left_superpop_df,10)
write.table(left_superpop_df, file = paste0("Phy/", chr, "/", chr, "_superpop.xls"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


node_cenhap_df <- left_df %>%
  select(left_node, sample_hap) %>%
  mutate(left_node = str_replace(left_node, " ", "_")) %>%
  inner_join(cenhap_df, by = "sample_hap") %>%
  group_by(left_node, Final_subcluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Final_subcluster, values_from = count, values_fill = 0)
head(node_cenhap_df,10)
dim(node_cenhap_df)
write.table(node_cenhap_df, file = paste0("Phy/", chr, "/", chr, "_cenhap_anno.xls"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

right_df <- right_orderdf %>%
  inner_join(left_orderdf, by = c("right_node" = "left_node")) %>%
  inner_join(node_df, by = c("right_node" = "node")) %>%
  inner_join(df_cen, by = "sample_hap") %>%
  filter(filterflag == 0) %>%
  inner_join(popdf, by = "sample")
head(right_df,10)
dim(right_df)
write.table(right_df, file = paste0("Phy/", chr, "/", chr, "_right_ordered_out.xls"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

right_trimed <- right_orderdf %>% 
  filter(! right_node %in% right_df$right_node) %>%
  select(right_node) %>%
  mutate(right_node = str_replace(right_node, " ", "_"))
head(right_trimed,10)
dim(right_trimed)
write.table(right_trimed, file = paste0("Phy/", chr, "/", chr, "_right_trimed.xls"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


####### update left and right phy order #########
left_updated_order <- read.csv(paste0("Phy/", chr, "/", chr, "_left_phy_updated_ordered.xls"), header=FALSE, sep="\t")
colnames(left_updated_order) <- c("left_node")
left_updated_order <- left_updated_order %>% mutate(left_y = n() - row_number())
head(left_updated_order)

right_updated_order <- read.csv(paste0("Phy/", chr, "/", chr, "_right_phy_updated_ordered.xls"), header=FALSE, sep="\t")
colnames(right_updated_order) <- c("right_node")
right_updated_order <- right_updated_order %>% mutate(right_y = n() - row_number())

####### HiCAT ############
infile<-paste0(chr, ".HiCAT.horstv.bed")
data <- read.csv(infile, header=FALSE, sep="\t")
colnames(data) <- c("sample_hap_chrom", "horstart", "horend", "hor_mn", "horclass", "horclass_color", "horstv_index", "horstv_color")
print(data[20000, ])

color_map <- c(
  "#dfc27d" = "#C5C8D9",
  "#B66E64" = "#1a9850",
  "#F25CA2" = "#8C3074",
  "#4a37b6" = "#F25CA2",
  "#f46d43" = "#058187",
  "#446b4c" = "#F2B33D",
  "#46c4c4" = "#0C4DA7",
  "#877DD4" = "#01325E",
  "#A50161" = "#D9D2BF"
)

recolor_file <- "/share/home/zhanglab/user/sunyanqing/vscode_scripts/plot_scripts/all.horstv.hicat.update.color.txt"
recolordf  <- read.csv(recolor_file, header=TRUE, sep="\t")
colnames(recolordf) <- c("horstv_index", "horstv_newcolor")
recolordf <- recolordf %>%
mutate(horstv_newcolor = recode(horstv_newcolor, !!!color_map))
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
  left_join(recolordf, by = "horstv_index") %>%
  mutate(horstv_newcolor = ifelse(horstv_color == '-', "#525252",
                                  ifelse(!is.na(horstv_newcolor), horstv_newcolor, horstv_color))) %>%
  select(sample_hap_chrom, sample_hap, horstart, horend, horstv_index, horstv_newcolor)
print(head(stv, 10))
print(unique(stv$sample_hap))
print(length(unique(stv$sample_hap_chrom)))


stv_min_hor_start <- stv %>%
  group_by(sample_hap_chrom) %>%
  summarise(min_hor_start = min(horstart, na.rm = TRUE)) %>%
  ungroup()
head(stv_min_hor_start,10)

left_stv_df <- stv %>%
  inner_join(left_df, by = c("sample_hap_chrom")) %>%
  filter(horstart >= cen_start & horend <= cen_end) %>%
  inner_join(stv_min_hor_start, by = c("sample_hap_chrom")) %>%
  mutate(
    modstart = if_else(min_hor_start < cen_start, horstart - cen_start + 1, horstart - min_hor_start + 1),
    modend = if_else(min_hor_start < cen_start, horend - cen_start + 1, horend - min_hor_start + 1)
  ) %>%
  inner_join(left_updated_order, by = c("left_node" = "left_node"))
head(left_stv_df,10)
print(length(unique(left_stv_df$sample_hap_chrom)))

left_te_df <- asat_te_df %>% 
  inner_join(node_df, by = "sample_hap") %>%
  inner_join(df_cen, by = "sample_hap") %>%
  filter(te_start >= cen_start & te_end <= cen_end) %>%
  inner_join(stv_min_hor_start, by = c("sample_hap_chrom")) %>%
    mutate(
    modstart = if_else(min_hor_start < cen_start, te_start - cen_start + 1, te_start - min_hor_start + 1),
    modend = if_else(min_hor_start < cen_start, te_end- cen_start + 1, te_end - min_hor_start + 1)
  ) %>%
  inner_join(left_updated_order, by = c("node" = "left_node"))
head(left_te_df,10)

stv_max_hor_end <- stv %>%
  group_by(sample_hap_chrom) %>%
  filter(horstv_index != "Mon") %>%
  summarise(max_hor_end = max(horend, na.rm = TRUE)) %>%
  ungroup()
head(stv_max_hor_end,10)


right_stv_df <- stv %>%
  inner_join(right_df, by = c("sample_hap_chrom")) %>%
  filter(horstart >= cen_start & horend <= cen_end) %>%
  inner_join(stv_max_hor_end, by = c("sample_hap_chrom")) %>%
  mutate(
    modstart = if_else(max_hor_end > cen_end, horstart - cen_end + 1, horstart - max_hor_end + 1),
    modend = if_else(max_hor_end > cen_end, horend - cen_end + 1, horend - max_hor_end + 1)
  ) %>%
  inner_join(right_updated_order, by = c("right_node" = "right_node"))
head(right_stv_df,10)
print(length(unique(right_stv_df$sample_hap_chrom)))

right_te_df <- asat_te_df %>% 
  inner_join(node_df, by = "sample_hap") %>%
  inner_join(df_cen, by = "sample_hap") %>%
  filter(te_start >= cen_start & te_end <= cen_end) %>%
  inner_join(stv_max_hor_end, by = c("sample_hap_chrom")) %>%
  mutate(
    modstart = if_else(max_hor_end > cen_end, te_start - cen_end + 1, te_start - max_hor_end + 1),
    modend = if_else(max_hor_end > cen_end, te_end - cen_end + 1, te_end - max_hor_end + 1)
  ) %>%
  inner_join(right_updated_order, by = c("node" = "right_node"))
head(right_te_df,10)

#######plot left ####
max_len <- max(left_stv_df$modend)
y_labels <- left_stv_df %>%distinct(sample_hap_chrom, left_y)

p1 <- ggplot() +
  geom_rect(
    data = left_stv_df,
    aes(xmin = modstart, xmax = modend, ymin = left_y - 0.3, ymax = left_y + 0.3, fill = horstv_newcolor),
    color = NA
  ) +
  scale_fill_identity() +
  new_scale("fill") +  
   geom_point(
    data = left_te_df,
    aes(
      x = (modstart + modend) / 2, 
      y = left_y + 0.4,  
      fill = te_color,    
      color = te_color 
    ),
    shape = 25,  
    size = 0.1,
    stroke = 0   
  ) +
  scale_color_identity() + 
  scale_fill_identity() +
  scale_x_continuous(
      breaks = seq(0, max_len, by = 2000000),
      labels = function(x) x / 1e6 #Mb
    ) +
  scale_y_continuous(breaks = y_labels$left_y, labels = y_labels$sample_hap_chrom) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 1),
    axis.title.y = element_text(size = 5)
  )

print(p1)


#######plot right ####
if (chr == "chr10") {
  right_stv_df <- right_stv_df %>%
    filter(horstart <= max_hor_end & horend <= max_hor_end)
}
max_len <- min(right_stv_df$modend)
y_labels <- right_stv_df %>%distinct(sample_hap_chrom, right_y)
print(max_len)

p2 <- ggplot() +
  geom_rect(
    data = right_stv_df,
    aes(xmin = modstart, xmax = modend, ymin = right_y - 0.3, ymax = right_y + 0.3, fill = horstv_newcolor),
    color = NA
  ) +
  scale_fill_identity() +
  new_scale("fill") +  
   geom_point(
    data = right_te_df,
    aes(
      x = (modstart + modend) / 2, 
      y = right_y  + 0.4,  
      fill = te_color,    
      color = te_color 
    ),
    shape = 25,  
    size = 0.1,
    stroke = 0   
  ) +
  scale_color_identity() + 
  scale_fill_identity() +
  scale_x_continuous(
      breaks = seq(ceiling(max_len / 1e6) * 1e6, 0, by = 2e6),
      labels = function(x) x / 1e6 #Mb
    ) +
  scale_y_continuous(breaks = y_labels$right_y, labels = y_labels$sample_hap_chrom,
      position = "right") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 1),
    axis.title.y = element_text(size = 5)
  )

print(p2)

## plot linked line ##
head(right_df,10)
head(left_df,10)

linked_df <- left_df %>%
  inner_join(left_updated_order, by = "left_node") %>%
  select(sample_hap, left_cluster, left_y) %>%
  inner_join(right_df %>% inner_join(right_updated_order, by = "right_node") %>% select(sample_hap, right_cluster, right_y), by = "sample_hap") 
head(linked_df,10)

sorted_linked_df <- linked_df %>%
  group_by(left_cluster, right_cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(count)
print(sorted_linked_df, n = 100)


if( chr == "chr17"){
sorted_linked_df <- sorted_linked_df %>%
  mutate(
    left_prefix = str_split(as.character(left_cluster), "\\.", simplify = TRUE)[, 1],
    right_prefix = str_split(as.character(right_cluster), "\\.", simplify = TRUE)[, 1],
    line_color = case_when(
      left_prefix == right_prefix & left_prefix == "1" ~ "#af8dc3",
      left_prefix == right_prefix & left_prefix == "2" ~ "#7fbf7b",
      TRUE ~ "#F46D43"
    )
  )%>%
    mutate(index = row_number())
  }

if (chr == "chr4") {
  sorted_linked_df <- sorted_linked_df %>%
    mutate(
      line_color = case_when(
        left_cluster == 2 & right_cluster == 1 ~ "#ABDDA4",
        left_cluster == 1 & right_cluster == 1 ~ "#ABDDA4",
        left_cluster == 4 & right_cluster == 4 ~ "#FEE08B",
        left_cluster == 3 & right_cluster == 2 ~ "#FEE08B",
        left_cluster == 4 & right_cluster == 3 ~ "#FEE08B",
        left_cluster == 4 & right_cluster == 'M2' ~ "#FEE08B",
        left_cluster == 4 & right_cluster == 2 ~ "#FEE08B",
        left_cluster == 3 & right_cluster == 4 ~ "#FEE08B",
        left_cluster == 'M' & right_cluster == 3 ~ "#FFFFBF",
        left_cluster == 4 & right_cluster == 1 ~ "#9E0142",
        left_cluster == 1 & right_cluster == 4 ~ "#F46D43",
        left_cluster == 1 & right_cluster == 3 ~ "#F46D43",
        TRUE ~ "grey"
      )
    ) %>%
    mutate(index = row_number())
}
head(sorted_linked_df, 10)  


if (chr == "chr5") {
  sorted_linked_df <- sorted_linked_df %>%
    mutate(
      line_color = case_when(
        left_cluster == 3 & right_cluster == 3 ~ "#3288BD",
        left_cluster == 2 & right_cluster == 2 ~ "#E6F598",
        left_cluster == 1 & right_cluster == 1 ~ "#E6F598",
        left_cluster == 'M2' & right_cluster == 'M1' ~ "#E6F598",
        left_cluster == 2 & right_cluster == 'M2' ~ "#E6F598",
        left_cluster == 2 & right_cluster == 1 ~ "#E6F598",
        left_cluster == 'M2' & right_cluster == 'M2' ~ "#E6F598",
        left_cluster == 1 & right_cluster == 'M2' ~ "#E6F598",
        left_cluster == 3 & right_cluster == 'M2' ~ "#FDAE61",
        left_cluster == 2 & right_cluster == 3 ~ "#9E0142",
        left_cluster == 'M1' & right_cluster == 3 ~ "#d53e4f",
        left_cluster == 1 & right_cluster == 3 ~ "#F46D43",
        TRUE ~ "grey"
      )
    ) %>%
    mutate(index = row_number())
}
print(sorted_linked_df, n = 100)


if (chr == "chr8") {
  sorted_linked_df <- sorted_linked_df %>%
    mutate(
      line_color = case_when(
        left_cluster == 3 & right_cluster == 4 ~ "#ABDDA4",
        left_cluster == 3 & right_cluster == 2 ~ "#E6F598",
        left_cluster == 3 & right_cluster == 3 ~ "#E6F598",
        left_cluster == 2 & right_cluster == 3 ~ "#E6F598",
        left_cluster == 1 & right_cluster == 1 ~ "#E6F598",
        left_cluster == 2 & right_cluster == 1 ~ "#E6F598",
        left_cluster == 2 & right_cluster == 'M' ~ "#E6F598",
        left_cluster == 2 & right_cluster == 4 ~ "#F46D43",
        left_cluster == 'M' & right_cluster == 2 ~ "#E6F598",
        left_cluster == 3 & right_cluster == 'M' ~ "#D5683F",
        left_cluster == 1 & right_cluster == 4 ~ "#d53e4f",
        left_cluster == 1 & right_cluster == 3 ~ "#EF949F",
        TRUE ~ "grey"
      )
    ) %>%
    mutate(index = row_number())
}
print(sorted_linked_df, n = 100)

if (chr == "chr10") {
  sorted_linked_df <- sorted_linked_df %>%
    mutate(
      line_color = case_when(
        left_cluster == 5 & right_cluster == 4 ~ "#3288BD",
        left_cluster == 4 & right_cluster == 4 ~ "#3288BD",
        left_cluster == 1 & right_cluster == 2 ~ "#66C2A5",
        left_cluster == 2 & right_cluster == 3 ~ "#66C2A5",
        left_cluster == 1 & right_cluster == 1 ~ "#5E4FA2",
        left_cluster == 1 & right_cluster == 3 ~ "#66C2A5",
        left_cluster == 4 & right_cluster == 3 ~ "#F46D43",
        left_cluster == 2 & right_cluster == 4 ~ "#D53E4F",
        left_cluster == 1 & right_cluster == 4 ~ "#9E0142",
        left_cluster == 5 & right_cluster == 3 ~ "#FEE08B",
        left_cluster == 4 & right_cluster == 2 ~ "#C55C10",
        TRUE ~ "grey"
      )
    ) %>%
    mutate(index = row_number())
}
print(sorted_linked_df, n = 100)

if (chr == "chr18") {
  sorted_linked_df <- sorted_linked_df %>%
    mutate(
      line_color = case_when(
        left_cluster == 4 & right_cluster == 3 ~ "#80B1D3",
        left_cluster == 4 & right_cluster == 2 ~ "#80B1D3",
        left_cluster == 3 & right_cluster == 3 ~ "#80B1D3",
        left_cluster == 3 & right_cluster == 2 ~ "#80B1D3",
        left_cluster == 2 & right_cluster == 1 ~ "#BEBADA",
        left_cluster == 1 & right_cluster == 1 ~ "#8DD3C7",
        left_cluster == 2 & right_cluster == 3 ~ "#FDB462",
        left_cluster == 4 & right_cluster == 1 ~ "#FB8072",
        left_cluster == 1 & right_cluster == 2 ~ "#E54C5E",
        left_cluster == 3 & right_cluster == 1 ~ "#EE822F",
        left_cluster == 2 & right_cluster == 2 ~ "#BEBADA",
        TRUE ~ "grey"
      )
    ) %>%
    mutate(index = row_number())
}
print(sorted_linked_df, n = 100)

if (chr == "chr19") {
  sorted_linked_df <- sorted_linked_df %>%
    mutate(
      line_color = case_when(
        left_cluster == 4 & right_cluster == 6 ~ "#80B1D3",
        left_cluster == 4 & right_cluster == 5 ~ "#80B1D3",
        left_cluster == 3 & right_cluster == 3 ~ "#80B1D3",
        left_cluster == 2 & right_cluster == 5 ~ "#FB8072",
        left_cluster == 1 & right_cluster == 1 ~ "#8DD3C7",
        left_cluster == 4 & right_cluster == 4 ~ "#80B1D3",
        left_cluster == 4 & right_cluster == 1 ~ "#B3DE69",
        left_cluster == 4 & right_cluster == 3 ~ "#BEBADA",
        left_cluster == 4 & right_cluster == 2 ~ "#BEBADA",
        left_cluster == 3 & right_cluster == 4 ~ "#2D529F",
        left_cluster == 3 & right_cluster == 5 ~ "#4874CB",
        left_cluster == 3 & right_cluster == 6 ~ "#30C0B4",
        left_cluster == 3 & right_cluster == 2 ~ "#E54C5E",
        left_cluster == 2 & right_cluster == 2 ~ "#FFFFB3",
        left_cluster == 2 & right_cluster == 3 ~ "#C55C10",
        left_cluster == 2 & right_cluster == 1 ~ "#44546A",
        TRUE ~ "grey"
      )
    ) %>%
    mutate(index = row_number())
}
print(sorted_linked_df, n = 100)

library(RColorBrewer)
if ( chr == "chr11") {
  # Define the Set3 palette with 11 colors
  set3_colors <- brewer.pal(11, "Spectral")
  
  # Assign colors based on the count
  sorted_linked_df <- sorted_linked_df %>%
    mutate(line_color = ifelse(count > 10, set3_colors[(row_number() %% 11) + 1], "grey")) %>%
    mutate(index = row_number())
}

# if ( chr == "chr19") {
#   set3_colors <- brewer.pal(7, "Set3")
  
#   # Assign colors based on the count
#   sorted_linked_df <- sorted_linked_df %>%
#     mutate(line_color = ifelse(count > 10, set3_colors[(row_number() %% 7) + 1], "grey")) %>%
#     mutate(index = row_number())
# }
print(sorted_linked_df, n = 100)

# Arrange by index
linked_df <- linked_df %>%
  inner_join(
    sorted_linked_df %>% select(left_cluster, right_cluster, line_color, index),
    by = c("left_cluster", "right_cluster")
  ) %>% arrange(index)
head(linked_df,10)

line_plot <- ggplot(linked_df) +
  geom_segment(
    aes(
      x = 1, xend = 2, 
      y = left_y, yend = right_y, 
      color = line_color
    ),
    size = 0.2
  ) +
  scale_color_identity() + 
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
  )

# Print the plot
print(line_plot)

library(cowplot)
combined_plot <- plot_grid(p1, line_plot, p2, ncol = 3, rel_widths = c(1.5, 1, 1.5))
print(combined_plot)

outpdf <- paste0("Phy/", chr, "/", chr, "_phy_HiCAT_cenhap.pdf")
print(outpdf)
ggsave(outpdf, plot = combined_plot, device = "pdf", width = 6, height = 10)



##### HORmon #####
graph_file <-paste0(chr, ".graph.hordecomposition.final.xls")
print(graph_file)
data <- read.csv(graph_file, header=FALSE, sep="\t")
colnames(data) <- c("sample_hap_chrom", "horstart", "horend", "hor_mn", "horclass", "horclass_color",  "horstv_color")
head(data)
print(data[30000, ])

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
  )
print(head(graphstv, 10))
print(unique(graphstv$sample_hap))
print(length(unique(graphstv$sample_hap_chrom)))


stv_min_hor_start <- graphstv %>%
  group_by(sample_hap_chrom) %>%
  summarise(min_hor_start = min(horstart, na.rm = TRUE)) %>%
  ungroup()
head(stv_min_hor_start,10)


left_graph_plot_stv_df <- graphstv %>%
  inner_join(left_df, by = c("sample_hap_chrom")) %>%
  filter(horstart >= cen_start & horend <= cen_end) %>%
  inner_join(stv_min_hor_start, by = c("sample_hap_chrom")) %>%
  mutate(
    modstart = if_else(min_hor_start < cen_start, horstart - cen_start + 1, horstart - min_hor_start + 1),
    modend = if_else(min_hor_start < cen_start, horend - cen_start + 1, horend - min_hor_start + 1)
  ) %>%
  inner_join(left_updated_order, by = c("left_node"))
head(left_graph_plot_stv_df,10)
print(length(unique(left_graph_plot_stv_df$sample_hap_chrom)))


left_graph_te_df <- asat_te_df %>% 
  inner_join(node_df, by = "sample_hap") %>%
  inner_join(df_cen, by = "sample_hap") %>%
  filter(te_start >= cen_start & te_end <= cen_end) %>%
  inner_join(stv_min_hor_start, by = c("sample_hap_chrom")) %>%
    mutate(
    modstart = if_else(min_hor_start < cen_start, te_start - cen_start + 1, te_start - min_hor_start + 1),
    modend = if_else(min_hor_start < cen_start, te_end- cen_start + 1, te_end - min_hor_start + 1)
  ) %>%
  inner_join(left_updated_order, by = c("node" = "left_node"))
head(left_te_df,10)



stv_max_hor_end <- graphstv %>%
  group_by(sample_hap_chrom) %>%
  filter(horclass_color != "#C3C3C3") %>%
  summarise(max_hor_end = max(horend, na.rm = TRUE)) %>%
  ungroup()
head(stv_max_hor_end,10)


right_graph_plot_stv_df <- graphstv %>%
  inner_join(right_df, by = c("sample_hap_chrom")) %>%
  filter(horstart >= cen_start & horend <= cen_end) %>%
  inner_join(stv_max_hor_end, by = c("sample_hap_chrom")) %>%
  mutate(
    modstart = if_else(max_hor_end > cen_end, horstart - cen_end + 1, horstart - max_hor_end + 1),
    modend = if_else(max_hor_end > cen_end, horend - cen_end + 1, horend - max_hor_end + 1)
  ) %>%
  inner_join(right_updated_order, by = c("right_node" = "right_node"))
head(right_graph_plot_stv_df,10)
print(length(unique(right_graph_plot_stv_df$sample_hap_chrom)))

right_graph_te_df <- asat_te_df %>% 
  inner_join(node_df, by = "sample_hap") %>%
  inner_join(df_cen, by = "sample_hap") %>%
  filter(te_start >= cen_start & te_end <= cen_end) %>%
  inner_join(stv_max_hor_end, by = c("sample_hap_chrom")) %>%
  mutate(
    modstart = if_else(max_hor_end > cen_end, te_start - cen_end + 1, te_start - max_hor_end + 1),
    modend = if_else(max_hor_end > cen_end, te_end - cen_end + 1, te_end - max_hor_end + 1)
  ) %>%
  inner_join(right_updated_order, by = c("node" = "right_node"))
head(right_te_df,10)


max_len <- max(left_graph_plot_stv_df$modend)
y_labels <- left_graph_plot_stv_df %>%distinct(sample_hap_chrom, left_y)

###plot left ###
p3 <- ggplot() +
  geom_rect(
    data = left_graph_plot_stv_df,
    aes(xmin = modstart, xmax = modend, ymin = left_y - 0.3, ymax =left_y + 0.3, fill = horstv_color),
    color = NA
  ) +
  scale_fill_identity() +
  new_scale("fill") +  
   geom_point(
    data = left_graph_te_df,
    aes(
      x = (modstart + modend) / 2, 
      y = left_y + 0.4,  
      fill = te_color,    
      color = te_color 
    ),
    shape = 25,  
    size = 0.1,
    stroke = 0   
  ) +
  scale_color_identity() + 
  scale_fill_identity() +
  scale_x_continuous(
      breaks = seq(0, max_len, by = 2000000),
      labels = function(x) x / 1e6 #Mb
    ) +
  scale_y_continuous(breaks = y_labels$left_y, labels = y_labels$sample_hap_chrom) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 1),
    axis.title.y = element_text(size = 5)
  )

print(p3)

###plot right ###
if (chr == "chr10") {
  right_graph_plot_stv_df <- right_graph_plot_stv_df %>%
    filter(horstart <= max_hor_end & horend <= max_hor_end)
}

max_len <- min(right_graph_plot_stv_df$modend)
y_labels <- right_graph_plot_stv_df %>%distinct(sample_hap_chrom, right_y)


p4 <- ggplot() +
  geom_rect(
    data = right_graph_plot_stv_df,
    aes(xmin = modstart, xmax = modend, ymin = right_y - 0.3, ymax =right_y + 0.3, fill = horstv_color),
    color = NA
  ) +
  scale_fill_identity() +
  new_scale("fill") +  
   geom_point(
    data = right_graph_te_df,
    aes(
      x = (modstart + modend) / 2, 
      y = right_y + 0.4,  
      fill = te_color,    
      color = te_color 
    ),
    shape = 25,  
    size = 0.1,
    stroke = 0   
  ) +
  scale_color_identity() + 
  scale_fill_identity() +
  scale_x_continuous(
      breaks = seq(ceiling(max_len / 1e6) * 1e6, 0, by = 2e6),
      labels = function(x) x / 1e6 #Mb
    ) +
  scale_y_continuous(breaks = y_labels$right_y, labels = y_labels$sample_hap_chrom,
  position = "right") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 1),
    axis.title.y = element_text(size = 5)
  )

print(p4)

library(cowplot)
combined_plot <- plot_grid(p3, line_plot, p4, ncol = 3, rel_widths = c(1.5, 1, 1.5))
print(combined_plot)

outpdf <- paste0("Phy/", chr, "/", chr, "_phy_HORmon_cenhap.pdf")
ggsave(outpdf, plot = combined_plot, device = "pdf", width = 6, height = 10)
