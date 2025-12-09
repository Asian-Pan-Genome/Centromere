setwd("/share/home/zhanglab/user/sunyanqing/human/alpha/vsearch/20241022/step-11/cumsum")
getwd()
rm(list = ls())

library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library(patchwork)

##load hicat cumsum file##
# hicat_df <- read.csv("hicat_apg_first_cumsum_ordered_count.xls", header=FALSE, sep="\t")
hicat_df <- read.csv("hicat_hprc_first_cumsum_ordered_count_update_0907.xls", header=FALSE, sep="\t")
colnames(hicat_df) <- c("index", "sample_hap", "hicat_union_num", "hicat_shared_num")
head(hicat_df)

# graph_df <- read.csv("graph_apg_first_cumsum_ordered_count.xls", header=FALSE, sep="\t")
graph_df <- read.csv("graph_hprc_first_cumsum_ordered_count_update_0907.xls", header=FALSE, sep="\t")
colnames(graph_df) <- c("index", "sample_hap", "graph_union_num", "graph_shared_num")
head(graph_df)

merged_df <- merge(hicat_df, graph_df, by = c("index", "sample_hap"))
head(merged_df)


##load project info##
info_df <- read.csv("cumsum.order.xls", header=TRUE, sep="\t")
colnames(info_df) <- c("sample_hap", "sample", "population", "superpopulation", "project", "Order", "cumsum_Order", "ref_apg_hprc_hgsvc", "ref_hprc_hgsvc_apg")
head(info_df)

merged_df <- merge(merged_df, info_df %>% select(sample_hap,project), by="sample_hap")  %>% arrange(index)
head(merged_df)

custom_colors <- c(
  "APG" = "#E85827",
  "REF" = "#969696",
  "HPRC" = "#008C8C",
  "HGSVC" = "#285943"
)


long_df <- merged_df %>%
  pivot_longer(
    cols = c(hicat_union_num, hicat_shared_num, graph_union_num, graph_shared_num),
    names_to = c("type", "metric"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(
    source = ifelse(type == "hicat", "HiCAT", "HORmon"),
    shape = ifelse(metric == "union", 21, 24)
  )
print(head(long_df))

union_long_df <- long_df %>%
  filter(metric == "union") 


##plot##

p <- ggplot(union_long_df) +
  geom_point(aes(x = index, y = value, color = project, shape = factor(type)),
             size = 1, fill = NA, stroke = 0.2) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = c(21, 22)) +
  scale_y_continuous(limits=c(0,6000), breaks = seq(0, 6000, by = 1000)) +
  scale_x_continuous(breaks = seq(0, 500, by = 100)) +
  labs(
    x = "Number of shared haploid assemblies",
    y = "Discovered HORs",
    color = "Project",
    shape = "Method"
  ) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    
    # 设置边框线
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # 边框颜色和粗细
    axis.line = element_line(color = "black", size = 0.5),  # 轴线
    axis.ticks = element_line(color = "black", size = 0.5)  # 轴刻度线
  )
print(p)

outpdf <- "cumsum_curve_apg_first_hicat_graph.pdf"
outpdf <- "union_cumsum_curve_hprc_first_hicat_graph.pdf"
print(outpdf)
ggsave(outpdf, plot = p, device = "pdf", width = 8 , height = 6)


####load hicat cumsum with pantype file##
hicat_cumsum_pantype_df <- read.csv("hicat_hprc_first_cumsum_ordered_count_pantype.xls", header=FALSE, sep="\t")
colnames(hicat_cumsum_pantype_df) <- c("index", "sample_hap", "pan_hor_type", "hicat_union_num", "hicat_shared_num" )
head(hicat_cumsum_pantype_df)

hicat_cumsum_pantype_df$pan_hor_type <- factor(hicat_cumsum_pantype_df$pan_hor_type, levels = rev(c("Core", "Dispensable", "Low", "Rare", "Singleton")))
head(hicat_cumsum_pantype_df)

cumsum_order_df <- read.csv("cumsum.order.xls", header=TRUE, sep="\t")
head(cumsum_order_df)

cumsum_order_plot <- cumsum_order_df %>%
  select(project, ref_hprc_hgsvc_apg) %>%
  group_by(project) %>%
  summarize(
    min_value = min(ref_hprc_hgsvc_apg, na.rm = TRUE),
    max_value = max(ref_hprc_hgsvc_apg, na.rm = TRUE)
  )
head(cumsum_order_plot)

print(hicat_cumsum_pantype_df$pan_hor_type %>% unique())

# pantype_color <- c(
#     "Core" = "#7F7340",
#     "Dispensable" = "#C6B28D",
#     "Low" = "#BB8642",
#     "Rare" = "#914532",
#     "Singleton" = "#511F20"
# )

# pantype_color <- c(
#     "Core" = "#67001f",
#     "Dispensable" = "#b2182b",
#     "Low" = "#d6604d",
#     "Rare" = "#f4a582",
#     "Singleton" = "#fddbc7"
# )

pantype_color <- c(
    "Core" = "#313695",
    "Dispensable" = "#4575b4",
    "Low frequency" = "#74add1",
    "Rare" = "#abd9e9",
    "Singleton" = "#e0f3f8"
)

custom_colors <- c(
  "APG" = "#E85827",
  "REF" = "#969696",
  "HPRC" = "#008C8C",
  "HGSVC" = "#285943"
)

library(ggnewscale) 

p <- ggplot(hicat_cumsum_pantype_df, aes(x = index, y = hicat_union_num, fill = pan_hor_type)) +
  geom_bar(stat = "identity", position = "stack", alpha=0.8) +
  labs(x = "Number of haploid assemblies", y = "Discovered HORs") +
  scale_y_continuous(breaks = seq(0, 6000, by = 1000)) + 
  scale_x_continuous(breaks = seq(0, 500, by = 100)) +
  scale_fill_manual(values = pantype_color) +
  new_scale_fill() +
  geom_rect(data = cumsum_order_plot,
            aes(xmin = min_value, xmax = max_value, 
                ymin = 5995, ymax = 5999, 
                fill = project),
            inherit.aes = FALSE) +  
  scale_fill_manual(
    name = "Super Population",
    values = custom_colors) +
  theme_classic()
print(p)
ggsave("cumsum_hicat_pantype.pdf",plot = p, width = 8,height = 6)


######load graph cumsum with pantype file##
graph_cumsum_pantype_df <- read.csv("graph_hprc_first_cumsum_ordered_count_pantype.xls", header=FALSE, sep="\t")
colnames(graph_cumsum_pantype_df) <- c("index", "sample_hap", "pan_hor_type", "hicat_union_num", "hicat_shared_num" )
head(graph_cumsum_pantype_df)

graph_cumsum_pantype_df$pan_hor_type <- factor(hicat_cumsum_pantype_df$pan_hor_type, levels = rev(c("Core", "Dispensable", "Low frequency", "Rare", "Singleton")))
head(graph_cumsum_pantype_df)

cumsum_order_df <- read.csv("cumsum.order.xls", header=TRUE, sep="\t")
head(cumsum_order_df)

cumsum_order_plot <- cumsum_order_df %>%
  select(project, ref_hprc_hgsvc_apg) %>%
  group_by(project) %>%
  summarize(
    min_value = min(ref_hprc_hgsvc_apg, na.rm = TRUE),
    max_value = max(ref_hprc_hgsvc_apg, na.rm = TRUE)
  )
head(cumsum_order_plot)

print(graph_cumsum_pantype_df$pan_hor_type %>% unique())

p <- ggplot(graph_cumsum_pantype_df, aes(x = index, y = hicat_union_num, fill = pan_hor_type)) +
  geom_bar(stat = "identity", position = "stack", alpha=0.8) +
  labs(x = "Number of haploid assemblies", y = "Discovered HORs") +
  scale_y_continuous(breaks = seq(0, 5000, by = 1000)) + 
  scale_x_continuous(breaks = seq(0, 500, by = 100)) +
  scale_fill_manual(values = pantype_color) +
  new_scale_fill() +
  geom_rect(data = cumsum_order_plot,
            aes(xmin = min_value, xmax = max_value, 
                ymin = 4995, ymax = 4999, 
                fill = project),
            inherit.aes = FALSE) +  
  scale_fill_manual(
    name = "Super Population",
    values = custom_colors) +
  theme_classic()
print(p)
ggsave("cumsum_graph_pantype.pdf",plot = p, width = 8,height = 6)
