setwd("D:/桌面/孙砚砚/科研/human-pan/analysis/Part3/figure/chr8")
rm(list = ls())


####plot_map_of cenhap ####
library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("scatterpie")

map <- map_data('world')
p <- ggplot() +
  geom_polygon(data = map, aes(x = long, y = lat, group = group), fill = 'gray80') + 
  theme_void() + 
  theme(panel.grid.minor = element_blank(), legend.background = element_blank()) 

print(p)

df <- read.csv("chr8_cenhap_map.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
head(df)

df_summary <- df %>%
  group_by(manual_classification, longitude, latitude) %>%
  summarise(total_count = n_distinct(sample_hap))
head(df_summary)

plot_df <- df %>%
  group_by(manual_classification, Final_subcluster, longitude, latitude) %>%
  summarise(count = n_distinct(sample_hap)) %>%
  left_join(df_summary, by=c("manual_classification", "longitude", "latitude")) 
head(plot_df)

plot_df_wide <- plot_df %>%
  pivot_wider(names_from = Final_subcluster, values_from = count, values_fill = 0)
  
head(plot_df_wide)


cenhap_color <- c(
  "A1" = "#1b7837",
  "A2" = "#00441b",
  "A3" = "#5aae61",
  "A4" = "#a6dba0",
  "A5" = "#d9f0d3",
  "B1" = "#40004b",
  "B2" = "#762a83",
  "B3" = "#965784",
  "B4" = "#9970ab",
  "B5" = "#FFD1FF",
  "C1" = "#bdbdbd",
  "C2" = "#737373"
)



####final version
cenhap_color <- c(
  "A1" = "#fdae61",
  "A2" = "#d73027",
  "A3" = "#f46d43",
  "A4" = "#a50026",
  "A5" = "#fee08b",
  "B1" = "#d9ef8b",
  "B2" = "#a6d96a",
  "B3" = "#66bd63",
  "B4" = "#1a9850",
  "B5" = "#006837",
  "C1" = "#6AAED6",
  "C2" = "#084A91"
)


cenhap_color <- c(
  "A1" = "#E8B522",
  "A2" = "#E6641C",
  "A3" = "#F4AA63",
  "A4" = "#a50026",
  "A5" = "#fee08b",
  "B1" = "#7DC6BA",
  "B2" = "#085B6F",
  "B3" = "#1AA4B4",
  "B4" = "#0554F2",
  "B5" = "#2793F2",
  "C1" = "#F280BF",
  "C2" = "#F230BF"
)


cenhap_color <- c(
  "A1" = "#E8B522",
  "A2" = "#E6641C",
  "A3" = "#F4AA63",
  "B1" = "#7DC6BA",
  "B2" = "#085B6F",
  "B3" = "#1AA4B4",
  "C1" = "#F280BF",
  "C2" = "#F230BF"
)


p <- ggplot() +
  geom_polygon(data = map, aes(x = long, y = lat, group = group), fill = 'gray80') + 
  theme_void() + 
  theme(panel.grid.minor = element_blank(), legend.background = element_blank()) +
  geom_scatterpie(data = plot_df_wide, 
                  aes(x = longitude, y = latitude, r = log(total_count+100)), 
                  cols = c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2"),
                  alpha = 0.7, color = NA) + 
  scale_fill_manual(values = cenhap_color) + 
  #geom_text(data = plot_df_wide, 
  #          aes(x = longitude, y = latitude + 0.5, label = manual_classification ), size = 1) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'cenhap')+ coord_equal()
print(p)
ggsave(p, filename = "chr8_cenhap_map_v5.pdf", width = 6, height = 8)

