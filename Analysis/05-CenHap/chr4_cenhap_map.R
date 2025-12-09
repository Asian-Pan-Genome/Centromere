setwd("D:/桌面/孙砚砚/科研/human-pan/analysis/Part3/figure/chr4")
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

df <- read.csv("chr4_cenhap_map.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
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
  "A1" = "#8992BF",
  "A2" = "#81539E",
  "A3" = "#7A1C75",
  "B1" = "#C0CB64",
  "B2" = "#77C4BC",
  "B3" = "#094022",
  "B4" = "#20924F"
)


p <- ggplot() +
  geom_polygon(data = map, aes(x = long, y = lat, group = group), fill = 'gray80') + 
  theme_void() + 
  theme(panel.grid.minor = element_blank(), legend.background = element_blank()) +
  geom_scatterpie(data = plot_df_wide, 
                  aes(x = longitude, y = latitude, r = log(total_count+100)), 
                  cols = c("A1", "A2", "A3", "B1", "B2", "B3", "B4"),
                  alpha = 0.7, color = NA) + 
  scale_fill_manual(values = cenhap_color) + 
  #geom_text(data = plot_df_wide, 
  #          aes(x = longitude, y = latitude + 0.5, label = manual_classification ), size = 1) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'cenhap')+ coord_equal()
print(p)
ggsave(p, filename = "chr4_cenhap_map_v3.pdf", width = 6, height = 8)
