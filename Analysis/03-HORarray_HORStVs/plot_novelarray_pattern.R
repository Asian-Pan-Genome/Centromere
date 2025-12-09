library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("patchwork")
library(cowplot)
library("ggplotify")
library(RColorBrewer)

info  <- read.csv("/share/home/zhanglab/user/sunyanqing/vscode_scripts/ancient_hor_array/Novel_horclass_info.xls", header=TRUE, sep="\t")
info <- info %>%
    mutate(horclass = paste0(horclass, "(", HORarray_reported_by_CHM13, ")"))%>%
    select(-HORarray_reported_by_CHM13)
head(info)

##count novel_array on chrom##
new_info <- info %>%
  mutate(chrom_manual_check_split = strsplit(as.character(chrom_manual_check), "_")) %>%
  unnest(chrom_manual_check_split) 
head(new_info,20)

novel_chrom_count <- new_info %>%
  group_by(chrom_manual_check_split) %>%
  summarise(unique_horclass_count = n_distinct(horclass)) %>%
  arrange(desc(unique_horclass_count))

head(novel_chrom_count, n=25)
###############################


df <- read.csv("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/horclass_length/novelarray_sample_count_mean_length.xls",  header=TRUE, sep="\t")
head(df)
print(df$length %>% median())
print(df$length %>% mean())
print(df$unique_sample_hap_count %>% median())

plot_df <- df %>%
    left_join(info, by = "horclass") 
head(plot_df)

print(plot_df$Max_repeat_number %>% min())
print(plot_df$Max_repeat_number %>% max())



# p <- ggplot(plot_df, aes(x = unique_sample_hap_count, y = length, size = HORstv, color = SF_hor)) +
#     geom_point() +
#     scale_size(range = c(3, 10)) +
#     labs(
#         x = "Unique Sample Hap Count",
#         y = "Length"       
#     ) +
#     theme_minimal()
# print(p)


SFcolor <- c(
    "Mixed" = "#000000",
    "SF01" =  "#1b9e77",
    "SF1" = "#7fc97f", 
    "SF4" = "#ffff99",
    "SF5" = "#386cb0"
)


p <- ggplot(plot_df, aes(x = unique_sample_hap_count, y = length/1000, size = Max_repeat_number, fill = SF_hor, color=SF_hor)) +
    geom_point(alpha=0.8, shape=21, stroke = 0) +
    scale_color_manual(values = SFcolor) +
    scale_fill_manual(values = SFcolor) + 
    # scale_size_continuous(
    #     range = c(1, 10),  # Set the range for point sizes (minimum and maximum size)
    #     breaks = c(5,10,15,20,25),  # Define the legend breaks
    #     name = "Max_repeat_number"  # Label for the size legend
    # ) +
    scale_x_continuous(
    breaks = c(0, 100, 200, 300, 400, 500)  # Set x-axis breaks
    ) +
    labs(
        x = "Number of haploid assemblies",
        y = "Mean length of HOR array (Kb)"
    ) +
    theme(
        panel.background = element_blank(),  # Remove background color
        panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add black border
        plot.background = element_blank(),  # Remove plot background color
        axis.line = element_line(color = "black")  # Add black axis lines
    )
print(p)
ggsave("/share/home/zhanglab/user/sunyanqing/vscode_scripts/ancient_hor_array/Novelarray_pattern.svg", plot = p, device = "svg", width = 8, height = 4)



library(ggExtra)

# Main scatter plot
p <- ggplot(plot_df, aes(x = unique_sample_hap_count, y = length / 1000, size = Max_repeat_number, fill = SF_hor, color = SF_hor)) +
    geom_point(alpha = 0.8, shape = 21, stroke = 0) +
    scale_color_manual(values = SFcolor) +
    scale_fill_manual(values = SFcolor) + 
    scale_size_continuous(
        range = c(1, 5),  # Set the range for point sizes (minimum and maximum size)
        breaks = c(5, 10, 15, 20, 25),  # Define the legend breaks
        name = "Max_repeat_number"  # Label for the size legend
    ) +
    scale_x_continuous(
        breaks = c(0, 100, 200, 300, 400, 500)  # Set x-axis breaks
    ) +
    labs(
        x = "Number of haploid assemblies",
        y = "Mean length of HOR array (Kb)"
    ) +
    theme(
        panel.background = element_blank(),  # Remove background color
        panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add black border
        plot.background = element_blank(),  # Remove plot background color
        axis.line = element_line(color = "black")  # Add black axis lines
    ) +
    geom_vline(xintercept = mean(plot_df$unique_sample_hap_count, na.rm = TRUE), linetype = "dashed", color = "blue") +  # Vertical line for mean x
    geom_hline(yintercept = mean(plot_df$length / 1000, na.rm = TRUE), linetype = "dashed", color = "blue")  # Horizontal line for mean y
print(p)
# Add marginal density plots
p_with_density <- ggMarginal(
    p,
    type = "density",  # Add density plots
    margins = "both",  # Add to both x and y axes
    size = 5,  # Adjust the size of the marginal plots
    colour = "black",  # Outline color for density plots
    fill = "#c7eae5",
)

# Print the plot
print(p_with_density)
ggsave("/share/home/zhanglab/user/sunyanqing/vscode_scripts/ancient_hor_array/Novelarray_pattern_update.svg", 
plot = p_with_density, device = "svg", width = 8, height = 4)
