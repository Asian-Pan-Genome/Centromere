setwd("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap")
rm(list = ls())

library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")

chr <- "chr17"
cenhap_file <- paste0("cenhap_classification/", chr, "_cenhap.xls")
print(cenhap_file)

if (chr=="chr8"){
  cenhap_df <- read.csv(cenhap_file, sep="\t", header=TRUE)
  colnames(cenhap_df) <- c("sample_hap", "Final_cluster", "CenHaps")
  cenhap_df <- cenhap_df %>%
      select(sample_hap, CenHaps)

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
}

if (chr=="chr4"){
  cenhap_df <- read.csv(cenhap_file, sep="\t", header=TRUE)
  colnames(cenhap_df) <- c("sample_hap", "Final_cluster", "CenHaps")
  cenhap_df <- cenhap_df %>%
      select(sample_hap, CenHaps)

  cenhap_color <- c(
  "A1" = "#8992BF",
  "A2" = "#81539E",
  "A3" = "#7A1C75",
  "B1" = "#C0CB64",
  "B2" = "#77C4BC",
  "B3" = "#094022",
  "B4" = "#20924F"
)
}

if (chr=="chr5"){
  cenhap_df <- read.csv(cenhap_file, sep="\t", header=TRUE)
  colnames(cenhap_df) <- c("sample_hap", "Final_cluster", "CenHaps")
  cenhap_df <- cenhap_df %>%
      select(sample_hap, CenHaps)

  cenhap_color <- c(
  "A" = "#BFA59B",
  "B" = "#8FADBF",
  "C" = "#797FF2",
  "D" = "#8C3545",
  "E" = "#F2B705",
  "F" = "#034AA6"
)
}

if (chr=="chr10"){
  cenhap_df <- read.csv(cenhap_file, sep="\t", header=TRUE)
  colnames(cenhap_df) <- c("sample_hap", "Final_cluster", "CenHaps")
  cenhap_df <- cenhap_df %>%
      select(sample_hap, CenHaps)

  cenhap_color <- c(
  "A1" = "#C57E7A",
  "A2" = "#CDAF83",
  "A3" = "#9C7465",
  "A4" = "#D0AFA5",
  "B1" = "#95B596",
  "B2" = "#768F7A",
  "B3" = "#AFCD99",
  "B4" = "#6AB384",
  "B5" = "#769D72",
  "C1" = "#8FA9C0",
  "C2" = "#5A7FAD"  
)
}


if (chr=="chr17"){
  cenhap_df <- read.csv(cenhap_file, sep="\t", header=TRUE)
  colnames(cenhap_df) <- c("sample_hap", "Final_cluster", "CenHaps", "insertion")
  cenhap_df <- cenhap_df %>%
      select(sample_hap, CenHaps, insertion)

  cenhap_color <- c(
  "A1" = "#EAD1DC",
  "A2" = "#C27BA0",
  "A3" = "#C90076",
  "B1" = "#6FA8DC",
  "B2" = "#146EA5"
)
}
head(cenhap_df)
dim(cenhap_df)

###basic info###
df_cen <- read.csv("cent_chrom.txt", header=TRUE, sep="\t")
df_cen <- df_cen %>% filter(chrom == chr)
head(df_cen,10)
complete_df_cen <- df_cen %>% filter(filterflag == 0)
print(dim(complete_df_cen))

df_pop <- read.csv("populaion.xls", header=TRUE, sep="\t")
df_stat <- df_cen  %>% 
    left_join(df_pop %>% select(sample, superpopulation), by = "sample") %>%
    mutate(sample_hap = str_replace(sample_hap, "CHM13_-", "CHM13")) %>%
    mutate(superpopulation = ifelse(superpopulation == "EAS-APG", "EAS", superpopulation)) %>%
    filter(filterflag == 0) %>%
    inner_join(cenhap_df, by = "sample_hap")
head(df_stat,10)
dim(df_stat)


####superpopulation ###
cenhap_pop <- df_stat %>%
    group_by(CenHaps, superpopulation) %>%
    summarise(count = n(), .groups = "drop")
head(cenhap_pop)

cenhap_num <- cenhap_pop %>%
    group_by(CenHaps) %>%
    summarise(total = sum(count))
head(cenhap_num)

superpopulation_colors <- c(
    "AFR" = "#319b62", 
    "AMR" = "#939393", 
    "EAS" = "#d22e77", 
    "EUR" = "#0070c0", 
    "SAS" = "#893f8b")

p <- ggplot(cenhap_pop, aes(x = CenHaps, y = count, fill = superpopulation)) +
    scale_fill_manual(values = superpopulation_colors) +
    geom_bar(stat = "identity") +
    geom_text(data = cenhap_num, aes(label = total, y = total), vjust = -0.5, size = 4, fontface = "bold", color = "black") +
    labs(x = "CenHaps", y = "Count") +
    theme_classic()
print(p)

p <- ggplot(cenhap_pop, aes(x = CenHaps, y = count, fill = superpopulation)) +
    scale_fill_manual(values = superpopulation_colors) +
    geom_bar(stat = "identity") +
    geom_text(data = cenhap_num, aes(x = CenHaps, y = total, label = total), 
              inherit.aes = FALSE,  # Prevent inheriting aesthetics from ggplot()
              vjust = -0.5, size = 4, color = "black") +
    labs(x = "CenHaps", y = "Count") +
    theme_classic()

print(p)
outpdf <- paste0("cenhap_classification/", chr, ".superpopulation.barplot.pdf")
print(outpdf)
ggsave(outpdf, plot = p, width = 6, height = 3)

### for chr17 superpopulation ###
cenhap_pop <- df_stat %>%
    group_by(CenHaps, insertion, superpopulation) %>%
    summarise(count = n(), .groups = "drop")
print(cenhap_pop, n=50)

cenhap_num <- cenhap_pop %>%
    group_by(CenHaps) %>%
    summarise(total = sum(count))
head(cenhap_num)


# Create a new variable in the data for the combination of superpopulation and insertion
cenhap_pop$fill_group <- paste(cenhap_pop$superpopulation, cenhap_pop$insertion, sep = "_")
head(cenhap_pop)

superpopulation_colors <- c(
    "AFR_0" = "#319b62", 
    "AMR_0" = "#939393", 
    "EAS_0" = "#d22e77", 
    "EUR_0" = "#0070c0", 
    "SAS_0" = "#893f8b",
    "AFR_1" = "#aadfc3", 
    "AMR_1" = "#df9c9c", 
    "EAS_1" = "#e290b5", 
    "EAS_2" = "#e290b5",
    "EUR_1" = "#81b6db", 
    "SAS_1" = "#e490e7")


# Update the ggplot code
p <- ggplot(cenhap_pop, aes(x = CenHaps, y = count, fill = fill_group)) +
    scale_fill_manual(values = superpopulation_colors) +
    geom_bar(stat = "identity") +
    geom_text(data = cenhap_num, aes(x = CenHaps, y = total, label = total), 
              inherit.aes = FALSE,  # Prevent inheriting aesthetics from ggplot()
              vjust = -0.5, size = 4, color = "black") +
    labs(x = "CenHaps", y = "Count") +
    theme_classic()
print(p)
outpdf <- paste0("cenhap_classification/", chr, ".superpopulation.barplot.pdf")
print(outpdf)
ggsave(outpdf, plot = p, width = 6, height = 3)


###alpha satellite length##
satfile <- "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/satellite_length/all.sat.length.xls"
satdf <- read.csv(satfile, sep="\t", header=TRUE)
head(satdf)

asat_df <- satdf %>%
    filter(satellite == "ASat") %>%
    select(sample_hap_chrom, length)
head(asat_df)

df_stat <- df_stat %>%
    inner_join(asat_df, by="sample_hap_chrom")
head(df_stat)

p2 <- ggplot(df_stat, aes(x = CenHaps, y = length/1000000, color = CenHaps)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +  # Boxplot with transparency and no outliers
    scale_color_manual(values = cenhap_color) +  
    geom_jitter(color = "grey", alpha = 0.3, width = 0.2) +  # Add jittered points
    labs(x = "CenHaps", y = "Length") +
    theme_classic()

print(p2)
outpdf <- paste0("cenhap_classification/", chr, ".asat.boxplot.pdf")
print(outpdf)
ggsave(outpdf, plot = p2, width = 6, height = 3)


p2 <- ggplot(df_stat, aes(x = CenHaps, y = len/1000000, color = CenHaps)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +  # Boxplot with transparency and no outliers
    scale_color_manual(values = cenhap_color) +  
    geom_jitter(color = "grey", alpha = 0.3, width = 0.2) +  # Add jittered points
    labs(x = "CenHaps", y = "Length (Mbp)") +
    theme_classic()

print(p2)
outpdf <- paste0("cenhap_classification/", chr, ".cent.boxplot.pdf")
print(outpdf)
ggsave(outpdf, plot = p2, width = 6, height = 3)



###plot chr4 hsat1 insertion distribution###
insertion_file <- "/share/home/zhanglab/user/sunyanqing/vscode_scripts/cenhap/chr4_hsat1_insertion.xls"
insertion_df <- read.csv(insertion_file, sep="\t", header = TRUE)
head(insertion_df)

df_stat <- df_stat %>%
    inner_join(insertion_df, by="sample_hap")
head(df_stat)    

p3 <- ggplot(df_stat, aes(x = CenHaps, y = HSat1A_ASat, color = CenHaps)) +
    geom_boxplot(aes(color = CenHaps)) +  # Add boxplot with color and fill
    scale_color_manual(values = cenhap_color) +  # Set fill colors
    scale_color_manual(values = cenhap_color) +  # Set border colors
    labs(x = "CenHaps", y = "HSat1A_ASat") +  # Add axis labels
    scale_y_continuous(breaks = seq(0, 10, by = 2)) +
    theme_classic()
# Print and save the plot

p3 <- ggplot(df_stat, aes(x = CenHaps, y = HSat1A_ASat, color = CenHaps)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +  # Boxplot with transparency and no outliers
    scale_color_manual(values = cenhap_color) +  
    geom_jitter(color = "grey", alpha = 0.3, width = 0.2) +  # Add jittered points
    labs(x = "CenHaps", y = "Number of HSat arrays") +
    scale_y_continuous(breaks = seq(0, 10, by = 2)) +
    theme_classic()
print(p3)
outpdf <- paste0("cenhap_classification/", chr, ".hsat1.boxplot.pdf")
print(outpdf)
ggsave(outpdf, plot = p3, width = 6, height = 6)
