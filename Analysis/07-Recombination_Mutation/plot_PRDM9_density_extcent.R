library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library(ggnewscale)

setwd("/share/home/zhanglab/user/sunyanqing/human/periphy/CHM13/PRDM9")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript plot_PRDM9_density.R <sample> <hap>")
}

sample <- args[1]
hap <- args[2]
chr <- "chr19"

# Print the input values for confirmation
cat("Sample:", sample, "\n")
cat("Hap:", hap, "\n")

# sample="C045-CHA-N05"
# hap="Mat"

if (sample == "CHM13") {
  cenanno <- "/share/home/zhanglab/user/sunyanqing/human/anno/CHM13/CHM13v2.full.merged.cenanno.bed"
  
} else if (sample == "CN1") {
  cenanno <- paste0("/share/home/zhanglab/user/sunyanqing/human/anno/", sample, "/", hap, "/", sample, "_", hap, ".v1.0.full.merged.cenanno.bed")
} else if (sample == "HG002") {
  cenanno <- paste0("/share/home/zhanglab/user/sunyanqing/human/anno/", sample, "/", hap, "/", sample, "_", hap, ".v1.0.1.full.merged.cenanno.bed")
} else if (sample == "YAO") {
  cenanno <- paste0("/share/home/zhanglab/user/sunyanqing/human/anno/", sample, "/", hap, "/", sample, "_", hap, ".full.merged.cenanno.bed")
} else if (grepl("HG", sample) || grepl("NA", sample)) {
  cenanno <- paste0("/share/home/zhanglab/user/sunyanqing/human/anno/", sample, "/", hap, "/", sample, "_", hap, ".full.merged.cenanno.bed")
} else {
  cenanno <- paste0("/share/home/zhanglab/user/sunyanqing/human/anno/", sample, "/", hap, "/", sample, "_", hap, ".v0.9.full.merged.cenanno.bed")
}

if (sample == "CHM13"){
  horstv <- "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/HiCAT_bhlayer_update_split/CHM13_-.HiCAT.horstv.bed"
} else{
  horstv = paste0("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/HiCAT_bhlayer_update_split/", sample, "_", hap, ".HiCAT.horstv.bed")
}


if (sample == "CHM13"){
  infile <- "CHM13/chm13v2.20kb.PRDM9.bed"
}else{
  infile <- paste0(sample, "_", hap, "/", sample, "_", hap, ".20kb.PRDM9.bed")
}
print(infile)


df <- read.csv(infile, sep = "\t", header = FALSE)
colnames(df) <- c("chrom", "start", "end", "chr", "motif_start", "motif_end", "motif", "score", "strand", "p-value", "q-value", "seq", "length")
head(df)


te_color_df <- read.csv("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/te_color.xls", header=T, sep="\t")
head(te_color_df)

tefile <- paste0("/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/TE/", chr, "_te.bed")
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

# motif_stat_df <- df %>%
#   group_by(chrom, start, end, motif) %>%
#   summarise(
#     count = n()
#   )
# head(motif_stat_df)

all_stat_df <- df %>%
  group_by(chrom, start, end) %>%
  summarise(
    count = n(),
  ) %>% mutate("motif" = "all")%>%
  select(chrom, start, end, motif, count)
head(all_stat_df)



# combined_df <- bind_rows(motif_stat_df, all_stat_df)
# head(combined_df)


if(sample == "CHM13"){
  for (target_chrom in unique(combined_df$chrom)) {
    if (grepl("chr19", target_chrom)){
      plot_df <- all_stat_df %>% filter(chrom == target_chrom)

      all_windows <- tibble(
      chrom = unique(plot_df$chrom),
      start = seq(min(plot_df$start), max(plot_df$start), by = 20000),
      end =seq(min(plot_df$end), max(plot_df$end), by = 20000),
      )

      plot_prdm9_all <- all_windows %>%
        left_join(plot_df, by = c("chrom", "start", "end")) %>%
        mutate(count = replace_na(count, 0))
      write.csv(plot_prdm9_all, file = paste0("CHM13/PRDM9_count_", target_chrom, ".csv"), row.names = FALSE, quote = FALSE, sep = "\t")
      print(plot_prdm9_all, n = 188)
      unique_motifs <- unique(plot_prdm9_all$motif)


      peri_region_files <- paste0("/share/home/zhanglab/user/sunyanqing/human/periphy/CHM13/IBS/",target_chrom,"/100kb_windows_with_snp_counts_filtered.bed")
      peri_regions <- read.csv(peri_region_files, sep = "\t", header = FALSE)
      colnames(peri_regions) <- c("chrom", "start", "end", "index", "snp_count")
      head(peri_regions)

      peri_regions_expanded <- peri_regions %>%
        tidyr::crossing(motif = unique_motifs)
      head(peri_regions_expanded)

      p <- ggplot(plot_prdm9_all, aes(x = start, y = count, color = chrom)) +
        geom_line() +
        scale_x_continuous(
          breaks = seq(0, max(plot_prdm9_all$start), by = 5000000),
          labels = seq(0, max(plot_prdm9_all$start), by = 5000000) / 1000000
        ) +
        labs(
          x = "Start (Mb)",
          y = "Count",
          title = paste0("PRDM9 Density Across Chromosomes: ", target_chrom)
        ) +
        geom_rect(peri_regions, mapping = aes(xmin = start, xmax = end, ymin = 0, ymax = 600), fill = "black", alpha = 0.2, inherit.aes = FALSE) +
        facet_wrap(~motif, ncol = 1) +
        theme_minimal() +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
        )
      print(p)
      ggsave(p, filename = paste0("PRDM9_density_", target_chrom, ".png"), width = 6, height = 10, dpi=600)
    }
  }
}else{
  # target_chrom <- "C045-CHA-N05#Mat#chr19"
  for (target_chrom in unique(all_stat_df$chrom)) {
    if (grepl("chr19", target_chrom) && !grepl("chrY-CEN", target_chrom)){
      plot_df <- all_stat_df %>% filter(chrom == target_chrom)

      all_windows <- tibble(
        chrom = unique(plot_df$chrom),
        start = seq(min(plot_df$start), max(plot_df$start), by = 20000),
        end = seq(min(plot_df$start) + 20000, max(plot_df$start) + 20000, by = 20000),
      )


      plot_prdm9_all <- all_windows %>%
        left_join(plot_df, by = c("chrom", "start", "end")) %>%
        mutate(count = replace_na(count, 0))
      # write.csv(plot_prdm9_all, file = paste0(sample, "_", hap, "/", "PRDM9_count_", target_chrom, ".csv"), row.names = FALSE, quote = FALSE, sep = "\t")
      print(plot_prdm9_all, n = 188)
      unique_motifs <- unique(plot_prdm9_all$motif)

      cenanno_df <- read.csv(cenanno, sep = "\t", header = FALSE)
      colnames(cenanno_df) <- c("chrom", "start", "end", "satellite", "score", "strand", "s", "e", "color")
      head(cenanno_df)

      cent_regions <- cenanno_df %>% filter(chrom == target_chrom)
      head(cent_regions)

      ##get peri_cent_start
      ext_cent_start <- floor((min(cent_regions$start) - 1000000) / 20000) * 20000
      ext_cent_end <- ceiling((max(cent_regions$end) +  1000000) / 20000) * 20000

      cat(ext_cent_start, ext_cent_end)

      horstv_df <- read.csv(horstv, sep = "\t", header = FALSE)
      colnames(horstv_df) <- c("chrom", "start", "end", "hor", "horclass", "horclass_color", "horstv_index", "horstv_color")
      head(horstv_df)

      recolor_file <- "/share/home/zhanglab/user/sunyanqing/human/anno/statistics/centhap/all.horstv.hicat.update.color.txt"
      recolordf  <- read.csv(recolor_file, header=TRUE, sep="\t")
      colnames(recolordf) <- c("horstv_index", "horstv_newcolor")
      print(head(recolordf,10))

      horstv_regions <- horstv_df %>% filter(chrom == target_chrom) %>%
        left_join(recolordf, by = "horstv_index") %>%
        mutate(horstv_newcolor = ifelse(horstv_color == '-', "#525252", 
                                  ifelse(!is.na(horstv_newcolor), horstv_newcolor, horstv_color)))
      head(horstv_regions)

      ext_cent_prdm9 <- plot_prdm9_all %>%
        filter(start >= ext_cent_start, 
               end <= ext_cent_end)
      tail(ext_cent_prdm9)

      plot_te_df <- asat_te_df %>% 
        filter(sample_hap == paste0(sample, "_", hap),
               te_start >= ext_cent_start,
               te_end <= ext_cent_end) 
      head(plot_te_df)
      str(plot_te_df$te_color)

      p <- ggplot(ext_cent_prdm9, aes(x = start, y = count)) +
      geom_line()+
      scale_x_continuous(
        breaks = seq(0, max(ext_cent_prdm9$start), by = 5000000),
        labels = seq(0, max(ext_cent_prdm9$start), by = 5000000) / 1000000
      ) +
      labs(
        x = "Start (Mb)",
        y = "PRDM9 motif hits",
        title = paste0("PRDM9 Density Across Chromosomes: ", target_chrom)
      ) +
      geom_rect(cent_regions, mapping = aes(xmin = start, xmax = end, ymin = 95, ymax = 100,  fill = color), inherit.aes = FALSE) +
      scale_fill_identity() + 
      new_scale("fill") +  
      geom_rect(horstv_regions, mapping = aes(xmin = start, xmax = end, ymin = 85, ymax = 90, fill = horstv_newcolor), inherit.aes = FALSE) +
      scale_fill_identity() + 
      new_scale("fill") +  
      geom_point(
        data = plot_te_df,
        mapping = aes(
          x = (te_start + te_end) / 2, 
          y = 91,  
          fill = te_color,    
          color = te_color 
        ),
        shape = 25,  
        size = 1,
        stroke = 0,
        inherit.aes = FALSE  
      ) +
      scale_color_identity() + 
      scale_fill_identity() +
      # facet_wrap(~motif, ncol = 1) +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4),
        plot.title = element_text(size = 6)
      )
    print(p)
    ggsave(p, filename = paste0(sample, "_", hap, "/", "PRDM9_density_cent_", target_chrom, ".pdf"), width = 6, height = 3)
    }
  }
}






