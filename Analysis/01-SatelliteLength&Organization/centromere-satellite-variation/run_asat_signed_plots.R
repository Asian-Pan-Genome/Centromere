suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

## load data (mirrors plot_mat_vs_pat_satellite_length.R prep) ---------------
satlen <- read.csv("data/all.sat.length.xls", header = TRUE, sep = "\t")
df_cen <- read.csv("data/cent_chrom.txt", header = TRUE, sep = "\t")

cen_flag <- df_cen %>%
  group_by(sample_hap_chrom) %>%
  summarise(filterflag = max(filterflag, na.rm = TRUE), .groups = "drop")

df <- satlen %>% left_join(cen_flag, by = "sample_hap_chrom") %>%
  filter(project == "APG", satellite == "ASat")

## pivot to Mat/Pat, keep pairs where both filterflag == 0 -------------------
asat_pair <- df %>%
  select(sample, chrom, satellite, hap, length, filterflag) %>%
  pivot_wider(names_from = hap, values_from = c(length, filterflag)) %>%
  filter(!is.na(length_Mat), !is.na(length_Pat),
         filterflag_Mat == 0, filterflag_Pat == 0) %>%
  mutate(diff_size = length_Mat - length_Pat)

## stats ---------------------------------------------------------------------
mean_signed_mb <- mean(asat_pair$diff_size, na.rm = TRUE) / 1e6
sd_signed_mb   <- sd(asat_pair$diff_size,   na.rm = TRUE) / 1e6
lower_2sd_mb   <- mean_signed_mb - 2 * sd_signed_mb
upper_2sd_mb   <- mean_signed_mb + 2 * sd_signed_mb

cat("n ASat pairs:", nrow(asat_pair), "\n")
cat("ASat signed (Mat - Pat) mean (Mb):", mean_signed_mb, "\n")
cat("ASat signed (Mat - Pat) sd   (Mb):", sd_signed_mb,   "\n")
cat("mean - 2SD (Mb):", lower_2sd_mb, "\n")
cat("mean + 2SD (Mb):", upper_2sd_mb, "\n")

## shared plot builder -------------------------------------------------------
build_asat_plot <- function(plot_type = c("density", "histogram")) {
  plot_type <- match.arg(plot_type)

  p <- ggplot(asat_pair, aes(x = diff_size / 1e6))

  if (plot_type == "density") {
    p <- p + geom_density(fill = "grey", alpha = 0.5)
  } else {
    p <- p + geom_histogram(fill = "grey", alpha = 0.5,
                            binwidth = 0.5, boundary = 0)
  }

  p +
    geom_vline(xintercept = lower_2sd_mb, linetype = "dashed", color = "black",
               linewidth = 0.6) +
    geom_vline(xintercept = upper_2sd_mb, linetype = "dashed", color = "black",
               linewidth = 0.6) +
    annotate("text", x = lower_2sd_mb, y = Inf,
             label = sprintf("%.2f(mean - 2s.d.)", lower_2sd_mb),
             hjust = 1.05, vjust = 1.5, size = 3.2, color = "black",
             family = "Arial") +
    annotate("text", x = upper_2sd_mb, y = Inf,
             label = sprintf("%.2f(mean + 2s.d.)", upper_2sd_mb),
             hjust = -0.05, vjust = 1.5, size = 3.2, color = "black",
             family = "Arial") +
    scale_x_continuous(breaks = seq(-8, 8, 4), limits = c(-8.5, 8.5)) +
    theme_classic(base_family = "Arial") +
    labs(
      title = "",
      x = "αSat length difference (Mat - Pat, Mbp)",
      y = if (plot_type == "density") "Density" else "Count"
    ) +
    theme(
      axis.title  = element_text(size = 12, color = "black"),
      axis.text   = element_text(size = 10, color = "black"),
      axis.ticks  = element_line(color = "black"),
      axis.line   = element_line(color = "black"),
      plot.margin = margin(5.5, 20, 5.5, 20)
    )
}

## density plot --------------------------------------------------------------
ggsave(build_asat_plot("density"),
       filename = "apg_mat_vs_pat_asat_length_signed_density.pdf",
       device = cairo_pdf, width = 6, height = 5)
cat("Saved: apg_mat_vs_pat_asat_length_signed_density.pdf\n")

## histogram plot ------------------------------------------------------------
ggsave(build_asat_plot("histogram"),
       filename = "apg_mat_vs_pat_asat_length_signed_histogram.pdf",
       device = cairo_pdf, width = 6, height = 5)
cat("Saved: apg_mat_vs_pat_asat_length_signed_histogram.pdf\n")

cat("\n=== All done ===\n")
