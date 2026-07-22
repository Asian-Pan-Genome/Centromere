suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

## load data ----------------------------------------------------------------
satlen <- read.csv("data/all.sat.length.xls", header = TRUE, sep = "\t")
df_cen <- read.csv("data/cent_chrom.txt", header = TRUE, sep = "\t")

## centromere filter flag (per sample_hap_chrom) -----------------------------
cen_flag <- df_cen %>%
  group_by(sample_hap_chrom) %>%
  summarise(filterflag = max(filterflag, na.rm = TRUE), .groups = "drop")

## HSat1 = HSat1A + HSat1B ---------------------------------------------------
hsat1 <- satlen %>%
  filter(project == "APG", satellite %in% c("HSat1A", "HSat1B")) %>%
  group_by(sample_hap_chrom, sample_hap, sample, hap, chrom, project) %>%
  summarise(length = sum(length, na.rm = TRUE), .groups = "drop") %>%
  mutate(satellite = "HSat1")

## HSat2 & HSat3 -------------------------------------------------------------
hsat23 <- satlen %>%
  filter(project == "APG", satellite %in% c("HSat2", "HSat3"))

## combine -------------------------------------------------------------------
df <- bind_rows(hsat1, hsat23) %>%
  left_join(cen_flag, by = "sample_hap_chrom")

## pivot to Mat/Pat, keep pairs where both filterflag == 0 -------------------
process_satellite <- function(data, sat_name) {
  pair <- data %>%
    filter(satellite == sat_name) %>%
    select(sample, chrom, satellite, hap, length, filterflag) %>%
    pivot_wider(names_from = hap, values_from = c(length, filterflag)) %>%
    filter(!is.na(length_Mat), !is.na(length_Pat),
           filterflag_Mat == 0, filterflag_Pat == 0) %>%
    mutate(diff_size = length_Mat - length_Pat)
  return(pair)
}

hsat1_pair <- process_satellite(df, "HSat1")
hsat2_pair <- process_satellite(df, "HSat2")
hsat3_pair <- process_satellite(df, "HSat3")

## helper: nice symmetric breaks from xlim -----------------------------------
nice_breaks <- function(limit) {
  step <- limit / 2  # one step per half-range, e.g. ±16 -> step=8
  seq(-limit, limit, by = step)
}

## density + histogram plot function -----------------------------------------
plot_signed <- function(pair_data, sat_label, x_limit, binwidth = 0.5,
                        plot_type = c("density", "histogram")) {
  plot_type <- match.arg(plot_type)

  mean_signed_mb <- mean(pair_data$diff_size, na.rm = TRUE) / 1e6
  sd_signed_mb   <- sd(pair_data$diff_size,   na.rm = TRUE) / 1e6
  lower_2sd_mb   <- mean_signed_mb - 2 * sd_signed_mb
  upper_2sd_mb   <- mean_signed_mb + 2 * sd_signed_mb

  cat("\n=== ", sat_label, " (", plot_type, ") ===\n")
  cat("n pairs:", nrow(pair_data), "\n")
  cat("Signed (Mat - Pat) mean (Mb):", mean_signed_mb, "\n")
  cat("Signed (Mat - Pat) sd   (Mb):", sd_signed_mb,   "\n")
  cat("mean - 2SD (Mb):", lower_2sd_mb, "\n")
  cat("mean + 2SD (Mb):", upper_2sd_mb, "\n")

  breaks <- nice_breaks(x_limit)

  p <- ggplot(pair_data, aes(x = diff_size / 1e6))

  if (plot_type == "density") {
    p <- p + geom_density(fill = "grey", alpha = 0.5)
  } else {
    p <- p + geom_histogram(fill = "grey", alpha = 0.5,
                            binwidth = binwidth, boundary = 0)
  }

  p <- p +
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
    scale_x_continuous(breaks = breaks, limits = c(-x_limit, x_limit)) +
    theme_classic(base_family = "Arial") +
    labs(
      title = "",
      x = paste0(sat_label, " length difference (Mat - Pat, Mbp)"),
      y = if (plot_type == "density") "Density" else "Count"
    ) +
    theme(
      axis.title  = element_text(size = 12, color = "black"),
      axis.text   = element_text(size = 10, color = "black"),
      axis.ticks  = element_line(color = "black"),
      axis.line   = element_line(color = "black"),
      plot.margin = margin(5.5, 20, 5.5, 20)
    )

  filename <- paste0("apg_mat_vs_pat_", tolower(sat_label), "_signed_",
                     plot_type, ".pdf")
  ggsave(p, filename = filename, device = cairo_pdf, width = 6, height = 5)
  cat("Saved:", filename, "\n")
  return(p)
}

## config: satellite -> (x_limit, binwidth) -----------------------------------
config <- list(
  HSat1 = list(pair = hsat1_pair, xlim = 8.5,  bw = 0.5),
  HSat2 = list(pair = hsat2_pair, xlim = 16,   bw = 1.0),
  HSat3 = list(pair = hsat3_pair, xlim = 14,   bw = 1.0)
)

for (sat in names(config)) {
  cfg <- config[[sat]]
  plot_signed(cfg$pair, sat, cfg$xlim, cfg$bw, "density")
  plot_signed(cfg$pair, sat, cfg$xlim, cfg$bw, "histogram")
}

cat("\n=== All done ===\n")
