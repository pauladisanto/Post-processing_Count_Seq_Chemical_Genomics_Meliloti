library(tidyverse)

# Set input and output directories
base_dir <- "/home/paula-di-santo/Documents/Post_analysis_meliloti"
deseq_dir <- file.path(base_dir, "DESeq2")
output_dir <- file.path(base_dir, "Normalized_count_plots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set the pools and page layout
pools <- paste0("POOL", 1:15)
features_per_page <- 6L
 
# Read sample metadata once for all pools
metadata <- read_csv(
  file.path(deseq_dir, "metadata_deseq2.csv"),
  show_col_types = FALSE
) %>%
  select(sample, Time, Compound, CultureID)

stopifnot(!anyDuplicated(metadata$sample))

plot_pool <- function(pool) {
  # Read all normalized counts and preserve sample names
  counts_pool <- read.csv(
    file.path(deseq_dir, paste0("normalized_counts_", pool, ".csv")),
    row.names = 1,
    check.names = FALSE
  )

  missing_samples <- setdiff(colnames(counts_pool), metadata$sample)
  if (length(missing_samples) > 0L) {
    stop("Samples missing from metadata in ", pool, ": ",
         paste(missing_samples, collapse = ", "))
  }

  # Keep T1–T3 and preserve the feature order from the count matrix
  features <- rownames(counts_pool)
  plot_data <- counts_pool %>%
    rownames_to_column("feature_ID") %>%
    pivot_longer(
      cols = -feature_ID,
      names_to = "sample",
      values_to = "normalized_count"
    ) %>%
    left_join(metadata, by = "sample") %>%
    filter(as.character(Time) %in% c("1", "2", "3")) %>%
    mutate(
      hours = 24 * as.integer(as.character(Time)),
      treatment = recode(as.character(Compound), C1 = "MMV688179"),
      feature_ID = factor(feature_ID, levels = features)
    )

  # Calculate mean and SD across cultures at each time point
  plot_summary <- plot_data %>%
    group_by(feature_ID, treatment, hours) %>%
    summarise(
      mean_count = mean(normalized_count),
      sd_count = sd(normalized_count),
      .groups = "drop"
    )

  # Split the feature panels into pages
  pages <- split(
    features,
    ceiling(seq_along(features) / features_per_page)
  )

  pdf_path <- file.path(
    output_dir,
    paste0(pool, "_all_features_mean_sd.pdf")
  )

  grDevices::pdf(pdf_path, width = 11, height = 7.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (page_number in seq_along(pages)) {
    page_features <- pages[[page_number]]

    p <- ggplot(
      filter(plot_summary, feature_ID %in% page_features),
      aes(x = hours, y = mean_count, color = treatment, group = treatment)
    ) +
      geom_point(
        data = filter(plot_data, feature_ID %in% page_features),
        aes(y = normalized_count),
        alpha = 0.35,
        size = 1.8
      ) +
      geom_line(linewidth = 0.8) +
      geom_errorbar(
        aes(ymin = mean_count - sd_count, ymax = mean_count + sd_count),
        width = 2,
        linewidth = 0.6
      ) +
      geom_point(size = 2.5) +
      facet_wrap(vars(feature_ID), ncol = 3, scales = "free_y") +
      scale_x_continuous(breaks = c(24, 48, 72)) +
      scale_color_manual(
        values = c(DMSO = "#0072B2", MMV688179 = "#D55E00")
      ) +
      labs(
        title = paste0(pool, " | Page ", page_number, " of ", length(pages)),
        x = "Time after treatment (h)",
        y = "DESeq2-normalized counts",
        color = "Treatment",
        caption = paste(
          "Lines connect means. Error bars show one SD above and below the mean.",
          "Faint points show individual cultures."
        )
      ) +
      theme_bw() +
      theme(legend.position = "bottom")

    print(p)
  }

  message(pool, ": ", length(features), " features, ", length(pages), " pages")
}

for (pool in pools) {
  plot_pool(pool)
}