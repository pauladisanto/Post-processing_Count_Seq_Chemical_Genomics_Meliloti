library(tidyverse)
library(ggplot2)

# =========================
# Importing data
# =========================
clean_table_annotated <- read_csv("clean_table_annotated.csv")

coldata <- read_csv("metadata_deseq2.csv") %>%
  column_to_rownames("sample")

size_factor_df <- read_csv("size_factors.csv")

plot_counts <- clean_table_annotated %>%
  mutate(
    sample = as.character(sample),
    gene_ID = as.character(gene_ID),
    transposon_set = as.character(transposon_set)
  ) %>%
  filter(!is.na(gene_ID), !is.na(transposon_set)) %>%
  group_by(sample, gene_ID, transposon_set) %>%
  summarise(
    raw_reads = sum(NumReads, na.rm = TRUE),
    .groups = "drop"
  )

plot_data <- plot_counts %>%
  inner_join(
    tibble(sample = rownames(coldata), coldata),
    by = "sample"
  ) %>%
  inner_join(
    size_factor_df,
    by = "sample"
  ) %>%
  mutate(
    normalized_reads = raw_reads / size_factor,
    Time = factor(Time, levels = sort(unique(as.numeric(as.character(Time))))),
    Compound = factor(Compound),
    transposon_set = factor(transposon_set),
    Replicate = factor(Replicate)
  )


# quick QC: expected around 3 replicates per gene/pool/condition
replicate_check <- plot_data %>%
  dplyr::distinct(sample, gene_ID, transposon_set, Time, Compound, Replicate) %>%
  dplyr::group_by(gene_ID, transposon_set, Time, Compound) %>%
  dplyr::summarise(
    n_replicates = dplyr::n(),
    .groups = "drop"
  )

write_csv(replicate_check, "plot_replicate_check_by_pool.csv")

# =========================
# Summarize data for plotting
# =========================
summary_plot_data <- plot_data %>%
  dplyr::group_by(gene_ID, transposon_set, Time, Compound) %>%
  dplyr::summarise(
    mean_reads = mean(normalized_reads, na.rm = TRUE),
    sd_reads = sd(normalized_reads, na.rm = TRUE),
    n = dplyr::n(),
    se_reads = ifelse(n > 1, sd_reads / sqrt(n), NA_real_),
    .groups = "drop"
  )

write_csv(summary_plot_data, "plot_summary_all_genes_by_pool.csv")


# =========================
# Duplicate baseline (Time 0, None) into C1 and DMSO
# =========================

# make sure Time behaves numerically for sorting/comparison
summary_plot_data <- summary_plot_data %>%
  mutate(Time_num = as.numeric(as.character(Time)))

plot_summary_lines <- summary_plot_data %>%
  filter(Compound %in% c("C1", "DMSO")) %>%
  mutate(LineGroup = Compound)

baseline_rows <- summary_plot_data %>%
  filter(Time_num == 0) %>%
  tidyr::crossing(LineGroup = c("C1", "DMSO")) %>%
  mutate(Compound = LineGroup)

plot_summary_lines <- bind_rows(plot_summary_lines, baseline_rows) %>%
  mutate(
    Time = factor(Time_num, levels = sort(unique(Time_num))),
    Compound = factor(Compound, levels = c("C1", "DMSO")),
    LineGroup = factor(LineGroup, levels = c("C1", "DMSO"))
  )



# =========================
# Prepare replicate-level plotting data
# =========================
plot_data_reps <- plot_data %>%
  mutate(
    Time_num = as.numeric(as.character(Time))
  )

rep_points_treat <- plot_data_reps %>%
  filter(Compound %in% c("C1", "DMSO")) %>%
  mutate(LineGroup = Compound)

rep_points_baseline <- plot_data_reps %>%
  filter(Time_num == 0) %>%
  tidyr::crossing(LineGroup = c("C1", "DMSO")) %>%
  mutate(Compound = LineGroup)

plot_data_reps_lines <- bind_rows(rep_points_treat, rep_points_baseline) %>%
  mutate(
    Time = factor(Time_num, levels = sort(unique(Time_num))),
    Compound = factor(Compound, levels = c("C1", "DMSO")),
    LineGroup = factor(LineGroup, levels = c("C1", "DMSO"))
  )


# =========================
# One plot per gene_ID per pool
# =========================
dir.create("gene_pool_plots", showWarnings = FALSE)

plot_keys <- plot_summary_lines %>%
  dplyr::distinct(gene_ID, transposon_set)

for (i in seq_len(nrow(plot_keys))) {
  
  this_gene <- plot_keys$gene_ID[i]
  this_pool <- plot_keys$transposon_set[i]
  
  df_sum <- plot_summary_lines %>%
    dplyr::filter(gene_ID == this_gene, transposon_set == this_pool)
  
  df_rep <- plot_data_reps_lines %>%
    dplyr::filter(gene_ID == this_gene, transposon_set == this_pool)
  
  p <- ggplot() +
    geom_point(
      data = df_rep,
      aes(x = Time, y = normalized_reads, color = LineGroup),
      alpha = 0.35,
      size = 2,
      position = position_jitter(width = 0.05, height = 0)
    ) +
    geom_line(
      data = df_sum,
      aes(x = Time, y = mean_reads, color = LineGroup, group = LineGroup),
      linewidth = 0.8
    ) +
    geom_point(
      data = df_sum,
      aes(x = Time, y = mean_reads, color = LineGroup),
      size = 2.5
    ) +
    geom_errorbar(
      data = df_sum,
      aes(
        x = Time,
        ymin = pmax(mean_reads - se_reads, 0),
        ymax = mean_reads + se_reads,
        color = LineGroup,
        group = LineGroup
      ),
      width = 0.12,
      na.rm = TRUE
    ) +
    theme_bw() +
    labs(
      title = paste0("Gene: ", this_gene, " | Pool: ", this_pool),
      x = "Time",
      y = "Normalized reads",
      color = "Condition"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(
    filename = file.path(
      "gene_pool_plots",
      paste0("plot_", this_gene, "_pool_", this_pool, ".pdf")
    ),
    plot = p,
    width = 6,
    height = 5
  )
}


