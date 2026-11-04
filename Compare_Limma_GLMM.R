# =========================
# Compare Limma and GLMM gene lists
# =========================

base_dir <- "/home/paula-di-santo/Documents/Post_analysis_meliloti"
limma_dir <- file.path(base_dir, "Limma", "Enrichmennt")
glmm_dir <- file.path(base_dir, "Enrichmennt")
output_dir <- file.path(base_dir, "Limma_vs_GLMM")

# Read unique gene IDs from a one-ID-per-line file
read_gene_ids <- function(path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }

  ids <- trimws(readLines(path, warn = FALSE))
  sort(unique(ids[nzchar(ids)]))
}

limma_enriched <- read_gene_ids(
  file.path(limma_dir, "limma_gene_ID_enriched.txt")
)
limma_depleted <- read_gene_ids(
  file.path(limma_dir, "limma_gene_ID_depleted.txt")
)
glmm_enriched <- read_gene_ids(
  file.path(glmm_dir, "gene_ID_enriched.txt")
)
glmm_depleted <- read_gene_ids(
  file.path(glmm_dir, "gene_ID_depleted.txt")
)

# Record list membership for each gene
gene_ids <- sort(unique(c(
  limma_enriched, limma_depleted, glmm_enriched, glmm_depleted
)))

comparison <- data.frame(
  gene_ID = gene_ids,
  limma_enriched = gene_ids %in% limma_enriched,
  limma_depleted = gene_ids %in% limma_depleted,
  glmm_enriched = gene_ids %in% glmm_enriched,
  glmm_depleted = gene_ids %in% glmm_depleted
)

comparison$shared_enriched <- with(comparison, limma_enriched & glmm_enriched)
comparison$shared_depleted <- with(comparison, limma_depleted & glmm_depleted)
comparison$opposite_directions <- with(comparison,
  (limma_enriched & glmm_depleted) | (limma_depleted & glmm_enriched)
)

# Summarize overlap separately for each direction
summary_table <- data.frame(
  direction = c("Enriched", "Depleted"),
  limma_total = c(length(limma_enriched), length(limma_depleted)),
  glmm_total = c(length(glmm_enriched), length(glmm_depleted)),
  shared = c(sum(comparison$shared_enriched), sum(comparison$shared_depleted))
)

summary_table$limma_only <- summary_table$limma_total - summary_table$shared
summary_table$glmm_only <- summary_table$glmm_total - summary_table$shared

# Export the comparison and shared gene lists
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(comparison, file.path(output_dir, "gene_list_comparison.csv"),
          row.names = FALSE)
write.csv(summary_table, file.path(output_dir, "overlap_summary.csv"),
          row.names = FALSE)
writeLines(comparison$gene_ID[comparison$shared_enriched],
           file.path(output_dir, "shared_enriched.txt"))
writeLines(comparison$gene_ID[comparison$shared_depleted],
           file.path(output_dir, "shared_depleted.txt"))
writeLines(comparison$gene_ID[comparison$opposite_directions],
           file.path(output_dir, "opposite_directions.txt"))

print(summary_table, row.names = FALSE)
cat("Genes listed in opposite directions:", sum(comparison$opposite_directions), "\n")



# =========================
# Plot gene-list overlap by direction
# =========================

plot_venn <- function(title, limma_only, shared, glmm_only) {
  plot.new()
  plot.window(xlim = c(-2.2, 2.2), ylim = c(-1.3, 1.6), asp = 1)

  # Draw schematic circles with transparent fills
  symbols(
    x = c(-0.65, 0.65), y = c(0, 0),
    circles = c(1, 1), inches = FALSE, add = TRUE,
    bg = adjustcolor(c("#0072B2", "#D55E00"), alpha.f = 0.30),
    fg = c("#0072B2", "#D55E00")
  )

  text(c(-0.9, 0.9), 1.2, c("Limma", "GLMM"), font = 2)
  text(c(-1.1, 0, 1.1), 0, c(limma_only, shared, glmm_only),
       cex = 1.5, font = 2)
  title(main = title)
}

pdf(file.path(output_dir, "limma_vs_glmm_venn.pdf"), width = 9, height = 4)
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))

for (i in seq_len(nrow(summary_table))) {
  row <- summary_table[i, ]

  plot_venn(
    title = paste(row$direction, "genes"),
    limma_only = row$limma_only,
    shared = row$shared,
    glmm_only = row$glmm_only
  )
}

dev.off()

