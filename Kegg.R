library(clusterProfiler)
library(httr)

base_dir <- "/home/paula-di-santo/Documents/Post_analysis_meliloti"

# -----------------------------
# Define output directory
# -----------------------------

out_dir <- file.path(base_dir, "Glmm", "KEGG")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load gene lists
# -----------------------------

enriched_genes <- readLines(
  file.path(base_dir, "Glmm", "gene_ID_enriched.txt")
)

depleted_genes <- readLines(
  file.path(base_dir, "Glmm", "gene_ID_depleted.txt")
)

# -----------------------------
# Clean gene lists
# -----------------------------

enriched_genes <- unique(trimws(enriched_genes))
depleted_genes <- unique(trimws(depleted_genes))

enriched_genes <- enriched_genes[enriched_genes != ""]
depleted_genes <- depleted_genes[depleted_genes != ""]
depleted_genes <- depleted_genes[depleted_genes != "intergenic"]

# -----------------------------
# KEGG format (for website)
# -----------------------------

enriched_genes_kegg <- paste0("sme:", enriched_genes)
depleted_genes_kegg <- paste0("sme:", depleted_genes)

writeLines(enriched_genes_kegg, file.path(out_dir, "enriched_genes_kegg_format.txt"))
writeLines(depleted_genes_kegg, file.path(out_dir, "depleted_genes_kegg_format.txt"))

# -----------------------------
# Download KEGG gene -> KO mapping
# -----------------------------

ko_txt <- httr::content(
  httr::GET("https://rest.kegg.jp/link/ko/sme"),
  as = "text",
  encoding = "UTF-8"
)

ko_map <- read.delim(
  text = ko_txt,
  header = FALSE,
  col.names = c("gene", "ko")
)

ko_map$gene <- sub("sme:", "", ko_map$gene)
ko_map$ko   <- sub("ko:", "", ko_map$ko)

# -----------------------------
# Match genes to KO IDs
# -----------------------------

enriched_ko <- ko_map[ko_map$gene %in% enriched_genes, ]
depleted_ko <- ko_map[ko_map$gene %in% depleted_genes, ]

# -----------------------------
# Save KO mappings
# -----------------------------

write.csv(enriched_ko, file = file.path(out_dir, "enriched_gene_to_KO_mapping.csv"), row.names = FALSE)
write.csv(depleted_ko, file = file.path(out_dir, "depleted_gene_to_KO_mapping.csv"), row.names = FALSE)

# -----------------------------
# Save KEGG Reconstruct input
# -----------------------------

write.table(
  enriched_ko,
  file = file.path(out_dir, "enriched_kegg_reconstruct_input.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  depleted_ko,
  file = file.path(out_dir, "depleted_kegg_reconstruct_input.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# -----------------------------
# Optional: pathway mapping
# -----------------------------

enriched_path <- bitr_kegg(
  enriched_genes,
  fromType = "kegg",
  toType = "Path",
  organism = "sme"
)

depleted_path <- bitr_kegg(
  depleted_genes,
  fromType = "kegg",
  toType = "Path",
  organism = "sme"
)

write.csv(enriched_path, file = file.path(out_dir, "enriched_gene_to_KEGG_pathway.csv"), row.names = FALSE)
write.csv(depleted_path, file = file.path(out_dir, "depleted_gene_to_KEGG_pathway.csv"), row.names = FALSE)

# -----------------------------
# Summary
# -----------------------------

cat("Output saved in:", out_dir, "\n\n")

cat("Enriched genes:", length(enriched_genes), "\n")
cat("Mapped to KO:", length(unique(enriched_ko$gene)), "\n\n")

cat("Depleted genes:", length(depleted_genes), "\n")
cat("Mapped to KO:", length(unique(depleted_ko$gene)), "\n")