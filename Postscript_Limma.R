library(tidyverse)
library(edgeR)
library(limma)
stopifnot(requireNamespace("statmod", quietly = TRUE))


# =========================
# 1. Load input files and create folders
# =========================

# Set input and output directories
base_dir <- "/home/paula-di-santo/Documents/Post_analysis_meliloti"
limma_dir <- file.path(base_dir, "Limma")
intermediate_dir <- file.path(limma_dir, "Intermediate")

dir.create(intermediate_dir, showWarnings = FALSE, recursive = TRUE)
message("Output directory: ", intermediate_dir)

# Load sample metadata and annotation
metadata <- read_csv(
  file.path(base_dir, "MetaDataMeliloti.csv"),
  show_col_types = FALSE
)

annotation <- read_csv(
  file.path(base_dir, "Pobigaylo_y_Serratia_corrected.csv"),
  show_col_types = FALSE
)

# Set the path to the Salmon quantification files
quant_path <- "/home/paula-di-santo/Documents/Meliloti_sequences/Test_complete_sequences/quants_merged/"

# =========================
# 2. Read Salmon .sf files
# =========================
# Find all .sf files
files <- list.files(
  path = quant_path,
  pattern = "\\.sf$",
  full.names = TRUE
)

# Read and merge all files
merged_counts <- files %>%
  set_names() %>%
  map_dfr(
    ~ read_tsv(.x, show_col_types = FALSE),
    .id = "file_path"
  ) %>%
  # Extract and clean sample names
  mutate(
    sample = basename(file_path),
    sample = str_remove(sample, "\\.sf$"),
    sample = str_remove(sample, "_trimmed.*")
  ) %>%
  # Keep only relevant columns
  select(sample, Name, NumReads) %>%
  filter(
    !str_starts(sample, "unknown_"), # Filter unwanted data
    NumReads > 0
  )

# Save the result
write_csv(
  merged_counts,
  file.path(intermediate_dir, "1_merged_counts.csv") 
)

# =========================
# 2b. Add Primer_name_F from metadata
# =========================
# Get the metadata information from the NGS name
merged_parsed <- merged_counts %>%
  extract(
    col = sample,
    into = c("Internal_label_2", "Internal_label_1", "NGS_name"),
    regex = "^(\\d+_Reverse)_(\\d+_Forward)_(NG-.*)$",
    remove = FALSE
  )

final_table <- merged_parsed %>%
  left_join(
    metadata %>%
      select(NGS_name, Internal_label_1, Internal_label_2, Primer_name_F),
    by = c("NGS_name", "Internal_label_1", "Internal_label_2")
  ) %>%
  select(sample, Name, NumReads, Primer_name_F)
# Removes the rows with NA from the unknown samples
clean_table <- final_table %>%
  filter(!is.na(Primer_name_F))

write_csv(
  final_table,
  file.path(intermediate_dir, "2_merged_counts_with_primer_all.csv")
)

write_csv(
  clean_table,
  file.path(intermediate_dir, "3_merged_counts_with_primer.csv")
)

# =========================
# 2c. Annotate clean_table with gene_ID
# =========================
# After running this block the annotation of the signatures will be done, using the pool number. 

# Selecting the metadata from Pobigaylo_y_Serratia_corrected.csv
annotation_clean <- annotation %>%
  mutate(
    tag_idfentifier = as.character(tag_idfentifier),
    transposon_set = as.character(transposon_set)
  )
# Selecting metadata from merged_counts_with_primer.csv
clean_table_annotated <- clean_table %>%
  mutate(
    Name = as.character(Name),
    Primer_name_F = as.character(Primer_name_F),
    transposon_set = str_extract(Primer_name_F, "\\d+") # Identifies the number in the string POOL 
  ) %>%
  left_join( # Joining both datasets.
    annotation_clean %>%
      select(tag_idfentifier, transposon_set, gene_ID),
    by = c("Name" = "tag_idfentifier", "transposon_set" = "transposon_set")
  )

# Saving the results
write_csv(
  clean_table_annotated,
  file.path(intermediate_dir, "4_merged_counts_with_primer_and_geneID.csv")
)

# =========================
# 2d. Create pool-aware feature ID
# =========================
clean_table_final <- clean_table_annotated %>%
  filter(
    !is.na(gene_ID),
    !is.na(transposon_set),
    !grepl("intergenic", gene_ID, ignore.case = TRUE)
  ) %>%
  mutate(
    feature_ID = paste(gene_ID, transposon_set, sep = "__pool_")
  ) %>%
  select(sample, feature_ID, gene_ID, transposon_set, NumReads)


# Saving the results
write_csv(
  clean_table_final,
  file.path(intermediate_dir, "5_clean_table_final.csv")
)

# =========================
# 3. Build count matrix
# =========================
# Build a count matrix
countdata <- clean_table_final %>%
  group_by(feature_ID, sample) %>% # Aggregate counts
  summarise(NumReads = sum(NumReads), .groups = "drop") %>%
  pivot_wider( # Convert from long → wide format
    names_from = sample,
    values_from = NumReads,
    values_fill = list(NumReads = 0)
  ) %>%
  as.data.frame()

rownames(countdata) <- countdata$feature_ID # Set row names
countdata$feature_ID <- NULL # Remove redundant column
countdata[] <- lapply(countdata, as.numeric) # Ensure numeric format

write_csv(
  tibble(feature_ID = rownames(countdata), countdata),
  file.path(intermediate_dir, "6_count_matrix.csv")
)

# =========================
# 4. Parse sample names from count matrix
# =========================
# Metadata from NGS sample name
sample_info <- tibble(sample = colnames(countdata)) %>%
  extract(
    col = sample,
    into = c("Internal_label_2", "Internal_label_1", "NGS_name"),
    regex = "^(\\d+_Reverse)_(\\d+_Forward)_(.*)$",
    remove = FALSE
  )

# =========================
# 5. Prepare metadata for join 
# =========================
metadata_joinable <- metadata %>%
  mutate(
    NGS_name = as.character(NGS_name),
    Internal_label_1 = as.character(Internal_label_1),
    Internal_label_2 = as.character(Internal_label_2)
  )
# =========================
# 6. Join samples to metadata
# =========================
metadata_matched <- sample_info %>%
  inner_join(
    metadata_joinable,
    by = c(
      "NGS_name",
      "Internal_label_1",
      "Internal_label_2"
    )
  ) %>%
  filter(
    !is.na(Replicate),
    !is.na(Compound),
    !is.na(Time)
  )

# Saving the results
write_csv(
  metadata_matched,
  file.path(intermediate_dir, "7_metadata_matched.csv")
)


# =========================
# 7. Build analysis metadata
# =========================
# Select only the metadata columns you need
coldata <- metadata_matched %>%
  dplyr::select(
    sample,
    NGS_name,
    Primer_name_F,
    Lid_label,
    Replicate,
    Compound,
    Time,
    Internal_label_1,
    Internal_label_2
  ) %>%
  distinct(sample, .keep_all = TRUE) %>%
  
  # NUEVO: identificar cada cultivo seguido a través del tiempo
  mutate(
    CultureID = if_else(
      as.character(Time) == "0",
      NA_character_,
      paste(Primer_name_F, Compound, Replicate, sep = "_")
    )
  ) %>%
  
  as.data.frame()

rownames(coldata) <- coldata$sample # Move sample names into row names
coldata$sample <- NULL

# Align the count matrix with the metadata
countdata <- countdata[, rownames(coldata), drop = FALSE]

coldata <- coldata[
  match(colnames(countdata), rownames(coldata)),
  ,
  drop = FALSE
]

coldata$NGS_name <- factor(coldata$NGS_name)
coldata$Replicate <- factor(coldata$Replicate)
coldata$Compound <- factor(coldata$Compound)
coldata$Time <- factor(coldata$Time)
coldata$Primer_name_F <- factor(coldata$Primer_name_F)

# NUEVO: convertir CultureID a factor
coldata$CultureID <- factor(coldata$CultureID)

write_csv(
  tibble(sample = rownames(coldata), coldata),
  file.path(intermediate_dir, "metadata_limma.csv")
)
# =========================
# 8. Filter features
# =========================

# Retain features with at least one count across the dataset
countdata_filtered <- countdata[
  rowSums(countdata) > 0,
  ,
  drop = FALSE
]

cat("Features before filtering:", nrow(countdata), "\n")
cat("Features after filtering:", nrow(countdata_filtered), "\n")
cat("Samples in countdata_filtered:", ncol(countdata_filtered), "\n")


# =========================
# 9. Prepare count matrices and metadata for limma
# =========================

limma_input <- list()
feature_pool <- sub(".*__pool_", "", rownames(countdata_filtered))

for (pool in unique(as.character(coldata$Primer_name_F))) {

  # Keep post-treatment samples from the current pool
  metadata_pool <- droplevels(coldata[
    coldata$Primer_name_F == pool &
      coldata$Time %in% c("1", "2", "3") &
      !is.na(coldata$CultureID),
    ,
    drop = FALSE
  ])

  # Define treatment-time groups
  metadata_pool$Group <- factor(
    paste0("T", metadata_pool$Time, "_", metadata_pool$Compound),
    levels = c(
      "T1_DMSO", "T1_C1",
      "T2_DMSO", "T2_C1",
      "T3_DMSO", "T3_C1"
    )
  )

  # Use the same rounded, unnormalized counts as in the GLMM
  pool_number <- stringr::str_extract(pool, "\\d+")
  counts_pool <- round(as.matrix(countdata_filtered[
    feature_pool == pool_number,
    rownames(metadata_pool),
    drop = FALSE
  ]))

  # Remove features with zero counts in every selected sample
  counts_pool <- counts_pool[rowSums(counts_pool) > 0, , drop = FALSE]

  limma_input[[pool]] <- list(
    countdata = counts_pool,
    metadata = metadata_pool
  )
}

# =========================
# 10. Run limma-voom with culture blocking
# =========================

limma_results_by_pool <- list()
limma_checks_by_pool <- list()

for (pool in names(limma_input)) {

  message("Running limma-voom for ", pool)

  x <- limma_input[[pool]]
  metadata_pool <- x$metadata
  counts_pool <- x$countdata[, rownames(metadata_pool), drop = FALSE]
  culture <- factor(metadata_pool$CultureID)

  # Fit a separate mean for each treatment-time combination
  design <- model.matrix(~ 0 + Group, data = metadata_pool)
  colnames(design) <- sub("^Group", "", colnames(design))

  if (anyNA(culture) ||
      !identical(rownames(design), colnames(counts_pool)) ||
      qr(design)$rank < ncol(design) ||
      nrow(design) < ncol(design) + 2L) {
    stop("Invalid design or culture metadata for ", pool)
  }

  if (any(colSums(counts_pool) == 0)) {
    stop("At least one sample has zero total counts in ", pool)
  }

  # Estimate TMM normalization factors within the current pool
  dge <- edgeR::DGEList(counts = counts_pool)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")

  # Estimate correlation between measurements of the same culture
  voom_initial <- limma::voom(dge, design, plot = FALSE)

  cor_initial <- limma::duplicateCorrelation(
    voom_initial,
    design,
    block = culture
  )$consensus.correlation

  if (!is.finite(cor_initial)) {
    stop("Initial culture correlation could not be estimated for ", pool)
  }

  # Recalculate precision weights using the estimated correlation
  voom_blocked <- limma::voom(
    dge,
    design,
    block = culture,
    correlation = cor_initial,
    plot = FALSE
  )

  cor_final <- limma::duplicateCorrelation(
    voom_blocked,
    design,
    block = culture
  )$consensus.correlation

  if (!is.finite(cor_final)) {
    stop("Final culture correlation could not be estimated for ", pool)
  }

  # Define treatment effects by time and their equally weighted average
  contrasts <- limma::makeContrasts(
    T1 = T1_C1 - T1_DMSO,
    T2 = T2_C1 - T2_DMSO,
    T3 = T3_C1 - T3_DMSO,
    Average = ((T1_C1 - T1_DMSO) +
               (T2_C1 - T2_DMSO) +
               (T3_C1 - T3_DMSO)) / 3,
    levels = design
  )

  # Fit each contrast directly as a coefficient of the blocked model
  limma_results_by_pool[[pool]] <- purrr::map_dfr(
    colnames(contrasts),
    function(contrast_name) {

      contrast_design <- limma::contrastAsCoef(
        design,
        contrasts[, contrast_name, drop = FALSE]
      )

      fit <- limma::lmFit(
        voom_blocked,
        contrast_design$design,
        block = culture,
        correlation = cor_final
      )
      fit <- limma::eBayes(fit)

      coefficient <- contrast_design$coef

      tibble::tibble(
        feature_ID = rownames(counts_pool),
        pool = pool,
        contrast = contrast_name,
        log2FC = fit$coefficients[, coefficient],
        pvalue = fit$p.value[, coefficient]
      )
    }
  )

  limma_checks_by_pool[[pool]] <- tibble::tibble(
    pool = pool,
    n_features = nrow(counts_pool),
    n_samples = ncol(counts_pool),
    n_cultures = nlevels(culture),
    culture_correlation = cor_final
  )
}

# Adjust time-specific tests jointly and average-effect tests separately
limma_results <- dplyr::bind_rows(limma_results_by_pool) %>%
  tidyr::separate(
    feature_ID,
    into = c("gene_ID", "transposon_set"),
    sep = "__pool_",
    remove = FALSE
  ) %>%
  dplyr::group_by(contrast == "Average") %>%
  dplyr::mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  dplyr::ungroup()

# Classify features using the average treatment contrast
limma_average <- limma_results %>%
  dplyr::filter(contrast == "Average") %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      !is.finite(log2FC) | !is.finite(padj) ~ "Not_evaluable",
      padj < 0.05 & log2FC > 0 ~ "Enriched_in_C1",
      padj < 0.05 & log2FC < 0 ~ "Depleted_in_C1",
      TRUE ~ "Insufficient_evidence"
    )
  )

limma_checks <- dplyr::bind_rows(limma_checks_by_pool)

# Save limma results and diagnostics
readr::write_csv(
  limma_results,
  file.path(limma_dir, "limma_contrasts.csv")
)

readr::write_csv(
  limma_average,
  file.path(limma_dir, "limma_average.csv")
)

readr::write_csv(
  limma_checks,
  file.path(limma_dir, "limma_pool_checks.csv")
)

print(limma_checks, n = Inf)

limma_average %>%
  dplyr::count(direction) %>%
  print()


# =========================
# 11. Export limma candidates and DAVID gene lists
# =========================

# Split the average-contrast results by direction
limma_enriched <- limma_average %>%
  dplyr::filter(direction == "Enriched_in_C1")

limma_depleted <- limma_average %>%
  dplyr::filter(direction == "Depleted_in_C1")

# Save complete feature-level results
readr::write_csv(
  limma_enriched,
  file.path(limma_dir, "limma_enriched_in_C1.csv")
)

readr::write_csv(
  limma_depleted,
  file.path(limma_dir, "limma_depleted_in_C1.csv")
)

# Prepare unique gene IDs for DAVID
clean_gene_ids <- function(ids) {
  ids <- trimws(as.character(ids))
  sort(unique(ids[
    !is.na(ids) & nzchar(ids) &
      !grepl("intergenic", ids, ignore.case = TRUE)
  ]))
}

enriched_genes <- clean_gene_ids(limma_enriched$gene_ID)
depleted_genes <- clean_gene_ids(limma_depleted$gene_ID)

# Save gene lists in the enrichment directory
enrichment_dir <- file.path(limma_dir, "Enrichmennt")
dir.create(enrichment_dir, showWarnings = FALSE, recursive = TRUE)

writeLines(
  enriched_genes,
  file.path(enrichment_dir, "limma_gene_ID_enriched.txt")
)

writeLines(
  depleted_genes,
  file.path(enrichment_dir, "limma_gene_ID_depleted.txt")
)

# Report feature and gene counts
cat("Enriched features:", nrow(limma_enriched), "\n")
cat("Depleted features:", nrow(limma_depleted), "\n")
cat("Unique enriched genes:", length(enriched_genes), "\n")
cat("Unique depleted genes:", length(depleted_genes), "\n")
cat("Genes in both lists:",
    length(intersect(enriched_genes, depleted_genes)), "\n")