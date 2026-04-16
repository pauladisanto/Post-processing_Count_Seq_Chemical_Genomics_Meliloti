library(tidyverse)
library(DESeq2)
library(glmmSeq)
library(dplyr)
library(tibble)
library(tidyr)

# =========================
# 1. Load input files
# =========================
metadata <- read_csv(
  "/home/paula-di-santo/Documents/Post_analysis_meliloti/MetaDataMeliloti.csv",
  show_col_types = FALSE
)

annotation <- read_csv(
  "/home/paula-di-santo/Documents/Post_analysis_meliloti/Pobigaylo_y_Serratia_corrected.csv",
  show_col_types = FALSE
)

quant_path <- "/home/paula-di-santo/Documents/Meliloti_sequences/Test_complete_sequences/quants_merged/"

# =========================
# 2. Read Salmon .sf files
# =========================
files <- list.files(
  path = quant_path,
  pattern = "\\.sf$",
  full.names = TRUE
)

merged_counts <- files %>%
  set_names() %>%
  map_dfr(
    ~ read_tsv(.x, show_col_types = FALSE),
    .id = "file_path"
  ) %>%
  mutate(
    sample = basename(file_path),
    sample = str_remove(sample, "\\.sf$"),
    sample = str_remove(sample, "_trimmed.*")
  ) %>%
  select(sample, Name, NumReads) %>%
  filter(
    !str_starts(sample, "unknown_"),
    NumReads > 0
  )

# Optional
# write_csv(merged_counts, "merged_counts.csv")

# =========================
# 2b. Add Primer_name_F from metadata
# =========================
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

clean_table <- final_table %>%
  filter(!is.na(Primer_name_F))

# Optional
# write_csv(final_table, "merged_counts_with_primer_all.csv")
# write_csv(clean_table, "merged_counts_with_primer.csv")

# =========================
# 2c. Annotate clean_table with gene_ID
# =========================
annotation_clean <- annotation %>%
  mutate(
    tag_idfentifier = as.character(tag_idfentifier),
    transposon_set = as.character(transposon_set)
  )

clean_table_annotated <- clean_table %>%
  mutate(
    Name = as.character(Name),
    Primer_name_F = as.character(Primer_name_F),
    transposon_set = str_extract(Primer_name_F, "\\d+")
  ) %>%
  left_join(
    annotation_clean %>%
      select(tag_idfentifier, transposon_set, gene_ID),
    by = c("Name" = "tag_idfentifier", "transposon_set" = "transposon_set")
  )

write_csv(clean_table_annotated, "merged_counts_with_primer_and_geneID.csv")

# =========================
# 2d. Create pool-aware feature ID
# =========================
clean_table_final <- clean_table_annotated %>%
  filter(!is.na(gene_ID), !is.na(transposon_set)) %>%
  mutate(
    feature_ID = paste(gene_ID, transposon_set, sep = "__pool_")
  ) %>%
  select(sample, feature_ID, gene_ID, transposon_set, NumReads)

# =========================
# 3. Build count matrix
# =========================
countdata <- clean_table_final %>%
  group_by(feature_ID, sample) %>%
  summarise(NumReads = sum(NumReads), .groups = "drop") %>%
  pivot_wider(
    names_from = sample,
    values_from = NumReads,
    values_fill = list(NumReads = 0)
  ) %>%
  as.data.frame()

rownames(countdata) <- countdata$feature_ID
countdata$feature_ID <- NULL
countdata[] <- lapply(countdata, as.numeric)

write_csv(
  tibble(feature_ID = rownames(countdata), countdata),
  "count_matrix.csv"
)

# =========================
# 4. Parse sample names from count matrix
# =========================
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
    Internal_label_1_chr = as.character(Internal_label_1),
    Internal_label_2_chr = as.character(Internal_label_2)
  )

# =========================
# 6. Join samples to metadata
# =========================
metadata_matched <- sample_info %>%
  inner_join(
    metadata_joinable,
    by = c(
      "NGS_name" = "NGS_name",
      "Internal_label_2" = "Internal_label_2_chr",
      "Internal_label_1" = "Internal_label_1_chr"
    )
  ) %>%
  filter(
    !is.na(Replicate),
    !is.na(Compound),
    !is.na(Time)
  )

# =========================
# 7. Build analysis metadata
# =========================
label1_col <- intersect(
  c("Internal_label_1.x", "Internal_label_1", "Internal_label_1.y", "Internal_label_1_chr"),
  colnames(metadata_matched)
)[1]

label2_col <- intersect(
  c("Internal_label_2.x", "Internal_label_2", "Internal_label_2.y", "Internal_label_2_chr"),
  colnames(metadata_matched)
)[1]

if (is.na(label1_col) || is.na(label2_col)) {
  stop("Could not find Internal_label_1 / Internal_label_2 columns in metadata_matched")
}

coldata <- metadata_matched %>%
  dplyr::select(
    sample,
    NGS_name,
    Lid_label,
    Replicate,
    Compound,
    Time,
    all_of(label1_col),
    all_of(label2_col)
  ) %>%
  distinct(sample, .keep_all = TRUE) %>%
  rename(
    Internal_label_1 = all_of(label1_col),
    Internal_label_2 = all_of(label2_col)
  ) %>%
  as.data.frame()

rownames(coldata) <- coldata$sample
coldata$sample <- NULL

countdata <- countdata[, rownames(coldata), drop = FALSE]
coldata <- coldata[match(colnames(countdata), rownames(coldata)), , drop = FALSE]

coldata$NGS_name  <- factor(coldata$NGS_name)
coldata$Replicate <- factor(coldata$Replicate)
coldata$Compound  <- factor(coldata$Compound)
coldata$Time      <- factor(coldata$Time)

write_csv(
  tibble(sample = rownames(coldata), coldata),
  "metadata_deseq2.csv"
)

# =========================
# 8. Filter features
# =========================
countdata <- countdata[rowSums(countdata) > 0, , drop = FALSE]

keep_features <- rowSums(countdata >= 10) >= 3
countdata_filtered <- countdata[keep_features, , drop = FALSE]

cat("Features before filtering:", nrow(countdata), "\n")
cat("Features after filtering:", nrow(countdata_filtered), "\n")
cat("Samples in countdata_filtered:", ncol(countdata_filtered), "\n")

if (nrow(countdata_filtered) == 0) {
  stop("No features left after filtering. Relax the filtering threshold.")
}

if (ncol(countdata_filtered) == 0) {
  stop("No samples left in countdata_filtered.")
}

if (!all(colnames(countdata_filtered) == rownames(coldata))) {
  stop("Sample names in countdata_filtered and coldata do not match.")
}

# =========================
# 9. DESeq2 normalization and dispersion
# =========================
coldata$Time <- droplevels(factor(coldata$Time))
coldata$Compound <- droplevels(factor(coldata$Compound))
coldata$Replicate <- droplevels(factor(coldata$Replicate))
coldata$NGS_name <- droplevels(factor(coldata$NGS_name))

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(countdata_filtered)),
  colData   = coldata,
  design    = ~ 1
)

dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- estimateDispersions(dds)

size_factors <- sizeFactors(dds)
dispersion_values <- dispersions(dds)
normalized_counts <- counts(dds, normalized = TRUE)

write.csv(normalized_counts, file = "normalized_counts_DESeq2.csv")

# =========================
# 10. Prepare data for glmmSeq
# =========================
dispersion_keep <- !is.na(dispersion_values)

countdata_glmm <- counts(dds, normalized = FALSE)[dispersion_keep, , drop = FALSE]
dispersion_values_glmm <- dispersion_values[dispersion_keep]
dispersion_values_glmm <- as.numeric(dispersion_values_glmm)
names(dispersion_values_glmm) <- rownames(countdata_glmm)

coldata_glmm <- coldata[, c("Time", "Compound", "NGS_name"), drop = FALSE]

coldata_glmm$Group <- factor(
  paste0("T", as.character(coldata_glmm$Time), "_", as.character(coldata_glmm$Compound)),
  levels = c(
    "T0_None",
    "T1_C1", "T1_DMSO",
    "T2_C1", "T2_DMSO",
    "T3_C1", "T3_DMSO"
  )
)

keep_samples <- !is.na(coldata_glmm$Group)

coldata_glmm <- droplevels(coldata_glmm[keep_samples, , drop = FALSE])
countdata_glmm <- countdata_glmm[, rownames(coldata_glmm), drop = FALSE]
size_factors_glmm <- size_factors[rownames(coldata_glmm)]

if (!all(colnames(countdata_glmm) == rownames(coldata_glmm))) {
  stop("Sample names in countdata_glmm and coldata_glmm do not match.")
}

# =========================
# 11. Run glmmSeq
# =========================
results_glmm <- glmmSeq(
  ~ Group + (1 | NGS_name),
  countdata   = as.data.frame(countdata_glmm),
  metadata    = coldata_glmm,
  sizeFactors = size_factors_glmm,
  dispersion  = dispersion_values_glmm,
  progress    = TRUE
)

# =========================
# 12. Save glmmSeq results
# =========================
coef_df   <- as.data.frame(results_glmm@stats$coef)   %>% tibble::rownames_to_column("feature_ID")
stderr_df <- as.data.frame(results_glmm@stats$stdErr) %>% tibble::rownames_to_column("feature_ID")
chisq_df  <- as.data.frame(results_glmm@stats$Chisq)  %>% tibble::rownames_to_column("feature_ID")
pval_df   <- as.data.frame(results_glmm@stats$pvals)  %>% tibble::rownames_to_column("feature_ID")
res_df    <- as.data.frame(results_glmm@stats$res)    %>% tibble::rownames_to_column("feature_ID")

colnames(coef_df)[-1]   <- paste0("coef_", colnames(coef_df)[-1])
colnames(stderr_df)[-1] <- paste0("stderr_", colnames(stderr_df)[-1])
colnames(chisq_df)[-1]  <- paste0("chisq_", colnames(chisq_df)[-1])
colnames(pval_df)[-1]   <- paste0("pval_", colnames(pval_df)[-1])
colnames(res_df)[-1]    <- paste0("stat_", colnames(res_df)[-1])

glmm_results_df <- res_df %>%
  left_join(coef_df, by = "feature_ID") %>%
  left_join(stderr_df, by = "feature_ID") %>%
  left_join(chisq_df, by = "feature_ID") %>%
  left_join(pval_df, by = "feature_ID") %>%
  tidyr::separate(
    feature_ID,
    into = c("gene_ID", "transposon_set"),
    sep = "__pool_",
    remove = FALSE
  )

write_csv(glmm_results_df, "glmmSeq_results.csv")

# =========================
# 13. Save objects needed for plotting script
# =========================
write_csv(clean_table_annotated, "clean_table_annotated.csv")

write_csv(
  tibble(sample = names(size_factors), size_factor = size_factors),
  "size_factors.csv"
)

write_csv(
  tibble(sample = rownames(coldata), coldata),
  "metadata_deseq2.csv"
)



results_final <- glmm_results_df %>%
  mutate(
    # Differences C1 vs DMSO
    diff_T1 = coef_GroupT1_C1 - coef_GroupT1_DMSO,
    diff_T2 = coef_GroupT2_C1 - coef_GroupT2_DMSO,
    diff_T3 = coef_GroupT3_C1 - coef_GroupT3_DMSO,

    # Standard errors
    se_T1 = sqrt(stderr_GroupT1_C1^2 + stderr_GroupT1_DMSO^2),
    se_T2 = sqrt(stderr_GroupT2_C1^2 + stderr_GroupT2_DMSO^2),
    se_T3 = sqrt(stderr_GroupT3_C1^2 + stderr_GroupT3_DMSO^2),

    # Z-scores
    z_T1 = diff_T1 / se_T1,
    z_T2 = diff_T2 / se_T2,
    z_T3 = diff_T3 / se_T3,

    # P-values
    pval_T1 = 2 * pnorm(-abs(z_T1)),
    pval_T2 = 2 * pnorm(-abs(z_T2)),
    pval_T3 = 2 * pnorm(-abs(z_T3))
  )


results_final <- results_final %>%
 mutate(
   score = abs(diff_T1) * (pval_T1 < 0.05) +
           abs(diff_T2) * (pval_T2 < 0.05) +
           abs(diff_T3) * (pval_T3 < 0.05)
  )


  top_mutants <- results_final %>%
  filter(score > 0) %>%   # keep only interesting ones
  arrange(desc(score))


  write_csv(top_mutants, "top_mutants_ranked.csv")

top_C1 <- top_mutants %>%
  filter(diff_T1 > 0 | diff_T2 > 0 | diff_T3 > 0)

write_csv(top_C1, "top_C1_enriched.csv")


top_sensitive <- top_mutants %>%
  filter(diff_T1 < 0 | diff_T2 < 0 | diff_T3 < 0)

write_csv(top_sensitive, "top_C1_sensitive.csv")