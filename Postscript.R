library(tidyverse)
library(DESeq2)
library(glmmSeq)
library(dplyr)
library(tibble)
library(tidyr)

# =========================
# 1. Load input files and create folders
# =========================
# ---- Create output directories ----

# Base directory 
base_dir <- "/home/paula-di-santo/Documents/Post_analysis_meliloti"

# Subdirectories
sub_dirs <- c("DESeq2", "Glmm", "Intermediate")

# Create base directory if it doesn't exist
if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive = TRUE)
}

# Create subdirectories
for (d in sub_dirs) {
  dir.create(file.path(base_dir, d), showWarnings = FALSE, recursive = TRUE)
}

message("Directories ready:")
print(file.path(base_dir, sub_dirs))


# ---- Load input data ----

metadata <- read_csv(
  file.path(base_dir, "MetaDataMeliloti.csv"),
  show_col_types = FALSE
)

annotation <- read_csv(
  file.path(base_dir, "Pobigaylo_y_Serratia_corrected.csv"),
  show_col_types = FALSE
)

# Quantification files path (kept separate since it's elsewhere)
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
  file.path(base_dir, "Intermediate", "1_merged_counts.csv")
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
  file.path(base_dir, "Intermediate", "2_merged_counts_with_primer_all.csv")
)

write_csv(
  clean_table,
  file.path(base_dir, "Intermediate", "3_merged_counts_with_primer.csv")
)

# =========================
# 2c. Annotate clean_table with gene_ID
# =========================
# After running this block hte annotation of the signatures will be done, using the pool number. 

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
  file.path(base_dir, "Intermediate", "4_merged_counts_with_primer_and_geneID.csv")
)

# =========================
# 2d. Create pool-aware feature ID
# =========================
clean_table_final <- clean_table_annotated %>%
  filter(!is.na(gene_ID), !is.na(transposon_set)) %>%
  mutate(
    feature_ID = paste(gene_ID, transposon_set, sep = "__pool_")
  ) %>%
  select(sample, feature_ID, gene_ID, transposon_set, NumReads)


# Saving the results
write_csv(
  clean_table_final,
  file.path(base_dir, "Intermediate", "5_clean_table_final.csv")
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
  file.path(base_dir, "Intermediate", "6_count_matrix.csv")
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
  file.path(base_dir, "Intermediate", "7_metadata_matched.csv")
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
  as.data.frame()

rownames(coldata) <- coldata$sample # Move sample names into row names
coldata$sample <- NULL
# Align the count matrix with the metadata
# Reorders and subsets the columns of countdata so they match the sample names in coldata
# Only samples found in coldata are kept
# The column order is forced to match rownames(coldata)
countdata <- countdata[, rownames(coldata), drop = FALSE] 
# This then reorders the rows of coldata so they match the columns of countdata.
coldata <- coldata[match(colnames(countdata), rownames(coldata)), , drop = FALSE]

coldata$NGS_name  <- factor(coldata$NGS_name)
coldata$Replicate <- factor(coldata$Replicate)
coldata$Compound  <- factor(coldata$Compound)
coldata$Time      <- factor(coldata$Time)
coldata$Primer_name_F  <- factor(coldata$Primer_name_F)


write_csv(
  tibble(sample = rownames(coldata), coldata),
  file.path(base_dir, "DESeq2", "metadata_deseq2.csv")
)
# =========================
# 8. Filter features
# =========================
countdata <- countdata[rowSums(countdata) > 0, , drop = FALSE]

keep_features <- rowSums(countdata >= 10) >= 3 # filtering for a minimun of reads 
countdata_filtered <- countdata[keep_features, , drop = FALSE]

cat("Features before filtering:", nrow(countdata), "\n")
cat("Features after filtering:", nrow(countdata_filtered), "\n")
cat("Samples in countdata_filtered:", ncol(countdata_filtered), "\n")


# =========================
# 9. DESeq2 normalization and dispersion
# =========================
# Each of these columns is converted to a factor (categorical variables)
coldata$Time <- droplevels(factor(coldata$Time))
coldata$Compound <- droplevels(factor(coldata$Compound))
coldata$Replicate <- droplevels(factor(coldata$Replicate))
coldata$NGS_name <- droplevels(factor(coldata$NGS_name))
coldata$Primer_name_F <- droplevels(factor(coldata$Primer_name_F))



dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(countdata_filtered)), # This converts your filtered count table into a matrix and rounds values to integers
  colData   = coldata,
  design    = ~ 1 # no experimental variable is included yet, only an intercept
)
# Size factors are used to normalize for differences in sequencing depth between samples.
dds <- estimateSizeFactors(dds, type = "poscounts") # type = "poscounts" This method is useful when your dataset has many zeros
# Dispersion measures how much counts vary across samples beyond what is expected from simple Poisson noise.
dds <- estimateDispersions(dds)
# This extracts the estimated size factor for each sample.
size_factors <- sizeFactors(dds)
# This extracts one dispersion estimate per feature. So each gene/feature gets a dispersion value
dispersion_values <- dispersions(dds)
# This returns the counts after normalization by size factors.
normalized_counts <- counts(dds, normalized = TRUE)

write.csv(
  normalized_counts,
  file = file.path(base_dir, "DESeq2", "normalized_counts_DESeq2.csv")
)


# =========================
# 10. Prepare data for glmmSeq
# =========================
# Keep only features with valid dispersion estimates. Because glmmSeq needs valid dispersions.
dispersion_keep <- !is.na(dispersion_values)
# This takes the raw counts from dds and keeps only the rows where dispersion is available.
countdata_glmm <- counts(dds, normalized = FALSE)[dispersion_keep, , drop = FALSE]

# subsets dispersions to the same kept genes
dispersion_values_glmm <- dispersion_values[dispersion_keep]
# ensures they are numeric
dispersion_values_glmm <- as.numeric(dispersion_values_glmm)
# names them with the corresponding feature IDs. # This step is important so the dispersions line up with the rows of countdata_glmm.
names(dispersion_values_glmm) <- rownames(countdata_glmm)
# Build metadata for glmmSeq
coldata_glmm <- coldata[, c("Time", "Compound", "NGS_name", "Primer_name_F"), drop = FALSE]
# Create a combined experimental group. Instead of modeling Time and Compound separately, you are treating each combination as its own group.
coldata_glmm$Group <- factor(
  paste0("T", as.character(coldata_glmm$Time), "_", as.character(coldata_glmm$Compound)),
  levels = c(
    "T0_None",
    "T1_C1", "T1_DMSO",
    "T2_C1", "T2_DMSO",
    "T3_C1", "T3_DMSO"
  )
)
# Remove samples with undefined groups
keep_samples <- !is.na(coldata_glmm$Group)
# keeps only metadata rows for valid samples
coldata_glmm <- droplevels(coldata_glmm[keep_samples, , drop = FALSE])
# keeps the same samples in the count matrix
countdata_glmm <- countdata_glmm[, rownames(coldata_glmm), drop = FALSE]
# keeps the matching size factors for those same samples
size_factors_glmm <- size_factors[rownames(coldata_glmm)]

if (!all(colnames(countdata_glmm) == rownames(coldata_glmm))) {
  stop("Sample names in countdata_glmm and coldata_glmm do not match.")
}

write.csv(
  countdata_glmm,
  file = file.path(base_dir, "Glmm", "countdata_glmm.csv")
)

# =========================
# 11. Run glmmSeq
# =========================

results_glmm <- glmmSeq(
  ~ Group + (1 | Primer_name_F),
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

write.csv(
  glmm_results_df,
  file = file.path(base_dir, "Glmm", "glmmSeq_results.csv")
)


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


# =========================
# 14. Save diff results
# =========================

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


#results_final <- results_final %>%
 #mutate(
 #  score = abs(diff_T1) * (pval_T1 < 0.05) +
 #          abs(diff_T2) * (pval_T2 < 0.05) +
 #          abs(diff_T3) * (pval_T3 < 0.05)
 # )



results_final <- results_final %>%
 mutate(
   score_directional =
     diff_T1*(pval_T1<0.05) +
     diff_T2*(pval_T2<0.05) +
     diff_T3*(pval_T3<0.05)
 )


results_final <- results_final %>%
 mutate(
   direction = case_when(
     score_directional > 0.5 ~ "Enriched_in_C1",
     score_directional < -0.5 ~ "Depleted_in_C1",
     TRUE ~ "Neutral"
   )
 )

results_enriched <- results_final %>%
  filter(direction == "Enriched_in_C1")

results_depleted <- results_final %>%
  filter(direction == "Depleted_in_C1")

results_neutral <- results_final %>%
  filter(direction == "Neutral")

  write.csv(
  results_enriched,
  file = file.path(base_dir, "Glmm", "GLMM_results_enriched.csv"),
  row.names = FALSE
)

write.csv(
  results_depleted,
  file = file.path(base_dir, "Glmm", "GLMM_results_depleted.csv"),
  row.names = FALSE
)

