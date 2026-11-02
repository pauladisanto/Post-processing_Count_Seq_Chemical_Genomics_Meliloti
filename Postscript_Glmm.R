library(tidyverse)
library(DESeq2)
library(glmmSeq)
library(dplyr)
library(tibble)
library(tidyr)
library(glmmTMB)

# =========================
# 1. Load input files and create folders
# =========================
# ---- Create output directories ----

# Base directory 
base_dir <- "/home/paula-di-santo/Documents/Post_analysis_meliloti"

# Subdirectories
sub_dirs <- c("DESeq2", "Glmm", "Intermediate", "Plotting")

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
  file.path(base_dir, "DESeq2", "metadata_deseq2.csv")
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
# 9. Normalize counts within each pool (this is for plotting))
# =========================

# Store one DESeq2 object per pool
dds_by_pool <- list()

# Extract the pool number from each feature ID
feature_pool <- sub(
  ".*__pool_",
  "",
  rownames(countdata_filtered)
)

for (pool in unique(as.character(coldata$Primer_name_F))) {

  # Select samples belonging to the current pool
  metadata_pool <- droplevels(
    coldata[coldata$Primer_name_F == pool, , drop = FALSE]
  )

  # Select features and samples belonging to the current pool
  pool_number <- str_extract(pool, "\\d+")

  counts_pool <- round(as.matrix(
    countdata_filtered[
      feature_pool == pool_number,
      rownames(metadata_pool),
      drop = FALSE
    ]
  ))

  # Remove features with zero counts in all samples after rounding
  counts_pool <- counts_pool[
    rowSums(counts_pool) > 0,
    ,
    drop = FALSE
  ]

  # Create a DESeq2 object for normalization
  dds_pool <- DESeqDataSetFromMatrix(
    countData = counts_pool,
    colData = metadata_pool,
    design = ~ 1
  )

  # Estimate size factors within the current pool
  dds_pool <- estimateSizeFactors(dds_pool, type = "poscounts")

  # Store the object for subsequent analysis
  dds_by_pool[[pool]] <- dds_pool

  # Save normalized counts
  write.csv(
    counts(dds_pool, normalized = TRUE),
    file = file.path(
      base_dir,
      "DESeq2",
      paste0("normalized_counts_", pool, ".csv")
    )
  )
}


# =========================
# 10. Prepare data for glmmTMB within each pool
# =========================

# Store model inputs for each pool
glmm_input <- list()

for (pool in names(dds_by_pool)) {

  dds_pool <- dds_by_pool[[pool]]

  # Extract sample metadata
  metadata_pool <- as.data.frame(colData(dds_pool))

  # Select post-treatment samples with a defined culture ID
  metadata_pool <- metadata_pool[
    metadata_pool$Time %in% c("1", "2", "3") &
      !is.na(metadata_pool$CultureID),
    ,
    drop = FALSE
  ]

  metadata_pool <- droplevels(metadata_pool)

  # Define treatment-time groups
  metadata_pool$Group <- factor(
    paste0("T", metadata_pool$Time, "_", metadata_pool$Compound),
    levels = c(
      "T1_DMSO", "T1_C1",
      "T2_DMSO", "T2_C1",
      "T3_DMSO", "T3_C1"
    )
  )

  # Extract raw counts in the same order as the metadata
  counts_pool <- counts(dds_pool, normalized = FALSE)[
    ,
    rownames(metadata_pool),
    drop = FALSE
  ]

  # Remove features with zero counts in all selected samples
  counts_pool <- counts_pool[
    rowSums(counts_pool) > 0,
    ,
    drop = FALSE
  ]

  # Extract the corresponding size factors
  size_factors_pool <- sizeFactors(dds_pool)[
    rownames(metadata_pool)
  ]

  # Store inputs for model fitting
  glmm_input[[pool]] <- list(
    countdata = counts_pool,
    metadata = metadata_pool,
    sizeFactors = size_factors_pool
  )
}



# =========================
# 11. Fit mixed models within each pool
# =========================

# Store fitted models by pool and feature
results_glmm_by_pool <- list()

for (pool in names(glmm_input)) {

  message("Fitting models for ", pool)

  x <- glmm_input[[pool]]

  # Prepare sample metadata and normalization offsets
  model_data <- x$metadata
  model_data$log_size_factor <- log(x$sizeFactors)

  pool_models <- list()

  for (feature in rownames(x$countdata)) {

    # Add raw counts for the current feature
    model_data$count <- as.numeric(x$countdata[feature, ])

    # Fit the model and retain any fitting errors
    pool_models[[feature]] <- tryCatch(
      glmmTMB(
        count ~ Group + (1 | CultureID) + offset(log_size_factor),
        data = model_data,
        family = nbinom2()
      ),
      error = function(e) e
    )
  }

  results_glmm_by_pool[[pool]] <- pool_models
}

# =========================
# 12. Save models and check numerical validity
# =========================

# Save complete models for subsequent contrasts
saveRDS(
  results_glmm_by_pool,
  file.path(base_dir, "Glmm", "fitted_models_by_pool.rds")
)

# Collect fitting errors and numerical diagnostics
model_checks <- purrr::imap_dfr(
  results_glmm_by_pool,
  function(pool_models, pool) {

    purrr::imap_dfr(pool_models, function(fit, feature) {

      # Record models that failed during fitting
      if (inherits(fit, "error")) {
        return(tibble::tibble(
          feature_ID = feature,
          pool = pool,
          convergence = NA_integer_,
          positive_hessian = NA,
          finite_fixed_effects = NA,
          finite_covariance = NA,
          status = "fit_error",
          message = conditionMessage(fit)
        ))
      }

      # Check optimizer convergence and Hessian
      converged <- isTRUE(fit$fit$convergence == 0)
      positive_hessian <- isTRUE(fit$sdr$pdHess)

      # Check fixed-effect estimates and their covariance matrix
      finite_fixed_effects <- all(
        is.finite(fixef(fit)$cond)
      )

      finite_covariance <- all(
        is.finite(as.matrix(vcov(fit)$cond))
      )

      tibble::tibble(
        feature_ID = feature,
        pool = pool,
        convergence = fit$fit$convergence,
        positive_hessian = positive_hessian,
        finite_fixed_effects = finite_fixed_effects,
        finite_covariance = finite_covariance,
        status = dplyr::if_else(
          converged & positive_hessian &
            finite_fixed_effects & finite_covariance,
          "passed_initial_checks",
          "needs_review"
        ),
        message = paste(fit$fit$message, collapse = "; ")
      )
    })
  }
)

# Recover gene and pool identifiers from feature IDs
model_checks <- model_checks %>%
  tidyr::separate(
    feature_ID,
    into = c("gene_ID", "transposon_set"),
    sep = "__pool_",
    remove = FALSE
  )

# Save model diagnostics
write_csv(
  model_checks,
  file.path(base_dir, "Glmm", "model_checks.csv")
)

# Summarize model status
model_checks %>%
  dplyr::count(status) %>%
  print()

# Select features that passed all initial numerical checks
features_for_contrasts <- model_checks %>%
  dplyr::filter(status == "passed_initial_checks") %>%
  dplyr::select(pool, feature_ID)

# =========================
# 12b. Retry problematic models using BFGS
# =========================

# Select models requiring review
models_to_retry <- model_checks %>%
  dplyr::filter(status == "needs_review")

# Preserve the original fits for comparison
original_problem_models <- list()
recovered <- 0L

for (i in seq_len(nrow(models_to_retry))) {

  pool <- models_to_retry$pool[i]
  feature <- models_to_retry$feature_ID[i]

  original_problem_models[[feature]] <-
    results_glmm_by_pool[[pool]][[feature]]

  # Reconstruct the data for this specific feature
  x <- glmm_input[[pool]]
  model_data <- x$metadata
  model_data$log_size_factor <- log(x$sizeFactors)
  model_data$count <- as.numeric(x$countdata[feature, ])

  # Refit the same model with a different optimizer
  retry_fit <- tryCatch(
    glmmTMB(
      count ~ Group + (1 | CultureID) + offset(log_size_factor),
      data = model_data,
      family = nbinom2(),
      control = glmmTMBControl(
        optimizer = optim,
        optArgs = list(method = "BFGS")
      )
    ),
    error = function(e) e
  )

  # Replace the original fit only if all initial numerical checks pass
  if (inherits(retry_fit, "glmmTMB") &&
      isTRUE(retry_fit$fit$convergence == 0) &&
      isTRUE(retry_fit$sdr$pdHess) &&
      all(is.finite(fixef(retry_fit)$cond)) &&
      all(is.finite(as.matrix(vcov(retry_fit)$cond)))) {

    results_glmm_by_pool[[pool]][[feature]] <- retry_fit
    recovered <- recovered + 1L
  }
}

message("Models recovered with BFGS: ", recovered)


# IMPORTANT: En una ejecución completa, el orden debe ser: 11 → 12 → 12b → 12 nuevamente. 
# Así, la selección final incluye los modelos recuperados con BFGS.
# =========================
# 12. Save models and check numerical validity
# =========================

# Save complete models for subsequent contrasts
saveRDS(
  results_glmm_by_pool,
  file.path(base_dir, "Glmm", "fitted_models_by_pool.rds")
)

# Collect fitting errors and numerical diagnostics
model_checks <- purrr::imap_dfr(
  results_glmm_by_pool,
  function(pool_models, pool) {

    purrr::imap_dfr(pool_models, function(fit, feature) {

      # Record models that failed during fitting
      if (inherits(fit, "error")) {
        return(tibble::tibble(
          feature_ID = feature,
          pool = pool,
          convergence = NA_integer_,
          positive_hessian = NA,
          finite_fixed_effects = NA,
          finite_covariance = NA,
          status = "fit_error",
          message = conditionMessage(fit)
        ))
      }

      # Check optimizer convergence and Hessian
      converged <- isTRUE(fit$fit$convergence == 0)
      positive_hessian <- isTRUE(fit$sdr$pdHess)

      # Check fixed-effect estimates and their covariance matrix
      finite_fixed_effects <- all(
        is.finite(fixef(fit)$cond)
      )

      finite_covariance <- all(
        is.finite(as.matrix(vcov(fit)$cond))
      )

      tibble::tibble(
        feature_ID = feature,
        pool = pool,
        convergence = fit$fit$convergence,
        positive_hessian = positive_hessian,
        finite_fixed_effects = finite_fixed_effects,
        finite_covariance = finite_covariance,
        status = dplyr::if_else(
          converged & positive_hessian &
            finite_fixed_effects & finite_covariance,
          "passed_initial_checks",
          "needs_review"
        ),
        message = paste(fit$fit$message, collapse = "; ")
      )
    })
  }
)

# Recover gene and pool identifiers from feature IDs
model_checks <- model_checks %>%
  tidyr::separate(
    feature_ID,
    into = c("gene_ID", "transposon_set"),
    sep = "__pool_",
    remove = FALSE
  )

# Save model diagnostics
write_csv(
  model_checks,
  file.path(base_dir, "Glmm", "model_checks.csv")
)

# Summarize model status
model_checks %>%
  dplyr::count(status) %>%
  print()

# Select features that passed all initial numerical checks
features_for_contrasts <- model_checks %>%
  dplyr::filter(status == "passed_initial_checks") %>%
  dplyr::select(pool, feature_ID)


# =========================
# 13. Calculate treatment contrasts within each time point
# =========================

contrast_results <- vector("list", nrow(features_for_contrasts))

for (i in seq_len(nrow(features_for_contrasts))) {

  pool <- features_for_contrasts$pool[i]
  feature <- features_for_contrasts$feature_ID[i]

  fit <- results_glmm_by_pool[[pool]][[feature]]

  # Extract fixed effects and their full covariance matrix
  beta <- fixef(fit)$cond
  V <- as.matrix(vcov(fit)$cond)
  V <- V[names(beta), names(beta), drop = FALSE]

  # Define C1 versus DMSO contrasts on the natural log scale
  # T1_DMSO is the reference group
  L <- matrix(
    0,
    nrow = 3,
    ncol = length(beta),
    dimnames = list(c("T1", "T2", "T3"), names(beta))
  )

  L["T1", "GroupT1_C1"] <- 1

  L["T2", "GroupT2_C1"] <- 1
  L["T2", "GroupT2_DMSO"] <- -1

  L["T3", "GroupT3_C1"] <- 1
  L["T3", "GroupT3_DMSO"] <- -1

  # Calculate effects and standard errors using full covariance
  effect <- as.numeric(L %*% beta)
  se <- sqrt(diag(L %*% V %*% t(L)))

  # Identify contrasts with finite estimates and positive standard errors
  valid <- is.finite(effect) & is.finite(se) & se > 0

  z <- rep(NA_real_, 3)
  z[valid] <- effect[valid] / se[valid]

  # Store one row per feature and time point
  contrast_results[[i]] <- tibble(
    feature_ID = feature,
    pool = pool,
    Time = rownames(L),
    log_effect = effect,
    SE = as.numeric(se),
    log2FC = effect / log(2),
    log2FC_lower = (effect - qnorm(0.975) * se) / log(2),
    log2FC_upper = (effect + qnorm(0.975) * se) / log(2),
    z = z,
    pvalue = 2 * pnorm(-abs(z)),
    contrast_valid = valid
  )
}

# Combine results across all pools
results_final <- bind_rows(contrast_results) %>%
  tidyr::separate(
    feature_ID,
    into = c("gene_ID", "transposon_set"),
    sep = "__pool_",
    remove = FALSE
  ) %>%
  mutate(
    # Adjust jointly across features and all three time points
    padj = p.adjust(pvalue, method = "BH")
  )

# Save treatment contrasts
write_csv(
  results_final,
  file.path(base_dir, "Glmm", "treatment_contrasts.csv")
)


# =========================
# 14. Combine time-point results for each feature
# =========================

results_trajectory <- results_final %>%
  dplyr::select(
    feature_ID,
    gene_ID,
    transposon_set,
    Time,
    log2FC,
    log2FC_lower,
    log2FC_upper,
    padj,
    contrast_valid
  ) %>%
  tidyr::pivot_wider(
    names_from = Time,
    values_from = c(
      log2FC,
      log2FC_lower,
      log2FC_upper,
      padj,
      contrast_valid
    ),
    names_glue = "{.value}_{Time}"
  )

# Save one row per feature with all time points
write_csv(
  results_trajectory,
  file.path(base_dir, "Glmm", "mutant_trajectories.csv")
)



# =========================
# 14b. Identify candidates using the average treatment effect
# =========================

average_results <- vector("list", nrow(features_for_contrasts))

for (i in seq_len(nrow(features_for_contrasts))) {

  pool <- features_for_contrasts$pool[i]
  feature <- features_for_contrasts$feature_ID[i]
  fit <- results_glmm_by_pool[[pool]][[feature]]

  # Extract fixed effects and their covariance matrix
  beta <- fixef(fit)$cond
  V <- as.matrix(vcov(fit)$cond)
  V <- V[names(beta), names(beta), drop = FALSE]

  # Define the average C1 versus DMSO contrast across T1-T3
  # T1_DMSO is the reference group
  weights <- setNames(rep(0, length(beta)), names(beta))

  weights[c("GroupT1_C1", "GroupT2_C1", "GroupT3_C1")] <- 1 / 3
  weights[c("GroupT2_DMSO", "GroupT3_DMSO")] <- -1 / 3

  # Calculate the average effect and its variance
  effect <- sum(weights * beta)
  variance <- as.numeric(t(weights) %*% V %*% weights)

  valid <- is.finite(effect) &&
    is.finite(variance) &&
    variance > 0

  se <- if (valid) sqrt(variance) else NA_real_
  pvalue <- if (valid) {
    2 * pnorm(-abs(effect / se))
  } else {
    NA_real_
  }

  average_results[[i]] <- tibble(
    feature_ID = feature,
    pool = pool,
    mean_log2FC = effect / log(2),
    SE_log2 = se / log(2),
    CI_lower = (effect - qnorm(0.975) * se) / log(2),
    CI_upper = (effect + qnorm(0.975) * se) / log(2),
    pvalue = pvalue,
    contrast_valid = valid
  )
}

# Adjust across all evaluated feature-level average contrasts
candidate_results <- bind_rows(average_results) %>%
  tidyr::separate(
    feature_ID,
    into = c("gene_ID", "transposon_set"),
    sep = "__pool_",
    remove = FALSE
  ) %>%
  mutate(
    padj = p.adjust(pvalue, method = "BH"),
    direction = case_when(
      !contrast_valid ~ "Not_evaluable",
      padj < 0.05 & mean_log2FC > 0 ~ "Enriched_in_C1",
      padj < 0.05 & mean_log2FC < 0 ~ "Depleted_in_C1",
      TRUE ~ "Insufficient_evidence"
    )
  )

# Flag changes in direction across time without altering selection
trajectory_flags <- results_final %>%
  group_by(feature_ID) %>%
  summarise(
    sign_change = any(log2FC > 0, na.rm = TRUE) &
      any(log2FC < 0, na.rm = TRUE),
    .groups = "drop"
  )

candidate_results <- candidate_results %>%
  left_join(trajectory_flags, by = "feature_ID")

# Extract candidate lists
results_enriched <- candidate_results %>%
  dplyr::filter(direction == "Enriched_in_C1")

results_depleted <- candidate_results %>%
  dplyr::filter(direction == "Depleted_in_C1")

# Save all results and candidate lists
write_csv(
  candidate_results,
  file.path(base_dir, "Glmm", "average_treatment_results.csv")
)

write_csv(
  results_enriched,
  file.path(base_dir, "Glmm", "candidate_mutants_enriched.csv")
)

write_csv(
  results_depleted,
  file.path(base_dir, "Glmm", "candidate_mutants_depleted.csv")
)

# Summarize candidate classifications
candidate_results %>%
  dplyr::count(direction) %>%
  print()


# =========================
# 16. Prepare gene lists for DAVID
# =========================

# Clean identifiers and retain one copy per gene
clean_gene_ids <- function(ids) {
  ids <- trimws(as.character(ids))
  sort(unique(ids[
    !is.na(ids) &
      nzchar(ids) &
      !grepl("intergenic", ids, ignore.case = TRUE)
  ]))
}

# Extract unique candidate genes
depleted_genes <- clean_gene_ids(results_depleted$gene_ID)
enriched_genes <- clean_gene_ids(results_enriched$gene_ID)

# Identify genes selected in opposite directions across pools
conflicting_genes <- intersect(depleted_genes, enriched_genes)

# Define the background from features with valid average contrasts
background_genes <- candidate_results %>%
  dplyr::filter(contrast_valid) %>%
  dplyr::pull(gene_ID) %>%
  clean_gene_ids()

# Create the enrichment output directory
enrichment_dir <- file.path(base_dir, "Enrichmennt")

dir.create(
  enrichment_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# Save one gene identifier per line, without headers
writeLines(
  depleted_genes,
  file.path(enrichment_dir, "gene_ID_depleted.txt")
)

writeLines(
  enriched_genes,
  file.path(enrichment_dir, "gene_ID_enriched.txt")
)

writeLines(
  background_genes,
  file.path(enrichment_dir, "gene_ID_background.txt")
)

writeLines(
  conflicting_genes,
  file.path(enrichment_dir, "gene_ID_conflicting_directions.txt")
)

# Summarize exported lists
cat("Unique depleted genes:", length(depleted_genes), "\n")
cat("Unique enriched genes:", length(enriched_genes), "\n")
cat("Background genes:", length(background_genes), "\n")
cat("Genes selected in both directions:", length(conflicting_genes), "\n")