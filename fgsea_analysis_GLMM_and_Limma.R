library(tidyverse)
library(fgsea)
library(ontologyIndex)

# =========================
# 1. Paths and analysis settings
# =========================

base_dir <- "/home/paula-di-santo/Documents/Post_analysis_meliloti"
annotation_dir <- file.path(base_dir, "GSEA_GLMM_vs_Limma", "Annotations")
go_dir <- file.path(annotation_dir, "GO_UniProt")
go_propagated_dir <- file.path(go_dir, "With_ancestors")
for (folder in c(annotation_dir, go_dir, go_propagated_dir)) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}

# Use the current complete average-contrast tables, not candidate-only tables
input_glmm <- file.path(base_dir, "Glmm", "average_treatment_results.csv")
input_limma <- file.path(base_dir, "Limma", "limma_contrasts.csv")

min_size <- 10L
max_size <- 500L
fdr_cutoff <- 0.05
random_seed <- 123L
run_direct_go <- TRUE

# Retain the complete evaluable feature set from each model
output_dir <- file.path(
  base_dir, "GSEA_GLMM_vs_Limma",
  paste0("All_features_minSize_", min_size, "_maxSize_", max_size)
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

annotation_files <- c(
  KEGG = file.path(annotation_dir, "KEGG_sme_gene_pathway_annotations.csv"),
  KEGG_names = file.path(annotation_dir, "KEGG_sme_pathway_names.csv"),
  GO_direct = file.path(annotation_dir, "GO_UniProt", "UniProt_1021_GO_annotations.csv"),
  GO_ancestors = file.path(
    annotation_dir, "GO_UniProt", "With_ancestors",
    "GO_annotations_with_ancestors.csv"
  )
)

# Check statistical input tables before preparing annotations
required_files <- c(input_glmm, input_limma)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing input files:\n", paste(missing_files, collapse = "\n"))
}

# Download only missing or empty files; promote completed downloads atomically
cached_download <- function(url, destination) {
  if (file.exists(destination) && file.info(destination)$size > 0) {
    message("Using existing file: ", basename(destination))
    return(invisible(destination))
  }

  temporary <- tempfile(tmpdir = dirname(destination))
  on.exit(unlink(temporary), add = TRUE)
  message("Downloading: ", basename(destination))

  download.file(url, destfile = temporary, mode = "wb")
  if (!file.exists(temporary) || file.info(temporary)$size == 0) {
    stop("Download produced an empty file: ", destination)
  }
  if (!file.rename(temporary, destination)) {
    stop("Could not save downloaded file: ", destination)
  }

  writeLines(
    c(paste("Download date:", Sys.Date()), paste("Source:", url)),
    paste0(destination, ".source.txt")
  )
  invisible(destination)
}


# =========================
# 2. Select evaluable average contrasts
# =========================

glmm <- read_csv(input_glmm, show_col_types = FALSE)
limma <- read_csv(input_limma, show_col_types = FALSE)

stopifnot(
  all(c("feature_ID", "gene_ID", "mean_log2FC", "SE_log2",
        "contrast_valid", "pvalue") %in% names(glmm)),
  all(c("feature_ID", "gene_ID", "contrast", "log2FC",
        "t_statistic", "pvalue") %in% names(limma))
)

# Calculate feature-level Wald Z statistics without significance filtering
glmm <- glmm %>%
  mutate(gene_ID = str_trim(gene_ID)) %>%
  filter(
    contrast_valid %in% TRUE,
    is.finite(mean_log2FC),
    is.finite(SE_log2), SE_log2 > 0,
    is.finite(pvalue),
    !is.na(gene_ID),
    !tolower(gene_ID) %in% c("", "na", "intergenic", "mapping")
  ) %>%
  mutate(wald_z = mean_log2FC / SE_log2) %>%
  filter(is.finite(wald_z))

# Use the same limma Average contrast and filtering as in the previous script
limma <- limma %>%
  mutate(gene_ID = str_trim(gene_ID)) %>%
  filter(
    contrast == "Average",
    is.finite(log2FC),
    is.finite(t_statistic),
    is.finite(pvalue),
    !is.na(gene_ID),
    !tolower(gene_ID) %in% c("", "na", "intergenic", "mapping")
  )

stopifnot(
  !anyDuplicated(glmm$feature_ID),
  !anyDuplicated(limma$feature_ID)
)

# Record feature overlap for information only; do not subset either model
common_features <- intersect(glmm$feature_ID, limma$feature_ID)
input_summary <- tibble(
  model = c("GLMM", "Limma"),
  evaluable_features = c(nrow(glmm), nrow(limma)),
  ranked_genes = c(n_distinct(glmm$gene_ID), n_distinct(limma$gene_ID)),
  shared_features = length(common_features),
  features_used = c(nrow(glmm), nrow(limma))
)
write_csv(input_summary, file.path(output_dir, "input_feature_summary.csv"))

# Check gene labels for shared features without restricting the analysis
feature_map <- inner_join(
  glmm %>% select(feature_ID, gene_GLMM = gene_ID),
  limma %>% select(feature_ID, gene_Limma = gene_ID),
  by = "feature_ID"
)
if (any(feature_map$gene_GLMM != feature_map$gene_Limma)) {
  stop("Gene identifiers disagree between models for shared features.")
}
write_csv(feature_map, file.path(output_dir, "feature_overlap_for_information.csv"))

# =========================
# 3. Build each model's complete gene-level ranking
# =========================

# Median statistics summarize genes; they are not combined gene-level tests
glmm_ranking <- glmm %>%
  group_by(gene_ID) %>%
  summarise(
    n_features = n(),
    ranking_score = median(wald_z),
    sign_change_between_pools = any(mean_log2FC > 0) & any(mean_log2FC < 0),
    .groups = "drop"
  ) %>%
  arrange(desc(ranking_score), gene_ID)

limma_ranking <- limma %>%
  group_by(gene_ID) %>%
  summarise(
    n_features = n(),
    ranking_score = median(t_statistic),
    sign_change_between_pools = any(log2FC > 0) & any(log2FC < 0),
    .groups = "drop"
  ) %>%
  arrange(desc(ranking_score), gene_ID)

rank_tables <- list(GLMM = glmm_ranking, Limma = limma_ranking)
ranks_by_model <- list()

for (model in names(rank_tables)) {
  model_dir <- file.path(output_dir, model)
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  tab <- rank_tables[[model]]
  ranks <- setNames(tab$ranking_score, tab$gene_ID)
  stopifnot(
    length(ranks) > 0,
    !anyDuplicated(names(ranks)),
    all(is.finite(ranks)),
    all(diff(ranks) <= 0) 
  )
  ranks_by_model[[model]] <- ranks
  write_csv(tab, file.path(model_dir, "gene_ranking.csv"))
  saveRDS(ranks, file.path(model_dir, "gene_ranking.rds"))
}

write_csv(glmm, file.path(output_dir, "GLMM_average_evaluable_features.csv"))
write_csv(limma, file.path(output_dir, "Limma_average_evaluable_features.csv"))

# Prepare annotation coverage over the union; each fgsea run uses its own full ranking
annotation_genes <- union(names(ranks_by_model$GLMM), names(ranks_by_model$Limma))

# =========================
# 4a. Reuse or download KEGG annotations
# =========================

kegg_annotation_file <- file.path(annotation_dir, "KEGG_sme_gene_pathway_annotations.csv")
kegg_names_file <- file.path(annotation_dir, "KEGG_sme_pathway_names.csv")

if (file.exists(kegg_annotation_file) && file.info(kegg_annotation_file)$size > 0) {
  kegg_annotations <- read_csv(kegg_annotation_file, show_col_types = FALSE)
} else {
  links_file <- file.path(annotation_dir, "KEGG_sme_links_raw.tsv")
  cached_download("https://rest.kegg.jp/link/pathway/sme", links_file)

  kegg_annotations <- read_tsv(
    links_file,
    col_names = c("kegg_gene_ID", "pathway_ID"),
    show_col_types = FALSE
  ) %>%
    mutate(
      gene_ID = str_remove(kegg_gene_ID, "^sme:"),
      gene_ID = str_replace(gene_ID, "^SM_([abc])", "SM\\1"),
      pathway_ID = str_remove(pathway_ID, "^path:")
    ) %>%
    distinct(pathway_ID, gene_ID, .keep_all = TRUE)

  write_csv(kegg_annotations, kegg_annotation_file)
}

if (file.exists(kegg_names_file) && file.info(kegg_names_file)$size > 0) {
  kegg_names <- read_csv(kegg_names_file, show_col_types = FALSE)
} else {
  names_file <- file.path(annotation_dir, "KEGG_sme_names_raw.tsv")
  cached_download("https://rest.kegg.jp/list/pathway/sme", names_file)

  kegg_names <- read_tsv(
    names_file,
    col_names = c("pathway_ID", "pathway_name"),
    show_col_types = FALSE
  ) %>%
    mutate(pathway_ID = str_remove(pathway_ID, "^path:"))

  write_csv(kegg_names, kegg_names_file)
}

writeLines(
  setdiff(annotation_genes, kegg_annotations$gene_ID),
  file.path(annotation_dir, "genes_without_KEGG_pathway_assignment.txt")
)

# =========================
# 4b. Reuse or download UniProt GO annotations
# =========================

uniprot_url <- paste0(
  "https://rest.uniprot.org/uniprotkb/stream?",
  "query=proteome%3AUP000001976",
  "&format=tsv&fields=accession,gene_names,go_p,go_f,go_c"
)
uniprot_file <- file.path(go_dir, "UniProt_1021_GO_raw.tsv")
cached_download(uniprot_url, uniprot_file)

uniprot_raw <- read_tsv(uniprot_file, show_col_types = FALSE)
stopifnot(ncol(uniprot_raw) == 5)
names(uniprot_raw) <- c("accession", "gene_names", "BP", "MF", "CC")

# Extract locus tags matching the identifiers used in the ranking
uniprot_genes <- uniprot_raw %>%
  mutate(
    gene_ID = stringr::str_extract_all(
      gene_names,
      "SM_?[abc]_?\\d+"
    )
  ) %>%
  tidyr::unnest_longer(gene_ID) %>%
  mutate(
    gene_ID = stringr::str_replace(
      gene_ID,
      "^SM_?([abc])_?",
      "SM\\1"
    )
  ) %>%
  distinct()

# Create one row per gene, GO term, and ontology
go_annotations <- uniprot_genes %>%
  pivot_longer(
    cols = c(BP, MF, CC),
    names_to = "ontology",
    values_to = "GO_annotation"
  ) %>%
  filter(!is.na(GO_annotation)) %>%
  separate_longer_delim(GO_annotation, delim = ";") %>%
  mutate(
    GO_ID = str_extract(GO_annotation, "GO:\\d{7}"),
    GO_name = str_remove(
      trimws(GO_annotation),
      "\\s*\\[GO:\\d{7}\\]$"
    )
  ) %>%
  filter(!is.na(GO_ID)) %>%
  select(accession, gene_ID, ontology, GO_ID, GO_name) %>%
  distinct()

# Restrict annotations to genes present in the ranking
go_ranked <- go_annotations %>%
  filter(gene_ID %in% annotation_genes)

# Summarize coverage for each GO ontology
go_coverage <- tibble(
  ontology = c("BP", "MF", "CC")
) %>%
  left_join(
    go_ranked %>%
      group_by(ontology) %>%
      summarise(
        annotated_genes = n_distinct(gene_ID),
        GO_terms = n_distinct(GO_ID),
        .groups = "drop"
      ),
    by = "ontology"
  ) %>%
  mutate(
    across(c(annotated_genes, GO_terms), ~ replace_na(.x, 0L)),
    ranked_genes = length(annotation_genes),
    coverage_percent = round(
      100 * annotated_genes / ranked_genes,
      1
    )
  )

# Save annotations, coverage, and download information
write_csv(
  go_annotations,
  file.path(go_dir, "UniProt_1021_GO_annotations.csv")
)

write_csv(
  go_coverage,
  file.path(go_dir, "GO_coverage.csv")
)

print(go_coverage)

# =========================
# 4c. Reuse the ontology and expand GO annotations
# =========================

obo_file <- file.path(go_propagated_dir, "go-basic.obo")
cached_download("https://purl.obolibrary.org/obo/go/go-basic.obo", obo_file)

# Reuse the parsed ontology when the underlying OBO file is unchanged
ontology_cache <- file.path(go_propagated_dir, "ontology_is_a_part_of.rds")
obo_md5 <- unname(tools::md5sum(obo_file))
cached_ontology <- if (file.exists(ontology_cache)) readRDS(ontology_cache) else NULL

if (!is.null(cached_ontology) && identical(cached_ontology$obo_md5, obo_md5)) {
  go <- cached_ontology$ontology
} else {
  go <- ontologyIndex::get_ontology(
    obo_file,
    propagate_relationships = c("is_a", "part_of"),
    extract_tags = "everything"
  )
  saveRDS(list(obo_md5 = obo_md5, ontology = go), ontology_cache)
}
rm(cached_ontology)

# Retain active ontology terms
active_ids <- unname(go$id[!go$obsolete])

# Map alternative identifiers to active primary identifiers
alternative_ids <- purrr::imap_dfr(
  go$alt_id,
  function(ids, primary_id) {
    tibble(
      original_GO_ID = as.character(ids),
      GO_ID = rep(primary_id, length(ids))
    )
  }
) %>%
  filter(GO_ID %in% active_ids)

id_map <- bind_rows(
  tibble(
    original_GO_ID = active_ids,
    GO_ID = active_ids
  ),
  alternative_ids
) %>%
  distinct()

stopifnot(!anyDuplicated(id_map$original_GO_ID))

# Match downloaded annotations to active ontology identifiers
direct_annotations <- go_annotations %>%
  filter(
    ontology %in% c("BP", "MF")
  ) %>%
  distinct(gene_ID, original_GO_ID = GO_ID) %>%
  left_join(id_map, by = "original_GO_ID")

# Record identifiers without an active primary or alternative match
unresolved_terms <- direct_annotations %>%
  filter(is.na(GO_ID)) %>%
  distinct(original_GO_ID)

write_csv(
  unresolved_terms,
  file.path(go_propagated_dir, "GO_unresolved_identifiers.csv")
)

# Expand each unique GO identifier once, then join to its genes
ancestor_map <- tibble(
  source_GO_ID = unique(na.omit(direct_annotations$GO_ID))
) %>%
  mutate(
    GO_ID = purrr::map(
      source_GO_ID,
      ~ unique(c(.x, go$ancestors[[.x]]))
    )
  ) %>%
  tidyr::unnest_longer(GO_ID) %>%
  filter(
    GO_ID %in% active_ids,
    !GO_ID %in% c("GO:0008150", "GO:0003674", "GO:0005575")
  )

expanded_annotations <- direct_annotations %>%
  filter(!is.na(GO_ID)) %>%
  distinct(gene_ID, source_GO_ID = GO_ID) %>%
  inner_join(ancestor_map, by = "source_GO_ID", relationship = "many-to-many") %>%
  distinct(gene_ID, GO_ID)

# Retrieve names and ontology membership from the downloaded ontology
term_information <- tibble(
  GO_ID = active_ids,
  GO_name = unname(go$name[active_ids]),
  namespace = vapply(
    active_ids,
    function(id) go$namespace[[id]][1],
    character(1)
  )
) %>%
  mutate(
    ontology = recode(
      namespace,
      biological_process = "BP",
      molecular_function = "MF",
      cellular_component = "CC",
      .default = NA_character_
    )
  )

go_expanded <- expanded_annotations %>%
  left_join(term_information, by = "GO_ID") %>%
  filter(ontology %in% c("BP", "MF")) %>%
  select(gene_ID, ontology, GO_ID, GO_name) %>%
  distinct()

# Save the expanded gene-term associations
write_csv(
  go_expanded,
  file.path(go_propagated_dir, "GO_annotations_with_ancestors.csv")
)

# Preserve ontology release information without changing the download date
obo_header <- readLines(obo_file, n = 40, warn = FALSE)
writeLines(
  c(
    "Propagation relationships: is_a, part_of",
    paste("OBO MD5:", obo_md5),
    obo_header[grepl("^(format-version|data-version|date):", obo_header)]
  ),
  file.path(go_propagated_dir, "GO_ontology_release.txt")
)
cat("Unresolved GO identifiers:", nrow(unresolved_terms), "\n")

# Include downloaded and processed annotations in the reproducibility manifest
required_files <- c(required_files, annotation_files)
if (!all(file.exists(required_files))) stop("Annotation preparation is incomplete.")

# Load the annotation collections prepared by this script

kegg <- read_csv(annotation_files[["KEGG"]], show_col_types = FALSE)
kegg_names <- read_csv(annotation_files[["KEGG_names"]], show_col_types = FALSE)
go_direct <- read_csv(annotation_files[["GO_direct"]], show_col_types = FALSE)
go_expanded <- read_csv(annotation_files[["GO_ancestors"]], show_col_types = FALSE)

collections <- list(
  KEGG = kegg %>%
    select(gene_ID, pathway_ID) %>%
    left_join(kegg_names, by = "pathway_ID") %>%
    transmute(set_ID = pathway_ID, gene_ID, term_name = pathway_name)
)

for (ont in c("BP", "MF")) {
  if (run_direct_go) {
    collections[[paste0("GO_", ont, "_direct")]] <- go_direct %>%
      filter(ontology == ont) %>%
      transmute(set_ID = GO_ID, gene_ID, term_name = GO_name)
  }
  collections[[paste0("GO_", ont, "_with_ancestors")]] <- go_expanded %>%
    filter(ontology == ont) %>%
    transmute(set_ID = GO_ID, gene_ID, term_name = GO_name)
}

# =========================
# 5. Run fgsea using the same settings for both models
# =========================

all_results <- list()
checks <- list()

for (collection in names(collections)) {
  # Reuse the same reference collection, intersecting memberships per model
  for (model in names(ranks_by_model)) {
    message("Running fgsea: ", model, " / ", collection)
    ranks <- ranks_by_model[[model]]

    annotations <- collections[[collection]] %>%
      filter(gene_ID %in% names(ranks)) %>%
      distinct(set_ID, gene_ID, term_name)

    term_names <- annotations %>% distinct(set_ID, term_name)
    stopifnot(!anyDuplicated(term_names$set_ID))
    gene_sets <- lapply(split(annotations$gene_ID, annotations$set_ID), unique)

    if (!any(lengths(gene_sets) >= min_size & lengths(gene_sets) <= max_size)) {
      stop("No eligible sets for ", model, " / ", collection)
    }

    set.seed(random_seed)

    result <- fgsea::fgseaMultilevel(
      pathways = gene_sets,
      stats = ranks,
      minSize = min_size,
      maxSize = max_size,
      scoreType = "std",
      eps = 0,
      nproc = 1
    ) %>%
      as_tibble() %>%
      left_join(term_names, by = c("pathway" = "set_ID")) %>%
      mutate(
        model = model,
        collection = collection,
        ranking_direction = case_when(
          NES > 0 ~ "Positive_C1_vs_DMSO",
          NES < 0 ~ "Negative_C1_vs_DMSO",
          TRUE ~ "Not_evaluable"
        )
      ) %>%
      arrange(padj, pval)

    destination <- file.path(output_dir, model, collection)
    dir.create(destination, recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, file.path(destination, "fgsea_results.rds"))

    export <- result %>%
      mutate(leadingEdge = map_chr(leadingEdge, ~ paste(.x, collapse = ";")))
    write_csv(export, file.path(destination, "fgsea_results.csv"))
    write_csv(
      export %>% filter(padj < fdr_cutoff),
      file.path(destination, "fgsea_significant.csv")
    )

    key <- paste(model, collection, sep = "__")
    all_results[[key]] <- export
    checks[[key]] <- tibble(
      model = model,
      collection = collection,
      ranked_genes = length(ranks),
      annotated_genes = n_distinct(annotations$gene_ID),
      tested_terms = nrow(result),
      evaluable_terms = sum(is.finite(result$padj)),
      significant_terms = sum(result$padj < fdr_cutoff, na.rm = TRUE)
    )
  }
}

# =========================
# 6. Export comparisons and reproducibility information
# =========================

combined <- bind_rows(all_results)
write_csv(combined, file.path(output_dir, "all_fgsea_results.csv"))

comparison <- combined %>%
  select(collection, pathway, term_name, model, size, NES, pval, padj) %>%
  pivot_wider(
    names_from = model,
    values_from = c(size, NES, pval, padj),
    names_glue = "{.value}_{model}"
  ) %>%
  mutate(
    evaluated_GLMM = is.finite(pval_GLMM),
    evaluated_Limma = is.finite(pval_Limma),
    significant_GLMM = padj_GLMM < fdr_cutoff,
    significant_Limma = padj_Limma < fdr_cutoff,
    same_NES_direction = sign(NES_GLMM) == sign(NES_Limma)
  )
write_csv(comparison, file.path(output_dir, "GLMM_vs_Limma_fgsea.csv"))

analysis_summary <- bind_rows(checks)
write_csv(analysis_summary, file.path(output_dir, "analysis_summary.csv"))

write_csv(
  tibble(
    minSize = min_size, maxSize = max_size,
    fdr_cutoff = fdr_cutoff, seed = random_seed,
    common_features = length(common_features),
    ranked_genes_GLMM = length(ranks_by_model$GLMM),
    ranked_genes_Limma = length(ranks_by_model$Limma),
    feature_selection = "All evaluable features independently per model",
    GLMM_ranking = "Median average-contrast Wald Z per gene",
    Limma_ranking = "Median Average-contrast moderated t per gene",
    p_adjustment = "BH separately within each model and annotation collection"
  ),
  file.path(output_dir, "analysis_parameters.csv")
)

write_csv(
  tibble(file = required_files, md5 = unname(tools::md5sum(required_files))),
  file.path(output_dir, "input_checksums.csv")
)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
print(input_summary)
print(analysis_summary, n = Inf)
message("Results saved in: ", output_dir)
