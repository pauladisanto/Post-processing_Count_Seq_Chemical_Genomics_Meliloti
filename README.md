# Post-processing_Count_Seq_Chemical_Genomics_Meliloti
Analysis pipeline for identifying treatment-specific mutant responses in time-course Tn-seq data using DESeq2 normalization and GLMM modeling.

# Overview
This repository contains an analysis pipeline for identifying treatment-specific mutant responses in time-course sequencing data.
The workflow processes Salmon quantifications, annotates mutants, normalizes counts, and applies a generalized linear mixed model (GLMM) to detect mutants that behave differently under treatment (C1) compared to control (DMSO).

# Biological objective
The goal is to identify mutants whose abundance:

changes over time
differs between C1 and DMSO conditions
reflects treatment-specific fitness effects

Mutants with similar behavior under C1 and DMSO are considered non-informative.

# Output files
count_matrix.csv → feature-level counts
metadata_deseq2.csv → sample metadata
glmmSeq_results.csv → model output
top_mutants_ranked.csv → ranked treatment-specific mutants
