#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

cat("\n")
cat("============================================================\n")
cat("CNV INPUT VALIDATION\n")
cat("============================================================\n\n")

project_dir <- "/home/katja/Cancer-Bioinformatics-Portfolio/03_SingleCell-RNAseq-Glioblastoma-Seurat"

input_dir <- file.path(
  project_dir,
  "results/cnv_input"
)

expr_file <- file.path(
  input_dir,
  "CNV_expression_matrix_normalized.rds"
)

gene_file <- file.path(
  input_dir,
  "CNV_gene_order_GENCODE_v47.csv"
)

metadata_file <- file.path(
  input_dir,
  "CNV_cell_metadata.csv"
)

# ------------------------------------------------------------
# LOAD FILES
# ------------------------------------------------------------

cat("Loading expression matrix...\n")

expr <- readRDS(expr_file)

cat("Expression matrix loaded.\n\n")

cat("Loading gene annotation...\n")

genes <- read.csv(
  gene_file,
  stringsAsFactors = FALSE
)

cat("Gene annotation loaded.\n\n")

cat("Loading cell metadata...\n")

metadata <- read.csv(
  metadata_file,
  row.names = 1,
  stringsAsFactors = FALSE
)

cat("Cell metadata loaded.\n\n")

# ------------------------------------------------------------
# BASIC DIMENSIONS
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("DIMENSIONS\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Expression:",
  nrow(expr),
  "genes x",
  ncol(expr),
  "cells\n"
)

cat(
  "Gene annotation:",
  nrow(genes),
  "genes\n"
)

cat(
  "Cell metadata:",
  nrow(metadata),
  "cells\n\n"
)

# ------------------------------------------------------------
# GENE ORDER CHECK
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("GENE ORDER CHECK\n")
cat("------------------------------------------------------------\n\n")

if (!identical(
  rownames(expr),
  genes$seurat_gene
)) {
  stop(
    "ERROR: Expression matrix genes do not match gene annotation order."
  )
}

cat("Gene order: PASSED\n")

# ------------------------------------------------------------
# CELL ORDER CHECK
# ------------------------------------------------------------

if (!identical(
  colnames(expr),
  rownames(metadata)
)) {
  stop(
    "ERROR: Expression matrix cells do not match metadata order."
  )
}

cat("Cell order: PASSED\n\n")

# ------------------------------------------------------------
# MISSING VALUES
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("MISSING VALUES\n")
cat("------------------------------------------------------------\n\n")

na_count <- sum(is.na(expr))

cat("NA values in expression matrix:", na_count, "\n")

if (na_count > 0) {
  warning(
    "Expression matrix contains NA values."
  )
} else {
  cat("Missing values: PASSED\n")
}

# ------------------------------------------------------------
# INFINITE VALUES
# ------------------------------------------------------------

inf_count <- sum(
  is.infinite(as.matrix(expr))
)

cat(
  "Infinite values in expression matrix:",
  inf_count,
  "\n"
)

if (inf_count > 0) {
  warning(
    "Expression matrix contains infinite values."
  )
} else {
  cat("Infinite values: PASSED\n")
}

# ------------------------------------------------------------
# EXPRESSION RANGE
# ------------------------------------------------------------

cat("\n")
cat("------------------------------------------------------------\n")
cat("EXPRESSION RANGE\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Minimum:",
  min(expr),
  "\n"
)

cat(
  "Maximum:",
  max(expr),
  "\n"
)

cat(
  "Mean:",
  mean(expr),
  "\n"
)

cat(
  "Median:",
  median(expr),
  "\n\n"
)

# ------------------------------------------------------------
# ZERO-EXPRESSION GENES
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("ZERO-EXPRESSION GENES\n")
cat("------------------------------------------------------------\n\n")

zero_genes <- rowSums(expr != 0) == 0

cat(
  "Genes with zero expression in all cells:",
  sum(zero_genes),
  "/",
  nrow(expr),
  "\n"
)

if (sum(zero_genes) > 0) {
  cat(
    "Percentage:",
    round(100 * mean(zero_genes), 2),
    "%\n"
  )
}

# ------------------------------------------------------------
# GENES WITH VERY LOW DETECTION
# ------------------------------------------------------------

gene_detection <- rowSums(expr > 0)

cat("\n")
cat("------------------------------------------------------------\n")
cat("GENE DETECTION\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Genes detected in >= 1 cell:",
  sum(gene_detection >= 1),
  "\n"
)

cat(
  "Genes detected in >= 10 cells:",
  sum(gene_detection >= 10),
  "\n"
)

cat(
  "Genes detected in >= 100 cells:",
  sum(gene_detection >= 100),
  "\n"
)

cat(
  "Genes detected in >= 1000 cells:",
  sum(gene_detection >= 1000),
  "\n"
)

# ------------------------------------------------------------
# CHROMOSOME CHECK
# ------------------------------------------------------------

cat("\n")
cat("------------------------------------------------------------\n")
cat("CHROMOSOME ORDER\n")
cat("------------------------------------------------------------\n\n")

print(
  table(
    factor(
      genes$chromosome,
      levels = c(
        paste0("chr", 1:22),
        "chrX",
        "chrY"
      )
    )
  )
)

# ------------------------------------------------------------
# CHECK MONOTONIC POSITIONS
# ------------------------------------------------------------

cat("\nChecking genomic positions within chromosomes...\n")

position_errors <- 0

for (chr in unique(genes$chromosome)) {

  idx <- which(
    genes$chromosome == chr
  )

  starts <- genes$start[idx]

  if (any(diff(starts) < 0)) {
    position_errors <- position_errors + 1

    cat(
      "Position order problem:",
      chr,
      "\n"
    )
  }
}

cat(
  "Chromosomes with position-order problems:",
  position_errors,
  "\n"
)

if (position_errors == 0) {
  cat("Genomic ordering: PASSED\n")
}

# ------------------------------------------------------------
# CELL METADATA
# ------------------------------------------------------------

cat("\n")
cat("------------------------------------------------------------\n")
cat("CELL METADATA\n")
cat("------------------------------------------------------------\n\n")

cat("Samples:\n")
print(table(metadata$sample_id))

cat("\nCell types:\n")
print(
  sort(
    table(metadata$cell_type),
    decreasing = TRUE
  )
)

# ------------------------------------------------------------
# SAMPLE x CELL TYPE
# ------------------------------------------------------------

cat("\n")
cat("Sample x cell type:\n")

print(
  table(
    metadata$sample_id,
    metadata$cell_type
  )
)

# ------------------------------------------------------------
# REFERENCE CANDIDATES\n
# ------------------------------------------------------------

cat("\n")
cat("------------------------------------------------------------\n")
cat("POTENTIAL NON-MALIGNANT REFERENCE POPULATIONS\n")
cat("------------------------------------------------------------\n\n")

reference_candidates <- c(
  "Astrocyte-like",
  "Oligodendrocyte",
  "Microglia-Macrophage",
  "Endothelial",
  "Mast-cell"
)

for (ct in reference_candidates) {

  n <- sum(
    metadata$cell_type == ct
  )

  cat(
    sprintf(
      "%-25s %5d cells\n",
      ct,
      n
    )
  )
}

# ------------------------------------------------------------
# SAVE VALIDATION SUMMARY
# ------------------------------------------------------------

summary_file <- file.path(
  input_dir,
  "CNV_input_validation_summary.txt"
)

sink(summary_file)

cat("CNV INPUT VALIDATION SUMMARY\n")
cat("============================\n\n")

cat("Genes:", nrow(expr), "\n")
cat("Cells:", ncol(expr), "\n")
cat("NA values:", na_count, "\n")
cat("Infinite values:", inf_count, "\n")
cat("Zero-expression genes:", sum(zero_genes), "\n")
cat("Position-order problems:", position_errors, "\n\n")

cat("Genes per chromosome:\n")
print(table(genes$chromosome))

cat("\nCells per sample:\n")
print(table(metadata$sample_id))

cat("\nCells per cell type:\n")
print(
  sort(
    table(metadata$cell_type),
    decreasing = TRUE
  )
)

sink()

# ------------------------------------------------------------
# FINISHED
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("CNV INPUT VALIDATION COMPLETED\n")
cat("============================================================\n\n")

cat(
  "Validation summary saved to:\n",
  summary_file,
  "\n\n"
)

cat("Next step:\n")
cat("Review validation results before CNV inference.\n")
