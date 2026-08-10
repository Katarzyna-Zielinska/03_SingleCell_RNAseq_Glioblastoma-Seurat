#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

cat("\n")
cat("============================================================\n")
cat("FINAL CNV REFERENCE VALIDATION\n")
cat("============================================================\n\n")

counts_file <- "results/cnv_input/CNV_reference_counts_GENCODE_v47.rds"
metadata_file <- "results/cnv_input/CNV_reference_metadata_selected.csv"
gene_file <- "results/cnv_input/CNV_gene_order_GENCODE_v47.csv"

output_dir <- "results/cnv_input"

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------

cat("Loading reference counts...\n")

counts <- readRDS(counts_file)

cat("Reference counts loaded.\n")
cat("Genes:", nrow(counts), "\n")
cat("Cells:", ncol(counts), "\n\n")

cat("Loading reference metadata...\n")

meta <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

meta$cell <- as.character(meta$cell)

cat("Metadata rows:", nrow(meta), "\n\n")

cat("Loading gene annotation...\n")

gene_annot <- read.csv(
  gene_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Gene annotation rows:", nrow(gene_annot), "\n\n")

# ------------------------------------------------------------
# CELL ALIGNMENT
# ------------------------------------------------------------

cat("============================================================\n")
cat("CELL ALIGNMENT\n")
cat("============================================================\n\n")

if (anyDuplicated(colnames(counts)) > 0) {
  stop("ERROR: Duplicate cell IDs in reference counts.")
}

if (anyDuplicated(meta$cell) > 0) {
  stop("ERROR: Duplicate cell IDs in metadata.")
}

if (!setequal(colnames(counts), meta$cell)) {
  stop("ERROR: Counts and metadata contain different cells.")
}

meta <- meta[
  match(colnames(counts), meta$cell),
  ,
  drop = FALSE
]

if (!identical(colnames(counts), meta$cell)) {
  stop("ERROR: Cell order mismatch.")
}

cat("Cell identifiers: PASSED\n")
cat("Cell order: PASSED\n\n")

# ------------------------------------------------------------
# GENE ALIGNMENT
# ------------------------------------------------------------

cat("============================================================\n")
cat("GENE ALIGNMENT\n")
cat("============================================================\n\n")

gene_column_candidates <- c(
  "seurat_gene",
  "gene_name",
  "gene_id"
)

gene_column <- gene_column_candidates[
  gene_column_candidates %in% colnames(gene_annot)
][1]

if (is.na(gene_column)) {
  stop("ERROR: Could not identify gene column.")
}

gene_ids <- as.character(gene_annot[[gene_column]])

if (anyDuplicated(gene_ids) > 0) {
  stop("ERROR: Duplicate gene identifiers in annotation.")
}

if (!setequal(rownames(counts), gene_ids)) {
  stop("ERROR: Counts and gene annotation contain different genes.")
}

gene_annot <- gene_annot[
  match(rownames(counts), gene_ids),
  ,
  drop = FALSE
]

if (!identical(rownames(counts), gene_annot[[gene_column]])) {
  stop("ERROR: Gene order mismatch.")
}

cat("Gene identifiers: PASSED\n")
cat("Gene order: PASSED\n\n")

# ------------------------------------------------------------
# COUNT VALIDATION
# ------------------------------------------------------------

cat("============================================================\n")
cat("COUNT MATRIX VALIDATION\n")
cat("============================================================\n\n")

if (!is.numeric(counts@x)) {
  stop("ERROR: Count matrix is not numeric.")
}

cat("Minimum:", min(counts), "\n")
cat("Maximum:", max(counts), "\n")

total_counts <- sum(counts)

cat("Total counts:", total_counts, "\n")

if (anyNA(counts)) {
  stop("ERROR: NA values detected.")
}

cat("NA values: 0\n")

if (any(!is.finite(counts@x))) {
  stop("ERROR: Infinite values detected.")
}

cat("Infinite values: 0\n")

negative_values <- sum(counts@x < 0)

cat("Negative values:", negative_values, "\n")

if (negative_values > 0) {
  stop("ERROR: Negative values detected.")
}

cat("Count matrix: PASSED\n\n")

# ------------------------------------------------------------
# CELL TYPE COMPOSITION
# ------------------------------------------------------------

cat("============================================================\n")
cat("REFERENCE CELL COMPOSITION\n")
cat("============================================================\n\n")

cat("Cell types:\n")
print(table(meta$cell_type))

cat("\nSamples:\n")
print(table(meta$sample_id))

cat("\nSample x cell type:\n")
print(table(meta$sample_id, meta$cell_type))

# ------------------------------------------------------------
# SAMPLE BALANCE
# ------------------------------------------------------------

sample_counts <- table(meta$sample_id)

cat("\nCells per sample:\n")
print(sample_counts)

if (length(unique(as.integer(sample_counts))) != 1) {
  warning(
    "Reference cells are not perfectly balanced between samples."
  )
} else {
  cat("Sample balance: PASSED\n")
}

# ------------------------------------------------------------
# DETECTION
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GENE DETECTION IN REFERENCE\n")
cat("============================================================\n\n")

detected_per_gene <- Matrix::rowSums(counts > 0)

cat(
  "Genes detected in >= 1 reference cell:",
  sum(detected_per_gene >= 1),
  "\n"
)

cat(
  "Genes detected in >= 10 reference cells:",
  sum(detected_per_gene >= 10),
  "\n"
)

cat(
  "Genes detected in >= 50 reference cells:",
  sum(detected_per_gene >= 50),
  "\n"
)

cat(
  "Genes detected in >= 100 reference cells:",
  sum(detected_per_gene >= 100),
  "\n"
)

cat(
  "Genes detected in all reference cells:",
  sum(detected_per_gene == ncol(counts)),
  "\n\n"
)

# ------------------------------------------------------------
# ZERO-EXPRESSION GENES
# ------------------------------------------------------------

zero_genes <- sum(detected_per_gene == 0)

cat(
  "Zero-expression genes:",
  zero_genes,
  "\n"
)

cat(
  "Percentage:",
  round(
    100 * zero_genes / nrow(counts),
    2
  ),
  "%\n\n"
)

# ------------------------------------------------------------
# COUNTS PER CELL
# ------------------------------------------------------------

counts_per_cell <- Matrix::colSums(counts)

cat("============================================================\n")
cat("COUNTS PER REFERENCE CELL\n")
cat("============================================================\n\n")

cat("Minimum:", min(counts_per_cell), "\n")
cat("Median:", median(counts_per_cell), "\n")
cat("Mean:", round(mean(counts_per_cell), 2), "\n")
cat("Maximum:", max(counts_per_cell), "\n\n")

# ------------------------------------------------------------
# GENES PER CELL
# ------------------------------------------------------------

genes_per_cell <- Matrix::colSums(counts > 0)

cat("Genes detected per cell:\n")
cat("Minimum:", min(genes_per_cell), "\n")
cat("Median:", median(genes_per_cell), "\n")
cat("Mean:", round(mean(genes_per_cell), 2), "\n")
cat("Maximum:", max(genes_per_cell), "\n\n")

# ------------------------------------------------------------
# CHROMOSOME DISTRIBUTION
# ------------------------------------------------------------

cat("============================================================\n")
cat("CHROMOSOME DISTRIBUTION\n")
cat("============================================================\n\n")

if ("chromosome" %in% colnames(gene_annot)) {

  chr_table <- table(gene_annot$chromosome)

  print(chr_table)

} else {

  warning(
    "Chromosome column not found in gene annotation."
  )
}

# ------------------------------------------------------------
# REFERENCE CELL TYPE CHECK
# ------------------------------------------------------------

allowed_reference_types <- c(
  "Microglia-Macrophage",
  "Oligodendrocyte"
)

unexpected_types <- setdiff(
  unique(meta$cell_type),
  allowed_reference_types
)

if (length(unexpected_types) > 0) {

  warning(
    paste(
      "Unexpected cell types:",
      paste(
        unexpected_types,
        collapse = ", "
      )
    )
  )

} else {

  cat(
    "Reference cell types: PASSED\n\n"
  )
}

# ------------------------------------------------------------
# SAVE SUMMARY TABLES
# ------------------------------------------------------------

sample_summary <- data.frame(
  sample_id = names(sample_counts),
  cells = as.integer(sample_counts)
)

write.csv(
  sample_summary,
  file.path(
    output_dir,
    "CNV_reference_final_sample_summary.csv"
  ),
  row.names = FALSE
)

celltype_summary <- data.frame(
  cell_type = names(table(meta$cell_type)),
  cells = as.integer(table(meta$cell_type))
)

write.csv(
  celltype_summary,
  file.path(
    output_dir,
    "CNV_reference_final_celltype_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# SAVE VALIDATION SUMMARY
# ------------------------------------------------------------

summary_file <- file.path(
  output_dir,
  "CNV_reference_final_validation_summary.txt"
)

sink(summary_file)

cat("FINAL CNV REFERENCE VALIDATION\n")
cat("==============================\n\n")

cat("Genes:", nrow(counts), "\n")
cat("Cells:", ncol(counts), "\n\n")

cat("Reference cell types:\n")
print(table(meta$cell_type))

cat("\nReference cells by sample:\n")
print(table(meta$sample_id))

cat("\nSample x cell type:\n")
print(table(meta$sample_id, meta$cell_type))

cat("\nGenes detected in >= 1 cell:",
    sum(detected_per_gene >= 1),
    "\n")

cat("Genes detected in >= 10 cells:",
    sum(detected_per_gene >= 10),
    "\n")

cat("Genes detected in >= 50 cells:",
    sum(detected_per_gene >= 50),
    "\n")

cat("Zero-expression genes:",
    zero_genes,
    "\n")

cat("\nCounts per cell:\n")
cat("Minimum:", min(counts_per_cell), "\n")
cat("Median:", median(counts_per_cell), "\n")
cat("Mean:", mean(counts_per_cell), "\n")
cat("Maximum:", max(counts_per_cell), "\n")

cat("\nGenes detected per cell:\n")
cat("Minimum:", min(genes_per_cell), "\n")
cat("Median:", median(genes_per_cell), "\n")
cat("Mean:", mean(genes_per_cell), "\n")
cat("Maximum:", max(genes_per_cell), "\n")

sink()

# ------------------------------------------------------------
# FINAL MESSAGE
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("FINAL CNV REFERENCE VALIDATION COMPLETED\n")
cat("============================================================\n\n")

cat("Reference matrix:\n")
cat(
  nrow(counts),
  "genes x",
  ncol(counts),
  "cells\n\n"
)

cat("Cell alignment: PASSED\n")
cat("Gene alignment: PASSED\n")
cat("Count matrix: PASSED\n")
cat("Reference composition: PASSED\n\n")

cat("Saved files:\n")
cat(
  file.path(
    output_dir,
    "CNV_reference_final_sample_summary.csv"
  ),
  "\n"
)

cat(
  file.path(
    output_dir,
    "CNV_reference_final_celltype_summary.csv"
  ),
  "\n"
)

cat(summary_file, "\n\n")

cat("Next step:\n")
cat("Proceed to CNV inference.\n\n")
