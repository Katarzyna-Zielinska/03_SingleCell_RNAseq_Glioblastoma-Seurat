#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

cat("\n")
cat("============================================================\n")
cat("PREPARING CNV INPUT\n")
cat("============================================================\n\n")

cat("Seurat version: ", as.character(packageVersion("Seurat")), "\n")
cat("SeuratObject version: ", as.character(packageVersion("SeuratObject")), "\n\n")

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

project_dir <- "/home/katja/Cancer-Bioinformatics-Portfolio/03_SingleCell-RNAseq-Glioblastoma-Seurat"

seurat_file <- file.path(
  project_dir,
  "results/seurat_annotation/GBM_3samples_Seurat_annotated.rds"
)

annotation_file <- file.path(
  project_dir,
  "results/cnv_annotation/GENCODE_v47_CNV_gene_annotation.csv"
)

output_dir <- file.path(
  project_dir,
  "results/cnv_input"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Seurat object:\n", seurat_file, "\n\n")
cat("CNV annotation:\n", annotation_file, "\n\n")
cat("Output directory:\n", output_dir, "\n\n")

# ------------------------------------------------------------
# LOAD SEURAT OBJECT
# ------------------------------------------------------------

cat("Loading Seurat object...\n")

obj <- readRDS(seurat_file)

cat("Object loaded successfully.\n\n")

cat("Genes:", nrow(obj), "\n")
cat("Cells:", ncol(obj), "\n\n")

# ------------------------------------------------------------
# CHECK ANNOTATION
# ------------------------------------------------------------

required_cols <- c(
  "sample_id",
  "cell_type"
)

missing_cols <- setdiff(
  required_cols,
  colnames(obj@meta.data)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------
# LOAD GENE ANNOTATION
# ------------------------------------------------------------

cat("Loading GENCODE gene annotation...\n")

gene_annot <- read.csv(
  annotation_file,
  stringsAsFactors = FALSE
)

cat("Annotated genes:", nrow(gene_annot), "\n\n")

# ------------------------------------------------------------
# CHECK REQUIRED COLUMNS
# ------------------------------------------------------------

required_annot_cols <- c(
  "chromosome",
  "start",
  "end",
  "gene_id",
  "gene_name",
  "gene_type",
  "seurat_gene",
  "match_type"
)

missing_annot_cols <- setdiff(
  required_annot_cols,
  colnames(gene_annot)
)

if (length(missing_annot_cols) > 0) {
  stop(
    "Missing annotation columns: ",
    paste(missing_annot_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------
# MATCH SEURAT GENES
# ------------------------------------------------------------

seurat_genes <- rownames(obj)

gene_annot <- gene_annot[
  gene_annot$seurat_gene %in% seurat_genes,
]

cat("Genes from Seurat represented in annotation:",
    nrow(gene_annot), "\n")

# ------------------------------------------------------------
# REMOVE DUPLICATES
# ------------------------------------------------------------

if (anyDuplicated(gene_annot$seurat_gene) > 0) {

  cat("Duplicated Seurat genes detected.\n")
  cat("Removing duplicated entries...\n")

  gene_annot <- gene_annot[
    !duplicated(gene_annot$seurat_gene),
  ]
}

cat("Genes after duplicate removal:",
    nrow(gene_annot), "\n\n")

# ------------------------------------------------------------
# STANDARD CHROMOSOMES
# ------------------------------------------------------------

standard_chr <- c(
  paste0("chr", 1:22),
  "chrX",
  "chrY"
)

gene_annot <- gene_annot[
  gene_annot$chromosome %in% standard_chr,
]

cat(
  "Genes on standard chromosomes:",
  nrow(gene_annot),
  "\n\n"
)

# ------------------------------------------------------------
# CHROMOSOME ORDER
# ------------------------------------------------------------

gene_annot$chromosome <- factor(
  gene_annot$chromosome,
  levels = standard_chr,
  ordered = TRUE
)

# ------------------------------------------------------------
# SORT GENES BY GENOMIC POSITION
# ------------------------------------------------------------

gene_annot <- gene_annot[
  order(
    gene_annot$chromosome,
    gene_annot$start,
    gene_annot$end
  ),
]

# ------------------------------------------------------------
# REMOVE INVALID POSITIONS
# ------------------------------------------------------------

gene_annot <- gene_annot[
  !is.na(gene_annot$start) &
  !is.na(gene_annot$end) &
  gene_annot$start > 0 &
  gene_annot$end >= gene_annot$start,
]

cat(
  "Genes after coordinate filtering:",
  nrow(gene_annot),
  "\n\n"
)

# ------------------------------------------------------------
# EXPRESSION MATRIX
# ------------------------------------------------------------

cat("Extracting normalized expression matrix...\n")

DefaultAssay(obj) <- "RNA"

# Seurat v5 may contain multiple layers.
# JoinLayers creates a single expression layer
# suitable for downstream matrix extraction.

cat("Joining RNA layers...\n")

obj <- JoinLayers(
  obj,
  assay = "RNA"
)

cat("RNA layers after JoinLayers:\n")
print(Layers(obj[["RNA"]]))
cat("\n")

# ------------------------------------------------------------
# GET NORMALIZED DATA
# ------------------------------------------------------------

cat("Extracting RNA expression data...\n")

expr <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "data"
)

cat(
  "Expression matrix:",
  nrow(expr),
  "genes x",
  ncol(expr),
  "cells\n\n"
)

# ------------------------------------------------------------
# MATCH EXPRESSION GENES
# ------------------------------------------------------------

common_genes <- intersect(
  gene_annot$seurat_gene,
  rownames(expr)
)

cat(
  "Genes present in both annotation and expression:",
  length(common_genes),
  "\n"
)

if (length(common_genes) < 10000) {
  stop(
    "Too few genes available for CNV analysis: ",
    length(common_genes)
  )
}

# ------------------------------------------------------------
# SUBSET + ORDER EXPRESSION
# ------------------------------------------------------------

gene_annot <- gene_annot[
  gene_annot$seurat_gene %in% common_genes,
]

gene_annot <- gene_annot[
  order(
    gene_annot$chromosome,
    gene_annot$start,
    gene_annot$end
  ),
]

ordered_genes <- gene_annot$seurat_gene

expr <- expr[
  ordered_genes,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# FINAL CONSISTENCY CHECK
# ------------------------------------------------------------

if (!identical(
  rownames(expr),
  gene_annot$seurat_gene
)) {
  stop(
    "ERROR: expression matrix and gene annotation are not in identical order."
  )
}

cat("Final consistency check: PASSED\n\n")

# ------------------------------------------------------------
# CELL METADATA
# ------------------------------------------------------------

cell_metadata <- obj@meta.data

cell_metadata <- cell_metadata[
  colnames(expr),
  ,
  drop = FALSE
]

if (!identical(
  rownames(cell_metadata),
  colnames(expr)
)) {
  stop(
    "ERROR: cell metadata and expression matrix are not aligned."
  )
}

cat("Cell metadata alignment: PASSED\n\n")

# ------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------

cat("============================================================\n")
cat("FINAL CNV INPUT\n")
cat("============================================================\n\n")

cat("Genes:", nrow(expr), "\n")
cat("Cells:", ncol(expr), "\n")

cat(
  "Chromosomes:",
  length(unique(as.character(gene_annot$chromosome))),
  "\n\n"
)

cat("Genes per chromosome:\n")
print(table(gene_annot$chromosome))

cat("\nCells per sample:\n")
print(table(cell_metadata$sample_id))

cat("\nCells per cell type:\n")
print(sort(table(cell_metadata$cell_type), decreasing = TRUE))

# ------------------------------------------------------------
# SAVE EXPRESSION MATRIX
# ------------------------------------------------------------

cat("\nSaving expression matrix...\n")

saveRDS(
  expr,
  file.path(
    output_dir,
    "CNV_expression_matrix_normalized.rds"
  )
)

# ------------------------------------------------------------
# SAVE GENE ORDER
# ------------------------------------------------------------

write.csv(
  gene_annot,
  file.path(
    output_dir,
    "CNV_gene_order_GENCODE_v47.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# SAVE CELL METADATA
# ------------------------------------------------------------

write.csv(
  cell_metadata,
  file.path(
    output_dir,
    "CNV_cell_metadata.csv"
  ),
  row.names = TRUE
)

# ------------------------------------------------------------
# SAVE GENE LIST
# ------------------------------------------------------------

write.table(
  ordered_genes,
  file.path(
    output_dir,
    "CNV_ordered_gene_list.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# ------------------------------------------------------------
# SAVE SUMMARY
# ------------------------------------------------------------

summary_file <- file.path(
  output_dir,
  "CNV_input_summary.txt"
)

sink(summary_file)

cat("CNV INPUT SUMMARY\n")
cat("=================\n\n")

cat("Reference: GENCODE v47\n")
cat("Genome: GRCh38\n")
cat("Seurat genes:", length(seurat_genes), "\n")
cat("Final CNV genes:", nrow(expr), "\n")
cat("Cells:", ncol(expr), "\n\n")

cat("Genes per chromosome:\n")
print(table(gene_annot$chromosome))

cat("\nCells per sample:\n")
print(table(cell_metadata$sample_id))

cat("\nCells per cell type:\n")
print(sort(table(cell_metadata$cell_type), decreasing = TRUE))

sink()

# ------------------------------------------------------------
# FINISHED
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("CNV INPUT PREPARATION COMPLETED\n")
cat("============================================================\n\n")

cat("Saved files:\n\n")

cat(
  file.path(
    output_dir,
    "CNV_expression_matrix_normalized.rds"
  ),
  "\n"
)

cat(
  file.path(
    output_dir,
    "CNV_gene_order_GENCODE_v47.csv"
  ),
  "\n"
)

cat(
  file.path(
    output_dir,
    "CNV_cell_metadata.csv"
  ),
  "\n"
)

cat(
  file.path(
    output_dir,
    "CNV_ordered_gene_list.txt"
  ),
  "\n"
)

cat(
  file.path(
    output_dir,
    "CNV_input_summary.txt"
  ),
  "\n\n"
)

cat("Next step:\n")
cat("Validate the CNV input matrix before CNV inference.\n")
