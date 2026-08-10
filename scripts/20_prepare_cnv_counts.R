#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

cat("\n")
cat("============================================================\n")
cat("PREPARING RAW COUNTS FOR CNV INFERENCE\n")
cat("============================================================\n\n")

cat("Seurat version: ",
    as.character(packageVersion("Seurat")),
    "\n")

cat("SeuratObject version: ",
    as.character(packageVersion("SeuratObject")),
    "\n\n")

project_dir <- "/home/katja/Cancer-Bioinformatics-Portfolio/03_SingleCell-RNAseq-Glioblastoma-Seurat"

seurat_file <- file.path(
  project_dir,
  "results/seurat_annotation/GBM_3samples_Seurat_annotated.rds"
)

gene_file <- file.path(
  project_dir,
  "results/cnv_annotation/GENCODE_v47_CNV_gene_annotation.csv"
)

output_dir <- file.path(
  project_dir,
  "results/cnv_input"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# LOAD SEURAT
# ------------------------------------------------------------

cat("Loading Seurat object...\n")

obj <- readRDS(seurat_file)

cat("Object loaded successfully.\n\n")

cat("Genes:", nrow(obj), "\n")
cat("Cells:", ncol(obj), "\n\n")

# ------------------------------------------------------------
# LOAD GENE ANNOTATION
# ------------------------------------------------------------

cat("Loading GENCODE annotation...\n")

gene_annot <- read.csv(
  gene_file,
  stringsAsFactors = FALSE
)

cat(
  "Annotated genes:",
  nrow(gene_annot),
  "\n\n"
)

# ------------------------------------------------------------
# JOIN LAYERS
# ------------------------------------------------------------

cat("Joining RNA layers...\n")

obj <- JoinLayers(
  obj,
  assay = "RNA"
)

cat("RNA layers:\n")
print(Layers(obj[["RNA"]]))
cat("\n")

# ------------------------------------------------------------
# EXTRACT RAW COUNTS
# ------------------------------------------------------------

cat("Extracting RAW counts...\n")

counts <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "counts"
)

cat(
  "Raw counts matrix:",
  nrow(counts),
  "genes x",
  ncol(counts),
  "cells\n\n"
)

# ------------------------------------------------------------
# MATCH GENES
# ------------------------------------------------------------

common_genes <- intersect(
  gene_annot$seurat_gene,
  rownames(counts)
)

cat(
  "Genes present in annotation and counts:",
  length(common_genes),
  "\n"
)

gene_annot <- gene_annot[
  gene_annot$seurat_gene %in% common_genes,
  ,
  drop = FALSE
]

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
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# CHROMOSOME ORDER
# ------------------------------------------------------------

gene_annot$chromosome <- factor(
  gene_annot$chromosome,
  levels = standard_chr,
  ordered = TRUE
)

# ------------------------------------------------------------
# GENOMIC ORDER
# ------------------------------------------------------------

gene_annot <- gene_annot[
  order(
    gene_annot$chromosome,
    gene_annot$start,
    gene_annot$end
  ),
  ,
  drop = FALSE
]

ordered_genes <- gene_annot$seurat_gene

# ------------------------------------------------------------
# SUBSET COUNTS
# ------------------------------------------------------------

counts <- counts[
  ordered_genes,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# CONSISTENCY CHECK
# ------------------------------------------------------------

if (!identical(
  rownames(counts),
  gene_annot$seurat_gene
)) {
  stop(
    "ERROR: counts and gene annotation are not aligned."
  )
}

cat("Gene order check: PASSED\n")

# ------------------------------------------------------------
# CELL METADATA
# ------------------------------------------------------------

metadata <- obj@meta.data

metadata <- metadata[
  colnames(counts),
  ,
  drop = FALSE
]

if (!identical(
  colnames(counts),
  rownames(metadata)
)) {
  stop(
    "ERROR: counts and metadata are not aligned."
  )
}

cat("Cell order check: PASSED\n\n")

# ------------------------------------------------------------
# BASIC COUNTS CHECK
# ------------------------------------------------------------

cat("Checking counts matrix...\n\n")

cat(
  "Minimum count:",
  min(counts),
  "\n"
)

cat(
  "Maximum count:",
  max(counts),
  "\n"
)

cat(
  "Total counts:",
  sum(counts),
  "\n\n"
)

# ------------------------------------------------------------
# DETECTION
# ------------------------------------------------------------

gene_detection <- Matrix::rowSums(
  counts > 0
)

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
  "\n\n"
)

# ------------------------------------------------------------
# SAMPLE COUNTS
# ------------------------------------------------------------

cat("Cells per sample:\n")

print(
  table(metadata$sample_id)
)

# ------------------------------------------------------------
# CELL TYPES
# ------------------------------------------------------------

cat("\nCells per cell type:\n")

print(
  sort(
    table(metadata$cell_type),
    decreasing = TRUE
  )
)

# ------------------------------------------------------------
# REFERENCE CANDIDATES
# ------------------------------------------------------------

reference_candidates <- c(
  "Astrocyte-like",
  "Oligodendrocyte",
  "Microglia-Macrophage",
  "Endothelial",
  "Mast-cell"
)

cat("\nPotential reference populations:\n\n")

for (ct in reference_candidates) {

  idx <- metadata$cell_type == ct

  cat(
    sprintf(
      "%-25s %5d\n",
      ct,
      sum(idx)
    )
  )
}

# ------------------------------------------------------------
# SAMPLE x REFERENCE
# ------------------------------------------------------------

cat("\nReference populations by sample:\n\n")

reference_idx <- metadata$cell_type %in% reference_candidates

print(
  table(
    metadata$sample_id[reference_idx],
    metadata$cell_type[reference_idx]
  )
)

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------

counts_file <- file.path(
  output_dir,
  "CNV_raw_counts_GENCODE_v47.rds"
)

gene_output <- file.path(
  output_dir,
  "CNV_gene_order_GENCODE_v47.csv"
)

metadata_output <- file.path(
  output_dir,
  "CNV_cell_metadata.csv"
)

cat("\nSaving raw counts...\n")

saveRDS(
  counts,
  counts_file
)

write.csv(
  gene_annot,
  gene_output,
  row.names = FALSE
)

write.csv(
  metadata,
  metadata_output,
  row.names = TRUE
)

# ------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------

summary_file <- file.path(
  output_dir,
  "CNV_raw_counts_summary.txt"
)

sink(summary_file)

cat("RAW COUNTS CNV INPUT\n")
cat("====================\n\n")

cat("Reference: GENCODE v47\n")
cat("Genome: GRCh38\n\n")

cat("Genes:", nrow(counts), "\n")
cat("Cells:", ncol(counts), "\n")

cat("\nGenes per chromosome:\n")
print(
  table(gene_annot$chromosome)
)

cat("\nCells per sample:\n")
print(
  table(metadata$sample_id)
)

cat("\nCells per cell type:\n")
print(
  sort(
    table(metadata$cell_type),
    decreasing = TRUE
  )
)

cat("\nReference populations by sample:\n")
print(
  table(
    metadata$sample_id[reference_idx],
    metadata$cell_type[reference_idx]
  )
)

sink()

# ------------------------------------------------------------
# FINISHED
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("RAW COUNTS PREPARATION COMPLETED\n")
cat("============================================================\n\n")

cat("Saved:\n")
cat(counts_file, "\n")
cat(gene_output, "\n")
cat(metadata_output, "\n")
cat(summary_file, "\n\n")

cat("Next step:\n")
cat("Evaluate reference populations before CNV inference.\n")
