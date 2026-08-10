#!/usr/bin/env Rscript

# ============================================================
# 09_seurat_filter.R
# Single-cell RNA-seq - Glioblastoma
# QC filtering of Seurat object
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

project_dir <- "/home/katja/Cancer-Bioinformatics-Portfolio/03_SingleCell-RNAseq-Glioblastoma-Seurat"

input_file <- file.path(
  project_dir,
  "results/seurat_qc/GBM_3samples_Seurat_raw.rds"
)

output_dir <- file.path(
  project_dir,
  "results/seurat_filter"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. QC thresholds
# ------------------------------------------------------------

min_nCount <- 1000
min_nFeature <- 500
max_percent.mt <- 20

# ------------------------------------------------------------
# 3. Load Seurat object
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("Seurat QC FILTERING\n")
cat("============================================================\n\n")

cat("Seurat version: ", as.character(packageVersion("Seurat")), "\n")
cat("SeuratObject version: ", as.character(packageVersion("SeuratObject")), "\n\n")

cat("Loading:\n")
cat(input_file, "\n\n")

seurat_obj <- readRDS(input_file)

cat("Object loaded successfully.\n\n")

# ------------------------------------------------------------
# 4. Check metadata
# ------------------------------------------------------------

required_columns <- c(
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt",
  "sample_id"
)

missing_columns <- setdiff(
  required_columns,
  colnames(seurat_obj@meta.data)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 5. BEFORE QC summary
# ------------------------------------------------------------

metadata_before <- seurat_obj@meta.data

cells_before <- ncol(seurat_obj)

cat("------------------------------------------------------------\n")
cat("BEFORE QC\n")
cat("------------------------------------------------------------\n")

cat("Total cells:", cells_before, "\n\n")

before_table <- as.data.frame(
  table(metadata_before$sample_id)
)

colnames(before_table) <- c(
  "sample_id",
  "cells_before"
)

print(before_table)

# ------------------------------------------------------------
# 6. Apply QC filters
# ------------------------------------------------------------

cat("\n")
cat("------------------------------------------------------------\n")
cat("QC FILTERS\n")
cat("------------------------------------------------------------\n")

cat("Minimum nCount_RNA:   ", min_nCount, "\n")
cat("Minimum nFeature_RNA: ", min_nFeature, "\n")
cat("Maximum percent.mt:   ", max_percent.mt, "%\n\n")

pass_qc <- (
  metadata_before$nCount_RNA >= min_nCount &
  metadata_before$nFeature_RNA >= min_nFeature &
  metadata_before$percent.mt <= max_percent.mt
)

cat("Cells passing QC:", sum(pass_qc), "\n")
cat("Cells removed:", sum(!pass_qc), "\n\n")

# ------------------------------------------------------------
# 7. Save list of removed cells
# ------------------------------------------------------------

removed_cells <- rownames(metadata_before)[!pass_qc]

removed_df <- data.frame(
  cell = removed_cells,
  sample_id = metadata_before[removed_cells, "sample_id"],
  nCount_RNA = metadata_before[removed_cells, "nCount_RNA"],
  nFeature_RNA = metadata_before[removed_cells, "nFeature_RNA"],
  percent.mt = metadata_before[removed_cells, "percent.mt"],
  stringsAsFactors = FALSE
)

write.csv(
  removed_df,
  file.path(output_dir, "removed_cells.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Filter Seurat object
# ------------------------------------------------------------

seurat_filtered <- subset(
  seurat_obj,
  cells = rownames(metadata_before)[pass_qc]
)

# ------------------------------------------------------------
# 9. AFTER QC summary
# ------------------------------------------------------------

metadata_after <- seurat_filtered@meta.data

cells_after <- ncol(seurat_filtered)

cat("------------------------------------------------------------\n")
cat("AFTER QC\n")
cat("------------------------------------------------------------\n")

cat("Total cells:", cells_after, "\n\n")

after_table <- as.data.frame(
  table(metadata_after$sample_id)
)

colnames(after_table) <- c(
  "sample_id",
  "cells_after"
)

print(after_table)

# ------------------------------------------------------------
# 10. Combined before/after table
# ------------------------------------------------------------

summary_table <- merge(
  before_table,
  after_table,
  by = "sample_id",
  all = TRUE
)

summary_table$cells_removed <-
  summary_table$cells_before -
  summary_table$cells_after

summary_table$percent_retained <-
  100 * summary_table$cells_after /
  summary_table$cells_before

summary_table$percent_removed <-
  100 * summary_table$cells_removed /
  summary_table$cells_before

write.csv(
  summary_table,
  file.path(output_dir, "qc_filter_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Save filtered Seurat object
# ------------------------------------------------------------

output_rds <- file.path(
  output_dir,
  "GBM_3samples_Seurat_QCfiltered.rds"
)

saveRDS(
  seurat_filtered,
  output_rds
)

cat("\n")
cat("Saved filtered Seurat object:\n")
cat(output_rds, "\n")

# ------------------------------------------------------------
# 12. QC plots BEFORE filtering
# ------------------------------------------------------------

p_before <- VlnPlot(
  seurat_obj,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "sample_id",
  pt.size = 0,
  ncol = 3
)

ggsave(
  file.path(output_dir, "QC_before_filtering.png"),
  p_before,
  width = 15,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 13. QC plots AFTER filtering
# ------------------------------------------------------------

p_after <- VlnPlot(
  seurat_filtered,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "sample_id",
  pt.size = 0,
  ncol = 3
)

ggsave(
  file.path(output_dir, "QC_after_filtering.png"),
  p_after,
  width = 15,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 14. Scatter plot AFTER filtering
# ------------------------------------------------------------

p_scatter <- FeatureScatter(
  seurat_filtered,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

ggsave(
  file.path(output_dir, "QC_nCount_vs_nFeature_after.png"),
  p_scatter,
  width = 7,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 15. Final summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("QC FILTERING COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat("Cells before QC:   ", cells_before, "\n")
cat("Cells after QC:    ", cells_after, "\n")
cat("Cells removed:     ", cells_before - cells_after, "\n")
cat(
  "Cells retained:    ",
  round(100 * cells_after / cells_before, 2),
  "%\n\n"
)

cat("QC thresholds:\n")
cat("  nCount_RNA   >= ", min_nCount, "\n")
cat("  nFeature_RNA >= ", min_nFeature, "\n")
cat("  percent.mt   <= ", max_percent.mt, "%\n\n")

cat("Output directory:\n")
cat(output_dir, "\n\n")

cat("Next step:\n")
cat("Normalization and identification of highly variable genes.\n\n")
