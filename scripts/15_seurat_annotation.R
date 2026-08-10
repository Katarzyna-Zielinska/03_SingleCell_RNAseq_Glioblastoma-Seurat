#!/usr/bin/env Rscript

# ============================================================
# Seurat - Cell Type Annotation
# Glioblastoma scRNA-seq
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

cat("\n============================================================\n")
cat("SEURAT CELL TYPE ANNOTATION\n")
cat("============================================================\n\n")

cat("Seurat version: ", as.character(packageVersion("Seurat")), "\n")
cat("SeuratObject version: ", as.character(packageVersion("SeuratObject")), "\n\n")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

project_dir <- "/home/katja/Cancer-Bioinformatics-Portfolio/03_SingleCell-RNAseq-Glioblastoma-Seurat"

input_file <- file.path(
  project_dir,
  "results/seurat_marker_validation/GBM_3samples_Seurat_marker_scored.rds"
)

output_dir <- file.path(
  project_dir,
  "results/seurat_annotation"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_rds <- file.path(
  output_dir,
  "GBM_3samples_Seurat_annotated.rds"
)

# ------------------------------------------------------------
# Load object
# ------------------------------------------------------------

cat("Loading Seurat object:\n")
cat(input_file, "\n\n")

obj <- readRDS(input_file)

cat("Object loaded successfully.\n\n")

cat("Genes:", nrow(obj), "\n")
cat("Cells:", ncol(obj), "\n\n")

# ------------------------------------------------------------
# Check clusters
# ------------------------------------------------------------

cat("Original clusters:\n\n")
print(table(Idents(obj)))

# ------------------------------------------------------------
# Define biological annotations
# ------------------------------------------------------------

cluster_annotation <- c(
  "0"  = "Neural-like",
  "1"  = "Hypoxic-Mesenchymal-like",
  "2"  = "Mesenchymal-like",
  "3"  = "Cycling-S-phase",
  "4"  = "Neural-like-mixed",
  "5"  = "OPC-like",
  "6"  = "Vascular-Mesenchymal-like",
  "7"  = "OPC-Neural-like",
  "8"  = "Microglia-Macrophage",
  "9"  = "Cycling-G2M",
  "10" = "Cycling-Mixed",
  "11" = "Astrocyte-like",
  "12" = "Highly-Cycling",
  "13" = "Endothelial",
  "14" = "Mast-cell",
  "15" = "Oligodendrocyte"
)

# ------------------------------------------------------------
# Verify all clusters have annotations
# ------------------------------------------------------------

clusters <- as.character(Idents(obj))

missing_annotations <- setdiff(
  unique(clusters),
  names(cluster_annotation)
)

if (length(missing_annotations) > 0) {
  stop(
    paste(
      "Missing annotations for clusters:",
      paste(missing_annotations, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# Add cell type annotation
# ------------------------------------------------------------

obj$cell_type <- unname(
  cluster_annotation[as.character(Idents(obj))]
)

# Make cell_type an ordered factor according to cluster order
obj$cell_type <- factor(
  obj$cell_type,
  levels = unique(cluster_annotation)
)

# ------------------------------------------------------------
# Print annotation table
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CELL TYPE ANNOTATION\n")
cat("============================================================\n\n")

annotation_table <- data.frame(
  cluster = names(cluster_annotation),
  cell_type = unname(cluster_annotation),
  stringsAsFactors = FALSE
)

print(annotation_table, row.names = FALSE)

# ------------------------------------------------------------
# Cell counts
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CELLS PER CELL TYPE\n")
cat("============================================================\n\n")

cells_per_type <- as.data.frame(
  table(obj$cell_type)
)

colnames(cells_per_type) <- c(
  "cell_type",
  "cells"
)

print(cells_per_type, row.names = FALSE)

write.csv(
  cells_per_type,
  file.path(output_dir, "cells_per_cell_type.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Cell type by sample
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CELL TYPE BY SAMPLE\n")
cat("============================================================\n\n")

celltype_sample <- as.data.frame(
  table(
    sample_id = obj$sample_id,
    cell_type = obj$cell_type
  )
)

print(celltype_sample, row.names = FALSE)

write.csv(
  celltype_sample,
  file.path(output_dir, "cell_type_by_sample.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Save annotation table
# ------------------------------------------------------------

write.csv(
  annotation_table,
  file.path(output_dir, "cluster_annotation.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# UMAP - cell types
# ------------------------------------------------------------

cat("\nGenerating cell type UMAP...\n")

p_umap <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  ggtitle("Glioblastoma scRNA-seq - Cell Type Annotation") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

ggsave(
  filename = file.path(
    output_dir,
    "UMAP_cell_type_annotation.png"
  ),
  plot = p_umap,
  width = 12,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# UMAP - samples
# ------------------------------------------------------------

cat("Generating sample UMAP...\n")

p_sample <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "sample_id",
  raster = FALSE
) +
  ggtitle("Glioblastoma scRNA-seq - Samples") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(
    output_dir,
    "UMAP_by_sample.png"
  ),
  plot = p_sample,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# UMAP - original clusters
# ------------------------------------------------------------

cat("Generating cluster UMAP...\n")

p_cluster <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  ggtitle("Glioblastoma scRNA-seq - Seurat Clusters") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(
    output_dir,
    "UMAP_seurat_clusters.png"
  ),
  plot = p_cluster,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Save annotated object
# ------------------------------------------------------------

saveRDS(
  obj,
  output_rds
)

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ANNOTATION COMPLETED\n")
cat("============================================================\n\n")

cat("Cells:", ncol(obj), "\n")
cat("Genes:", nrow(obj), "\n")
cat("Cell types:", length(unique(obj$cell_type)), "\n\n")

cat("Saved files:\n")
cat(output_rds, "\n")
cat(file.path(output_dir, "cluster_annotation.csv"), "\n")
cat(file.path(output_dir, "cells_per_cell_type.csv"), "\n")
cat(file.path(output_dir, "cell_type_by_sample.csv"), "\n")
cat(file.path(output_dir, "UMAP_cell_type_annotation.png"), "\n")
cat(file.path(output_dir, "UMAP_by_sample.png"), "\n")
cat(file.path(output_dir, "UMAP_seurat_clusters.png"), "\n\n")

cat("Important:\n")
cat("These annotations are marker-based and preliminary.\n")
cat("Malignant versus non-malignant status has NOT yet been established.\n")
cat("Do not interpret Neural-like/OPC-like/Mesenchymal-like as definitive tumor states yet.\n\n")

cat("Next step:\n")
cat("Evaluate malignant versus non-malignant populations using CNV/state analysis.\n\n")
