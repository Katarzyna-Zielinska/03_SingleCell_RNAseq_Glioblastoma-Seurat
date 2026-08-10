#!/usr/bin/env Rscript

# ============================================================
# 12_seurat_clustering_umap.R
# Single-cell RNA-seq - Glioblastoma
#
# PCA -> Neighbors -> Clustering -> UMAP
#
# Input:
#   GBM_3samples_Seurat_PCA.rds
#
# Main parameters:
#   PCs: 1:20
#   Resolution: 0.5
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
})


cat("\n")
cat("============================================================\n")
cat("Seurat CLUSTERING + UMAP\n")
cat("============================================================\n\n")

cat(
  "Seurat version: ",
  as.character(packageVersion("Seurat")),
  "\n"
)

cat(
  "SeuratObject version: ",
  as.character(packageVersion("SeuratObject")),
  "\n\n"
)


# ------------------------------------------------------------
# 2. Project paths
# ------------------------------------------------------------

project_dir <- "/home/katja/Cancer-Bioinformatics-Portfolio/03_SingleCell-RNAseq-Glioblastoma-Seurat"

input_file <- file.path(
  project_dir,
  "results/seurat_pca/GBM_3samples_Seurat_PCA.rds"
)

output_dir <- file.path(
  project_dir,
  "results/seurat_clustering"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Parameters
# ------------------------------------------------------------

dims_use <- 1:20
resolution_use <- 0.5


# ------------------------------------------------------------
# 4. Load object
# ------------------------------------------------------------

cat("Loading Seurat object:\n")
cat(input_file, "\n\n")

seurat_obj <- readRDS(input_file)

cat("Object loaded successfully.\n\n")


# ------------------------------------------------------------
# 5. Basic information
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("INPUT OBJECT\n")
cat("------------------------------------------------------------\n")

cat(
  "Genes: ",
  nrow(seurat_obj),
  "\n",
  sep = ""
)

cat(
  "Cells: ",
  ncol(seurat_obj),
  "\n\n",
  sep = ""
)

cat("Cells per sample:\n")
print(table(seurat_obj$sample_id))
cat("\n")


# ------------------------------------------------------------
# 6. Check PCA
# ------------------------------------------------------------

if (!"pca" %in% Reductions(seurat_obj)) {
  stop("ERROR: PCA reduction not found.")
}

available_pcs <- ncol(
  Embeddings(
    seurat_obj,
    reduction = "pca"
  )
)

cat(
  "Available PCs: ",
  available_pcs,
  "\n",
  sep = ""
)

cat(
  "PCs used downstream: 1:",
  max(dims_use),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Find neighbors
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("NEIGHBOR GRAPH\n")
cat("------------------------------------------------------------\n\n")

cat("Running FindNeighbors()...\n")

seurat_obj <- FindNeighbors(
  seurat_obj,
  reduction = "pca",
  dims = dims_use,
  verbose = TRUE
)

cat("\nFindNeighbors completed successfully.\n\n")


# ------------------------------------------------------------
# 8. Find clusters
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("CLUSTERING\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Resolution: ",
  resolution_use,
  "\n\n",
  sep = ""
)

cat("Running FindClusters()...\n")

seurat_obj <- FindClusters(
  seurat_obj,
  resolution = resolution_use,
  verbose = TRUE
)

cat("\nFindClusters completed successfully.\n\n")


# ------------------------------------------------------------
# 9. Cluster summary
# ------------------------------------------------------------

cluster_column <- grep(
  "^RNA_snn_res\\.",
  colnames(seurat_obj@meta.data),
  value = TRUE
)

if (length(cluster_column) == 0) {
  stop("ERROR: Cluster metadata column not found.")
}

cluster_column <- cluster_column[length(cluster_column)]

cat("Cluster metadata column:\n")
cat(cluster_column, "\n\n")

cat("Cells per cluster:\n")

cluster_counts <- table(
  seurat_obj@meta.data[[cluster_column]]
)

print(cluster_counts)

cat("\n")


# ------------------------------------------------------------
# 10. Cluster proportions
# ------------------------------------------------------------

cluster_proportions <- prop.table(
  cluster_counts
) * 100

cluster_summary <- data.frame(
  cluster = names(cluster_counts),
  cells = as.integer(cluster_counts),
  percent = as.numeric(cluster_proportions)
)

write.csv(
  cluster_summary,
  file.path(
    output_dir,
    "cluster_sizes.csv"
  ),
  row.names = FALSE
)

cat("Cluster summary saved.\n\n")


# ------------------------------------------------------------
# 11. Cluster composition by sample
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("CLUSTER COMPOSITION BY SAMPLE\n")
cat("------------------------------------------------------------\n\n")

cluster_sample_table <- table(
  seurat_obj@meta.data[[cluster_column]],
  seurat_obj$sample_id
)

print(cluster_sample_table)

write.csv(
  as.data.frame(cluster_sample_table),
  file.path(
    output_dir,
    "cluster_by_sample_counts.csv"
  ),
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 12. UMAP
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("UMAP\n")
cat("------------------------------------------------------------\n\n")

cat("Running RunUMAP()...\n")

seurat_obj <- RunUMAP(
  seurat_obj,
  reduction = "pca",
  dims = dims_use,
  verbose = TRUE
)

cat("\nUMAP completed successfully.\n\n")


# ------------------------------------------------------------
# 13. UMAP by cluster
# ------------------------------------------------------------

cat("Generating UMAP by cluster...\n")

p_umap_cluster <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = cluster_column,
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    paste0(
      "UMAP - Clusters (resolution ",
      resolution_use,
      ")"
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "UMAP_clusters.png"
  ),
  plot = p_umap_cluster,
  width = 10,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------
# 14. UMAP by sample
# ------------------------------------------------------------

cat("Generating UMAP by sample...\n")

p_umap_sample <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "sample_id"
) +
  ggtitle(
    "UMAP - Sample of origin"
  )

ggsave(
  filename = file.path(
    output_dir,
    "UMAP_by_sample.png"
  ),
  plot = p_umap_sample,
  width = 10,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------
# 15. UMAP split by sample
# ------------------------------------------------------------

cat("Generating split UMAP...\n")

p_umap_split <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = cluster_column,
  split.by = "sample_id",
  label = TRUE,
  repel = TRUE,
  ncol = 3
) +
  ggtitle(
    "UMAP - Clusters split by sample"
  )

ggsave(
  filename = file.path(
    output_dir,
    "UMAP_clusters_split_by_sample.png"
  ),
  plot = p_umap_split,
  width = 15,
  height = 5,
  dpi = 300
)


# ------------------------------------------------------------
# 16. Save metadata
# ------------------------------------------------------------

metadata_file <- file.path(
  output_dir,
  "cell_metadata_clustering.csv"
)

write.csv(
  seurat_obj@meta.data,
  metadata_file,
  row.names = TRUE
)

cat("Cell metadata saved:\n")
cat(metadata_file, "\n\n")


# ------------------------------------------------------------
# 17. Save Seurat object
# ------------------------------------------------------------

output_file <- file.path(
  output_dir,
  "GBM_3samples_Seurat_clustered_UMAP.rds"
)

saveRDS(
  seurat_obj,
  output_file
)

cat("Saved clustered Seurat object:\n")
cat(output_file, "\n\n")


# ------------------------------------------------------------
# 18. Final summary
# ------------------------------------------------------------

cat("============================================================\n")
cat("CLUSTERING + UMAP COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat(
  "Cells: ",
  ncol(seurat_obj),
  "\n",
  sep = ""
)

cat(
  "Genes: ",
  nrow(seurat_obj),
  "\n",
  sep = ""
)

cat(
  "PCs used: 1:",
  max(dims_use),
  "\n",
  sep = ""
)

cat(
  "Clustering resolution: ",
  resolution_use,
  "\n",
  sep = ""
)

cat(
  "Number of clusters: ",
  length(cluster_counts),
  "\n\n",
  sep = ""
)

cat("Output directory:\n")
cat(output_dir, "\n\n")

cat("Important:\n")
cat(
  "Clusters are NOT biologically annotated yet.\n",
  "The next step is marker-gene identification and ",
  "biological annotation of each cluster.\n\n",
  sep = ""
)

cat("============================================================\n")
cat("DONE\n")
cat("============================================================\n")
