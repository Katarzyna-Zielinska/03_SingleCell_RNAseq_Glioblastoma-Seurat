#!/usr/bin/env Rscript

# ============================================================
# 11_seurat_pca.R
# Single-cell RNA-seq - Glioblastoma
#
# Step:
#   Principal Component Analysis (PCA)
#
# Input:
#   GBM_3samples_Seurat_normalized.rds
#
# Output:
#   GBM_3samples_Seurat_PCA.rds
#   PCA plots
#   PCA variance table
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
cat("Seurat PCA ANALYSIS\n")
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
  "results/seurat_normalization/GBM_3samples_Seurat_normalized.rds"
)

output_dir <- file.path(
  project_dir,
  "results/seurat_pca"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Load normalized object
# ------------------------------------------------------------

cat("Loading Seurat object:\n")
cat(input_file, "\n\n")

seurat_obj <- readRDS(input_file)

cat("Object loaded successfully.\n\n")


# ------------------------------------------------------------
# 4. Basic information
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
# 5. Default assay
# ------------------------------------------------------------

DefaultAssay(seurat_obj) <- "RNA"

cat("Default assay:\n")
cat(DefaultAssay(seurat_obj), "\n\n")


# ------------------------------------------------------------
# 6. Check variable features
# ------------------------------------------------------------

variable_genes <- VariableFeatures(seurat_obj)

cat("------------------------------------------------------------\n")
cat("VARIABLE FEATURES\n")
cat("------------------------------------------------------------\n")

cat(
  "Number of variable genes: ",
  length(variable_genes),
  "\n\n",
  sep = ""
)

if (length(variable_genes) == 0) {
  stop(
    "ERROR: No variable features found. ",
    "Run the normalization/HVG step first."
  )
}


# ------------------------------------------------------------
# 7. Run PCA
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("PRINCIPAL COMPONENT ANALYSIS\n")
cat("------------------------------------------------------------\n\n")

cat("Running RunPCA()...\n")

seurat_obj <- RunPCA(
  seurat_obj,
  features = variable_genes,
  npcs = 50,
  verbose = TRUE
)

cat("\nPCA completed successfully.\n\n")


# ------------------------------------------------------------
# 8. PCA information
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("PCA INFORMATION\n")
cat("------------------------------------------------------------\n\n")

print(
  seurat_obj[["pca"]]
)

cat("\n")


# ------------------------------------------------------------
# 9. PCA loadings
# ------------------------------------------------------------

cat("Top genes contributing to the first PCs:\n\n")

for (pc in 1:5) {

  cat("PC", pc, ":\n")

  print(
    head(
      Loadings(
        seurat_obj[["pca"]]
      )[, pc, drop = FALSE],
      10
    )
  )

  cat("\n")
}


# ------------------------------------------------------------
# 10. PCA standard deviations
# ------------------------------------------------------------

pca_stdev <- Stdev(
  seurat_obj[["pca"]]
)

pca_variance <- data.frame(
  PC = seq_along(pca_stdev),
  StandardDeviation = pca_stdev,
  Variance = pca_stdev^2,
  PercentVariance = 100 *
    pca_stdev^2 /
    sum(pca_stdev^2)
)

pca_variance$CumulativeVariance <- cumsum(
  pca_variance$PercentVariance
)


# ------------------------------------------------------------
# 11. Save PCA variance table
# ------------------------------------------------------------

write.csv(
  pca_variance,
  file.path(
    output_dir,
    "PCA_variance.csv"
  ),
  row.names = FALSE
)

cat("PCA variance table saved.\n\n")


# ------------------------------------------------------------
# 12. Elbow plot
# ------------------------------------------------------------

cat("Generating ElbowPlot...\n")

p_elbow <- ElbowPlot(
  seurat_obj,
  ndims = 50
)

ggsave(
  filename = file.path(
    output_dir,
    "PCA_ElbowPlot.png"
  ),
  plot = p_elbow,
  width = 9,
  height = 6,
  dpi = 300
)

cat("Elbow plot saved.\n\n")


# ------------------------------------------------------------
# 13. PCA plot: PC1 vs PC2
# ------------------------------------------------------------

cat("Generating PCA plots...\n")

p_pca_12 <- DimPlot(
  seurat_obj,
  reduction = "pca",
  dims = c(1, 2),
  group.by = "sample_id"
) +
  ggtitle(
    "PCA: PC1 vs PC2 by sample"
  )

ggsave(
  filename = file.path(
    output_dir,
    "PCA_PC1_PC2_by_sample.png"
  ),
  plot = p_pca_12,
  width = 9,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# 14. PCA plot: PC2 vs PC3
# ------------------------------------------------------------

p_pca_23 <- DimPlot(
  seurat_obj,
  reduction = "pca",
  dims = c(2, 3),
  group.by = "sample_id"
) +
  ggtitle(
    "PCA: PC2 vs PC3 by sample"
  )

ggsave(
  filename = file.path(
    output_dir,
    "PCA_PC2_PC3_by_sample.png"
  ),
  plot = p_pca_23,
  width = 9,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# 15. PCA plot without sample grouping
# ------------------------------------------------------------

p_pca_uncolored <- DimPlot(
  seurat_obj,
  reduction = "pca",
  dims = c(1, 2),
  group.by = "sample_id"
) +
  ggtitle(
    "PCA: sample distribution"
  )

ggsave(
  filename = file.path(
    output_dir,
    "PCA_sample_distribution.png"
  ),
  plot = p_pca_uncolored,
  width = 9,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# 16. Save PCA object
# ------------------------------------------------------------

output_file <- file.path(
  output_dir,
  "GBM_3samples_Seurat_PCA.rds"
)

saveRDS(
  seurat_obj,
  output_file
)

cat("Saved PCA Seurat object:\n")
cat(output_file, "\n\n")


# ------------------------------------------------------------
# 17. Final summary
# ------------------------------------------------------------

cat("============================================================\n")
cat("PCA COMPLETED SUCCESSFULLY\n")
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
  "Variable genes used for PCA: ",
  length(variable_genes),
  "\n",
  sep = ""
)

cat(
  "Number of PCs calculated: 50\n\n"
)

cat("Output directory:\n")
cat(output_dir, "\n\n")

cat("Important:\n")
cat(
  "The number of PCs for downstream clustering/UMAP ",
  "will be selected after inspecting the ElbowPlot.\n\n",
  sep = ""
)

cat("Next step:\n")
cat(
  "Inspect PCA variance and select the number of PCs ",
  "for neighborhood graph construction.\n\n",
  sep = ""
)

cat("============================================================\n")
cat("DONE\n")
cat("============================================================\n")
