#!/usr/bin/env Rscript

# ============================================================
# 10_seurat_normalization.R
# Single-cell RNA-seq - Glioblastoma
#
# Step:
#   Normalization
#   Identification of highly variable genes
#   Scaling
#
# Input:
#   GBM_3samples_Seurat_QCfiltered.rds
#
# Output:
#   GBM_3samples_Seurat_normalized.rds
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
cat("Seurat NORMALIZATION AND VARIABLE FEATURES\n")
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
  "results/seurat_filter/GBM_3samples_Seurat_QCfiltered.rds"
)

output_dir <- file.path(
  project_dir,
  "results/seurat_normalization"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Load filtered object
# ------------------------------------------------------------

cat("Loading Seurat object:\n")
cat(input_file, "\n\n")

seurat_obj <- readRDS(input_file)

cat("Object loaded successfully.\n\n")


# ------------------------------------------------------------
# 4. Basic object information
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

print(
  table(seurat_obj$sample_id)
)

cat("\n")


# ------------------------------------------------------------
# 5. Check RNA assay
# ------------------------------------------------------------

DefaultAssay(seurat_obj) <- "RNA"

cat("Default assay:\n")
cat(DefaultAssay(seurat_obj), "\n\n")


# ------------------------------------------------------------
# 6. Normalize RNA expression
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("NORMALIZATION\n")
cat("------------------------------------------------------------\n\n")

cat("Running NormalizeData()...\n")

seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

cat("\nNormalization completed.\n\n")


# ------------------------------------------------------------
# 7. Identify highly variable genes
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("HIGHLY VARIABLE GENES\n")
cat("------------------------------------------------------------\n\n")

cat("Running FindVariableFeatures()...\n")

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)

cat("\nVariable feature selection completed.\n\n")


# ------------------------------------------------------------
# 8. Extract variable genes
# ------------------------------------------------------------

variable_genes <- VariableFeatures(
  seurat_obj
)

cat(
  "Number of variable genes: ",
  length(variable_genes),
  "\n\n",
  sep = ""
)

cat("Top 20 variable genes:\n\n")

print(
  head(variable_genes, 20)
)

cat("\n")


# ------------------------------------------------------------
# 9. Save variable genes
# ------------------------------------------------------------

write.table(
  variable_genes,
  file = file.path(
    output_dir,
    "highly_variable_genes.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


# ------------------------------------------------------------
# 10. Plot highly variable genes
# ------------------------------------------------------------

cat("Generating variable-feature plot...\n")

p_hvg <- VariableFeaturePlot(
  seurat_obj
)

p_hvg_labeled <- LabelPoints(
  plot = p_hvg,
  points = head(variable_genes, 10),
  repel = TRUE
)

ggsave(
  filename = file.path(
    output_dir,
    "highly_variable_genes.png"
  ),
  plot = p_hvg_labeled,
  width = 10,
  height = 7,
  dpi = 300
)

cat("Variable-feature plot saved.\n\n")


# ------------------------------------------------------------
# 11. Scale data
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("SCALING\n")
cat("------------------------------------------------------------\n\n")

cat("Running ScaleData()...\n")

seurat_obj <- ScaleData(
  seurat_obj,
  features = rownames(seurat_obj),
  verbose = TRUE
)

cat("\nScaling completed.\n\n")


# ------------------------------------------------------------
# 12. Check assay layers
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("RNA ASSAY LAYERS\n")
cat("------------------------------------------------------------\n\n")

print(
  Layers(seurat_obj[["RNA"]])
)

cat("\n")


# ------------------------------------------------------------
# 13. Save normalized object
# ------------------------------------------------------------

output_file <- file.path(
  output_dir,
  "GBM_3samples_Seurat_normalized.rds"
)

saveRDS(
  seurat_obj,
  output_file
)

cat("Saved normalized Seurat object:\n")
cat(output_file, "\n\n")


# ------------------------------------------------------------
# 14. Final summary
# ------------------------------------------------------------

cat("============================================================\n")
cat("NORMALIZATION COMPLETED SUCCESSFULLY\n")
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
  "Highly variable genes: ",
  length(VariableFeatures(seurat_obj)),
  "\n\n",
  sep = ""
)

cat("Normalization:\n")
cat("  Method: LogNormalize\n")
cat("  Scale factor: 10000\n\n")

cat("Variable feature selection:\n")
cat("  Method: vst\n")
cat("  Number of features: 2000\n\n")

cat("Output directory:\n")
cat(output_dir, "\n\n")

cat("Next step:\n")
cat("PCA and dimensionality reduction.\n\n")

cat("============================================================\n")
cat("DONE\n")
cat("============================================================\n")
