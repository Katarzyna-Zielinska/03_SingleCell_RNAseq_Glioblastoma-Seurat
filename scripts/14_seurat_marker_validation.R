#!/usr/bin/env Rscript

# ============================================================
# Seurat marker validation - Glioblastoma scRNA-seq
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(pheatmap)
})

cat("\n")
cat("============================================================\n")
cat("SEURAT MARKER VALIDATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. PROJECT DIRECTORIES
# ------------------------------------------------------------

project_dir <- "/home/katja/Cancer-Bioinformatics-Portfolio/03_SingleCell-RNAseq-Glioblastoma-Seurat"

input_file <- file.path(
  project_dir,
  "results/seurat_clustering/GBM_3samples_Seurat_clustered_UMAP.rds"
)

output_dir <- file.path(
  project_dir,
  "results/seurat_marker_validation"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Seurat version: ", as.character(packageVersion("Seurat")), "\n")
cat(
  "SeuratObject version: ",
  as.character(packageVersion("SeuratObject")),
  "\n\n"
)

# ------------------------------------------------------------
# 2. LOAD OBJECT
# ------------------------------------------------------------

cat("Loading Seurat object:\n")
cat(input_file, "\n\n")

obj <- readRDS(input_file)

cat("Object loaded successfully.\n\n")

cat("Genes: ", nrow(obj), "\n")
cat("Cells: ", ncol(obj), "\n\n")

cat("Clusters:\n")
print(table(Idents(obj)))

cat("\n")

# ------------------------------------------------------------
# 3. CANONICAL MARKER SETS
# ------------------------------------------------------------

marker_sets <- list(

  Astrocyte = c(
    "GFAP",
    "AQP4",
    "SLC1A3",
    "ALDH1L1",
    "SOX9",
    "S100B",
    "GLUL",
    "HOPX",
    "AQP1"
  ),

  OPC = c(
    "PDGFRA",
    "CSPG4",
    "OLIG1",
    "OLIG2",
    "SOX10",
    "GPR17",
    "BCAS1",
    "VCAN",
    "NKX2-2"
  ),

  Oligodendrocyte = c(
    "MBP",
    "MOG",
    "MAG",
    "MOBP",
    "PLP1",
    "CLDN11",
    "ERMN",
    "CNP",
    "OLIG2"
  ),

  Neuron = c(
    "RBFOX3",
    "MAP2",
    "SYT1",
    "SNAP25",
    "SYN1",
    "NEFL",
    "NEFM",
    "TUBB3",
    "ELAVL3",
    "ELAVL4"
  ),

  Microglia = c(
    "C1QA",
    "C1QB",
    "C1QC",
    "TMEM119",
    "P2RY12",
    "CX3CR1",
    "TREM2",
    "AIF1",
    "TYROBP"
  ),

  Macrophage = c(
    "CD68",
    "CTSD",
    "FCER1G",
    "LST1",
    "TYROBP",
    "CTSB",
    "LGALS3",
    "APOE",
    "CST3"
  ),

  Endothelial = c(
    "PECAM1",
    "VWF",
    "CLDN5",
    "KDR",
    "ESAM",
    "EMCN",
    "ENG",
    "RAMP2",
    "ADGRL4"
  ),

  Pericyte = c(
    "RGS5",
    "PDGFRB",
    "CSPG4",
    "MCAM",
    "COL4A1",
    "COL4A2",
    "ACTA2",
    "DES",
    "ABCC9"
  ),

  T_cell = c(
    "CD3D",
    "CD3E",
    "CD3G",
    "TRBC1",
    "TRBC2",
    "IL7R",
    "LTB",
    "CCL5",
    "LST1"
  ),

  NK_cell = c(
    "NKG7",
    "GNLY",
    "PRF1",
    "GZMB",
    "GZMH",
    "KLRD1",
    "FCGR3A"
  ),

  Mast_cell = c(
    "TPSAB1",
    "TPSB2",
    "KIT",
    "MS4A2",
    "CPA3",
    "HDC",
    "SLC18A2"
  ),

  Cycling = c(
    "MKI67",
    "TOP2A",
    "PCNA",
    "MCM2",
    "MCM4",
    "MCM6",
    "CCNB1",
    "CCNB2",
    "CDC20",
    "UBE2C",
    "AURKB",
    "CENPF"
  ),

  Hypoxia = c(
    "CA9",
    "HILPDA",
    "VEGFA",
    "ANGPTL4",
    "ADM",
    "DDIT4",
    "EGLN3",
    "BNIP3",
    "SLC2A1"
  ),

  Mesenchymal = c(
    "CHI3L1",
    "CHI3L2",
    "SERPINE1",
    "LGALS3",
    "TNC",
    "FN1",
    "COL1A1",
    "COL1A2",
    "VIM",
    "SPP1"
  )
)

# ------------------------------------------------------------
# 4. CHECK WHICH MARKERS ARE PRESENT
# ------------------------------------------------------------

all_markers <- unique(unlist(marker_sets))

genes_present <- intersect(all_markers, rownames(obj))
genes_missing <- setdiff(all_markers, rownames(obj))

cat("Canonical markers requested: ", length(all_markers), "\n")
cat("Markers present: ", length(genes_present), "\n")
cat("Markers missing: ", length(genes_missing), "\n\n")

if (length(genes_missing) > 0) {
  cat("Missing markers:\n")
  print(genes_missing)
  cat("\n")
}

write.csv(
  data.frame(
    gene = genes_missing
  ),
  file.path(output_dir, "missing_markers.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. SCORE EACH CELL
# ------------------------------------------------------------

cat("Calculating module scores...\n\n")

score_names <- character(0)

for (cell_type in names(marker_sets)) {

  markers <- intersect(marker_sets[[cell_type]], rownames(obj))

  if (length(markers) < 2) {
    warning(
      "Skipping ", cell_type,
      ": fewer than 2 markers present."
    )
    next
  }

  score_name <- paste0("score_", cell_type)

  obj <- AddModuleScore(
    object = obj,
    features = list(markers),
    name = score_name,
    search = FALSE
  )

  # AddModuleScore creates score_name + "1"
  actual_name <- paste0(score_name, "1")

  if (actual_name %in% colnames(obj@meta.data)) {
    colnames(obj@meta.data)[
      colnames(obj@meta.data) == actual_name
    ] <- score_name

    score_names <- c(score_names, score_name)
  }
}

cat("Scores calculated:\n")
print(score_names)

# ------------------------------------------------------------
# 6. MEAN SCORE PER CLUSTER
# ------------------------------------------------------------

cat("\nCalculating mean marker scores per cluster...\n\n")

score_data <- obj@meta.data[, score_names, drop = FALSE]

score_data$cluster <- as.character(Idents(obj))

cluster_scores <- aggregate(
  score_data[, score_names, drop = FALSE],
  by = list(cluster = score_data$cluster),
  FUN = mean
)

cluster_scores <- cluster_scores[
  order(as.numeric(as.character(cluster_scores$cluster))),
]

write.csv(
  cluster_scores,
  file.path(output_dir, "marker_scores_by_cluster.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. SAVE HIGHEST SCORING CELL TYPE PER CLUSTER
# ------------------------------------------------------------

score_matrix <- as.matrix(
  cluster_scores[, score_names, drop = FALSE]
)

rownames(score_matrix) <- cluster_scores$cluster

best_annotation <- apply(
  score_matrix,
  1,
  function(x) names(x)[which.max(x)]
)

best_score <- apply(
  score_matrix,
  1,
  max
)

annotation_table <- data.frame(
  cluster = rownames(score_matrix),
  predicted_marker_class = best_annotation,
  score = best_score,
  stringsAsFactors = FALSE
)

write.csv(
  annotation_table,
  file.path(output_dir, "preliminary_marker_annotation.csv"),
  row.names = FALSE
)

cat("Preliminary marker-based annotation:\n\n")
print(annotation_table)

# ------------------------------------------------------------
# 8. HEATMAP
# ------------------------------------------------------------

cat("\nGenerating marker score heatmap...\n")

pdf(
  file.path(
    output_dir,
    "MarkerScores_by_cluster_heatmap.pdf"
  ),
  width = 12,
  height = 9
)

pheatmap(
  score_matrix,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  border_color = NA,
  main = "Canonical marker scores by cluster"
)

dev.off()

png(
  file.path(
    output_dir,
    "MarkerScores_by_cluster_heatmap.png"
  ),
  width = 1800,
  height = 1400,
  res = 180
)

pheatmap(
  score_matrix,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  border_color = NA,
  main = "Canonical marker scores by cluster"
)

dev.off()

# ------------------------------------------------------------
# 9. DOTPLOT
# ------------------------------------------------------------

cat("Generating DotPlot...\n")

dotplot_genes <- unique(
  unlist(
    lapply(
      marker_sets,
      function(x) intersect(x, rownames(obj))
    )
  )
)

pdf(
  file.path(
    output_dir,
    "CanonicalMarkers_DotPlot.pdf"
  ),
  width = 18,
  height = 12
)

print(
  DotPlot(
    obj,
    features = dotplot_genes,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis() +
    ggtitle("Canonical marker expression across clusters")
)

dev.off()

png(
  file.path(
    output_dir,
    "CanonicalMarkers_DotPlot.png"
  ),
  width = 2600,
  height = 1700,
  res = 180
)

print(
  DotPlot(
    obj,
    features = dotplot_genes,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis() +
    ggtitle("Canonical marker expression across clusters")
)

dev.off()

# ------------------------------------------------------------
# 10. SAVE UPDATED OBJECT
# ------------------------------------------------------------

output_rds <- file.path(
  output_dir,
  "GBM_3samples_Seurat_marker_scored.rds"
)

saveRDS(
  obj,
  output_rds
)

cat("\n")
cat("============================================================\n")
cat("MARKER VALIDATION COMPLETED\n")
cat("============================================================\n\n")

cat("Saved files:\n")
cat(output_rds, "\n")
cat(file.path(output_dir, "marker_scores_by_cluster.csv"), "\n")
cat(file.path(output_dir, "preliminary_marker_annotation.csv"), "\n")
cat(file.path(output_dir, "MarkerScores_by_cluster_heatmap.png"), "\n")
cat(file.path(output_dir, "CanonicalMarkers_DotPlot.png"), "\n")

cat("\nImportant:\n")
cat("The marker-based labels are preliminary.\n")
cat("Do NOT rename clusters yet.\n")
cat("Next step: inspect the marker scores and DotPlot before final annotation.\n\n")
