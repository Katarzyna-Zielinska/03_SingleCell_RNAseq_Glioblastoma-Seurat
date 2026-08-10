#!/usr/bin/env Rscript

# ============================================================
# 13_seurat_markers.R
# Single-cell RNA-seq - Glioblastoma
#
# Step:
#   Join Seurat v5 assay layers
#   Identify marker genes for all clusters
#
# Input:
#   GBM_3samples_Seurat_clustered_UMAP.rds
#
# Output:
#   all_cluster_markers.csv
#   top20_markers_per_cluster.csv
#   cluster_marker_heatmap.png
#   broad_cell_type_markers_DotPlot.png
#   GBM_3samples_Seurat_markers.rds
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
cat("Seurat MARKER GENE ANALYSIS\n")
cat("============================================================\n\n")

cat(
  "Seurat version: ",
  as.character(packageVersion("Seurat")),
  "\n",
  sep = ""
)

cat(
  "SeuratObject version: ",
  as.character(packageVersion("SeuratObject")),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 2. Project paths
# ------------------------------------------------------------

project_dir <- "/home/katja/Cancer-Bioinformatics-Portfolio/03_SingleCell-RNAseq-Glioblastoma-Seurat"

input_file <- file.path(
  project_dir,
  "results/seurat_clustering/GBM_3samples_Seurat_clustered_UMAP.rds"
)

output_dir <- file.path(
  project_dir,
  "results/seurat_markers"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Parameters
# ------------------------------------------------------------

min_pct_use <- 0.25
logfc_threshold_use <- 0.25

top_n_markers <- 20
heatmap_top_n <- 10


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

cat("Cells per cluster:\n")
print(table(Idents(seurat_obj)))
cat("\n")


# ------------------------------------------------------------
# 6. Set RNA assay
# ------------------------------------------------------------

DefaultAssay(seurat_obj) <- "RNA"

cat("Default assay:\n")
cat(DefaultAssay(seurat_obj), "\n\n")


# ------------------------------------------------------------
# 7. Check clusters
# ------------------------------------------------------------

clusters <- levels(Idents(seurat_obj))

if (length(clusters) < 2) {
  stop(
    "ERROR: Fewer than two clusters detected."
  )
}

cat("Clusters detected:\n")
print(clusters)
cat("\n")


# ------------------------------------------------------------
# 8. Inspect RNA layers BEFORE joining
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("RNA ASSAY LAYERS BEFORE JoinLayers()\n")
cat("------------------------------------------------------------\n\n")

layers_before <- Layers(
  seurat_obj[["RNA"]]
)

print(layers_before)

cat("\n")


# ------------------------------------------------------------
# 9. Join Seurat v5 layers
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("JOINING RNA ASSAY LAYERS\n")
cat("------------------------------------------------------------\n\n")

cat("Running JoinLayers()...\n")

seurat_obj <- JoinLayers(
  seurat_obj,
  assay = "RNA"
)

cat("JoinLayers() completed.\n\n")


# ------------------------------------------------------------
# 10. Inspect RNA layers AFTER joining
# ------------------------------------------------------------

cat("RNA ASSAY LAYERS AFTER JoinLayers():\n")

layers_after <- Layers(
  seurat_obj[["RNA"]]
)

print(layers_after)

cat("\n")


# ------------------------------------------------------------
# 11. Verify required layers
# ------------------------------------------------------------

if (!"data" %in% layers_after) {

  stop(
    "ERROR: Joined RNA assay does not contain a 'data' layer."
  )

}

if (!"counts" %in% layers_after) {

  warning(
    "WARNING: Joined RNA assay does not contain a 'counts' layer."
  )

}


# ------------------------------------------------------------
# 12. Re-confirm cluster identities
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("CLUSTER IDENTITIES\n")
cat("------------------------------------------------------------\n\n")

cat("Active identities:\n")
print(levels(Idents(seurat_obj)))

cat("\nCells per cluster:\n")
print(table(Idents(seurat_obj)))

cat("\n")


# ------------------------------------------------------------
# 13. Find all markers
# ------------------------------------------------------------

cat("============================================================\n")
cat("FINDING MARKER GENES\n")
cat("============================================================\n\n")

cat(
  "Minimum percent expressed: ",
  min_pct_use,
  "\n",
  sep = ""
)

cat(
  "Minimum log2 fold-change: ",
  logfc_threshold_use,
  "\n\n",
  sep = ""
)

cat("Running FindAllMarkers()...\n")
cat("This may take several minutes.\n\n")


markers_all <- FindAllMarkers(
  object = seurat_obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = min_pct_use,
  logfc.threshold = logfc_threshold_use,
  test.use = "wilcox",
  verbose = TRUE
)


cat("\n")
cat("FindAllMarkers completed.\n\n")


# ------------------------------------------------------------
# 14. Check marker results
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("MARKER SUMMARY\n")
cat("------------------------------------------------------------\n\n")

if (nrow(markers_all) == 0) {

  stop(
    paste(
      "ERROR: FindAllMarkers returned zero markers.",
      "Check RNA layers and assay configuration."
    )
  )

}

cat(
  "Total marker rows: ",
  nrow(markers_all),
  "\n\n",
  sep = ""
)

cat("Columns returned by FindAllMarkers:\n")
print(colnames(markers_all))

cat("\n")


# ------------------------------------------------------------
# 15. Detect fold-change column
# ------------------------------------------------------------

if ("avg_log2FC" %in% colnames(markers_all)) {

  fc_column <- "avg_log2FC"

} else if ("avg_logFC" %in% colnames(markers_all)) {

  fc_column <- "avg_logFC"

} else {

  stop(
    "ERROR: No log fold-change column found in marker results."
  )

}


cat(
  "Fold-change column: ",
  fc_column,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 16. Markers per cluster
# ------------------------------------------------------------

cat("Markers per cluster:\n")

marker_counts <- table(
  markers_all$cluster
)

print(marker_counts)

cat("\n")


# ------------------------------------------------------------
# 17. Save ALL markers
# ------------------------------------------------------------

all_markers_file <- file.path(
  output_dir,
  "all_cluster_markers.csv"
)

write.csv(
  markers_all,
  all_markers_file,
  row.names = FALSE
)

cat("All markers saved:\n")
cat(all_markers_file, "\n\n")


# ------------------------------------------------------------
# 18. Select top markers per cluster
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("TOP MARKERS PER CLUSTER\n")
cat("------------------------------------------------------------\n\n")

markers_top <- do.call(
  rbind,
  lapply(
    split(markers_all, markers_all$cluster),
    function(df) {

      df <- df[
        order(
          df[[fc_column]],
          decreasing = TRUE
        ),
        ,
        drop = FALSE
      ]

      head(
        df,
        top_n_markers
      )

    }
  )
)

rownames(markers_top) <- NULL


# ------------------------------------------------------------
# 19. Save top markers
# ------------------------------------------------------------

top_markers_file <- file.path(
  output_dir,
  "top20_markers_per_cluster.csv"
)

write.csv(
  markers_top,
  top_markers_file,
  row.names = FALSE
)

cat("Top ", top_n_markers, " markers per cluster saved:\n", sep = "")
cat(top_markers_file, "\n\n")


# ------------------------------------------------------------
# 20. Print TOP 10 markers
# ------------------------------------------------------------

cat("============================================================\n")
cat("TOP 10 MARKERS FOR EACH CLUSTER\n")
cat("============================================================\n\n")

for (cluster_id in clusters) {

  cat("\n")
  cat("CLUSTER ", cluster_id, "\n", sep = "")
  cat("--------------------------------------------\n")

  cluster_markers <- markers_top[
    markers_top$cluster == cluster_id,
    ,
    drop = FALSE
  ]

  top10 <- head(
    cluster_markers,
    10
  )

  print(
    top10[
      ,
      c(
        "gene",
        fc_column,
        "pct.1",
        "pct.2",
        "p_val_adj"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )

}


# ------------------------------------------------------------
# 21. Heatmap genes
# ------------------------------------------------------------

cat("\n")
cat("------------------------------------------------------------\n")
cat("MARKER HEATMAP\n")
cat("------------------------------------------------------------\n\n")


heatmap_genes <- unique(
  unlist(
    lapply(
      split(markers_top, markers_top$cluster),
      function(df) {
        head(
          df$gene,
          heatmap_top_n
        )
      }
    )
  )
)


heatmap_genes <- heatmap_genes[
  heatmap_genes %in% rownames(seurat_obj)
]


cat(
  "Genes used in heatmap: ",
  length(heatmap_genes),
  "\n",
  sep = ""
)


if (length(heatmap_genes) > 0) {

  p_heatmap <- DoHeatmap(
    seurat_obj,
    features = heatmap_genes,
    group.by = "seurat_clusters"
  ) +
    ggtitle(
      "Top cluster marker genes"
    )

  ggsave(
    filename = file.path(
      output_dir,
      "cluster_marker_heatmap.png"
    ),
    plot = p_heatmap,
    width = 16,
    height = 14,
    dpi = 300
  )

  cat("Marker heatmap saved.\n\n")

} else {

  warning(
    "No genes available for marker heatmap."
  )

}


# ------------------------------------------------------------
# 22. Broad cell-type marker DotPlot
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("BROAD CELL-TYPE MARKERS\n")
cat("------------------------------------------------------------\n\n")


broad_markers <- c(

  # Astrocytes
  "GFAP",
  "AQP4",
  "SLC1A3",
  "ALDH1L1",

  # Oligodendrocytes
  "MBP",
  "MOG",
  "MOBP",
  "PLP1",

  # OPC
  "PDGFRA",
  "CSPG4",
  "VCAN",

  # Neurons
  "RBFOX3",
  "SYT1",
  "SNAP25",
  "STMN2",

  # Microglia / macrophages
  "C1QA",
  "C1QB",
  "C1QC",
  "TYROBP",
  "AIF1",
  "FCER1G",

  # Endothelial
  "CLDN5",
  "VWF",
  "PECAM1",
  "KDR",
  "ESAM",

  # Pericytes
  "RGS5",
  "PDGFRB",
  "MCAM",

  # Malignant / glioma-associated
  "EGFR",
  "SOX2",
  "OLIG2",
  "NES",
  "PDGFRA",

  # Proliferation
  "MKI67",
  "TOP2A",
  "UBE2C",
  "CENPF"
)


broad_markers_present <- unique(
  broad_markers[
    broad_markers %in% rownames(seurat_obj)
  ]
)


cat(
  "Markers found in dataset: ",
  length(broad_markers_present),
  " / ",
  length(unique(broad_markers)),
  "\n",
  sep = ""
)


if (length(broad_markers_present) > 0) {

  p_dotplot <- DotPlot(
    seurat_obj,
    assay = "RNA",
    features = broad_markers_present,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis() +
    ggtitle(
      "Broad cell-type marker expression"
    )

  ggsave(
    filename = file.path(
      output_dir,
      "broad_cell_type_markers_DotPlot.png"
    ),
    plot = p_dotplot,
    width = 16,
    height = 9,
    dpi = 300
  )

  cat("Broad marker DotPlot saved.\n\n")

}


# ------------------------------------------------------------
# 23. Save marker-enhanced object
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("SAVING MARKER OBJECT\n")
cat("------------------------------------------------------------\n\n")


output_file <- file.path(
  output_dir,
  "GBM_3samples_Seurat_markers.rds"
)


saveRDS(
  seurat_obj,
  output_file
)


cat("Saved Seurat object:\n")
cat(output_file, "\n\n")


# ------------------------------------------------------------
# 24. Final summary
# ------------------------------------------------------------

cat("============================================================\n")
cat("MARKER ANALYSIS COMPLETED SUCCESSFULLY\n")
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
  "Clusters: ",
  length(clusters),
  "\n",
  sep = ""
)

cat(
  "Total marker rows: ",
  nrow(markers_all),
  "\n",
  sep = ""
)

cat(
  "Marker threshold: min.pct = ",
  min_pct_use,
  ", logFC >= ",
  logfc_threshold_use,
  "\n\n",
  sep = ""
)

cat("Output directory:\n")
cat(output_dir, "\n\n")

cat("Output files:\n")
cat("  all_cluster_markers.csv\n")
cat("  top20_markers_per_cluster.csv\n")
cat("  cluster_marker_heatmap.png\n")
cat("  broad_cell_type_markers_DotPlot.png\n")
cat("  GBM_3samples_Seurat_markers.rds\n\n")

cat("Next step:\n")
cat("Biological annotation of clusters based on marker genes.\n\n")

cat("============================================================\n")
cat("DONE\n")
cat("============================================================\n")
