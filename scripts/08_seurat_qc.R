#!/usr/bin/env Rscript

# ============================================================
# Project: Single-cell RNA-seq Analysis of Human Glioblastoma
# Dataset: GSE139448
# Step: 08 - Seurat import and quality control
# Author: Katarzyna Zielinska
#
# Description:
# Import STARsolo Gene/filtered matrices for three samples,
# create Seurat objects, calculate QC metrics and generate
# diagnostic plots.
#
# IMPORTANT:
# No cell filtering is performed at this stage.
# Filtering thresholds will be selected after inspecting QC.
#
# Input format from STARsolo:
#   matrix.mtx
#   barcodes.tsv
#   features.tsv
#
# Seurat version: 5.x
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(Matrix)
    library(ggplot2)
})

cat("==========================================\n")
cat("Seurat QC - Glioblastoma scRNA-seq\n")
cat("==========================================\n\n")

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

project_dir <- path.expand(
    "~/Cancer-Bioinformatics-Portfolio/03_SingleCell-RNAseq-Glioblastoma-Seurat"
)

input_dir <- file.path(
    project_dir,
    "results",
    "starsolo"
)

output_dir <- file.path(
    project_dir,
    "results",
    "seurat_qc"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

cat("Project directory:\n")
cat(project_dir, "\n\n")

cat("STARsolo input directory:\n")
cat(input_dir, "\n\n")

cat("Seurat QC output directory:\n")
cat(output_dir, "\n\n")


# ------------------------------------------------------------
# 3. Sample IDs
# ------------------------------------------------------------

samples <- c(
    "SRR10353960",
    "SRR10353961",
    "SRR10353962"
)

cat("Samples:\n")
print(samples)
cat("\n")


# ------------------------------------------------------------
# 4. Create list for Seurat objects
# ------------------------------------------------------------

seurat_objects <- list()


# ------------------------------------------------------------
# 5. Import each STARsolo matrix
# ------------------------------------------------------------

for (sample_id in samples) {

    cat("\n")
    cat("==========================================\n")
    cat("Processing: ", sample_id, "\n", sep = "")
    cat("==========================================\n")

    # STARsolo output directory
    matrix_dir <- file.path(
        input_dir,
        paste0(sample_id, "_Solo.out"),
        "Gene",
        "filtered"
    )

    cat("\nMatrix directory:\n")
    cat(matrix_dir, "\n")

    # --------------------------------------------------------
    # Check directory
    # --------------------------------------------------------

    if (!dir.exists(matrix_dir)) {
        stop(
            "\nERROR: Matrix directory does not exist:\n",
            matrix_dir
        )
    }

    # --------------------------------------------------------
    # Define input files
    # --------------------------------------------------------

    matrix_file <- file.path(
        matrix_dir,
        "matrix.mtx"
    )

    barcode_file <- file.path(
        matrix_dir,
        "barcodes.tsv"
    )

    feature_file <- file.path(
        matrix_dir,
        "features.tsv"
    )

    # --------------------------------------------------------
    # Check input files
    # --------------------------------------------------------

    required_files <- c(
        matrix_file,
        barcode_file,
        feature_file
    )

    missing_files <- required_files[
        !file.exists(required_files)
    ]

    if (length(missing_files) > 0) {

        stop(
            "\nERROR: Missing STARsolo files for ",
            sample_id,
            ":\n",
            paste(
                missing_files,
                collapse = "\n"
            )
        )
    }

    cat("\nInput files found:\n")
    cat("  matrix.mtx   :", file.exists(matrix_file), "\n")
    cat("  barcodes.tsv :", file.exists(barcode_file), "\n")
    cat("  features.tsv :", file.exists(feature_file), "\n")


    # --------------------------------------------------------
    # Read STARsolo matrix
    # --------------------------------------------------------

    cat("\nReading matrix...\n")

    counts <- ReadMtx(
        mtx = matrix_file,
        cells = barcode_file,
        features = feature_file,
        feature.column = 2,
        cell.column = 1,
        unique.features = TRUE
    )

    cat("Matrix successfully loaded.\n\n")


    # --------------------------------------------------------
    # Matrix dimensions
    # --------------------------------------------------------

    cat("Matrix dimensions:\n")
    cat("  Genes :", nrow(counts), "\n")
    cat("  Cells :", ncol(counts), "\n")


    # --------------------------------------------------------
    # Create Seurat object
    # --------------------------------------------------------

    cat("\nCreating Seurat object...\n")

    seu <- CreateSeuratObject(
        counts = counts,
        project = sample_id,
        min.cells = 0,
        min.features = 0
    )


    # --------------------------------------------------------
    # Add sample metadata
    # --------------------------------------------------------

    seu$sample_id <- sample_id


    # --------------------------------------------------------
    # Calculate mitochondrial percentage
    # --------------------------------------------------------

    cat("Calculating mitochondrial percentage...\n")

    seu[["percent.mt"]] <- PercentageFeatureSet(
        seu,
        pattern = "^MT-"
    )


    # --------------------------------------------------------
    # Store object
    # --------------------------------------------------------

    seurat_objects[[sample_id]] <- seu

    cat("\nSeurat object created:\n")
    cat("  Cells :", ncol(seu), "\n")
    cat("  Genes :", nrow(seu), "\n")

    cat("\n")
}


# ------------------------------------------------------------
# 6. Save individual Seurat objects
# ------------------------------------------------------------

cat("\n")
cat("==========================================\n")
cat("Saving individual Seurat objects\n")
cat("==========================================\n")

for (sample_id in samples) {

    output_file <- file.path(
        output_dir,
        paste0(
            sample_id,
            "_seurat_raw.rds"
        )
    )

    saveRDS(
        seurat_objects[[sample_id]],
        file = output_file
    )

    cat(
        "Saved: ",
        output_file,
        "\n",
        sep = ""
    )
}


# ------------------------------------------------------------
# 7. Merge all three samples
# ------------------------------------------------------------

cat("\n")
cat("==========================================\n")
cat("Merging samples\n")
cat("==========================================\n")

seu_merged <- merge(
    x = seurat_objects[[1]],
    y = seurat_objects[2:3],
    add.cell.ids = samples,
    project = "GBM_scRNAseq"
)

cat("\nMerged object:\n")
cat("  Genes :", nrow(seu_merged), "\n")
cat("  Cells :", ncol(seu_merged), "\n")


# ------------------------------------------------------------
# 8. Cells per sample
# ------------------------------------------------------------

cat("\n")
cat("==========================================\n")
cat("Cells per sample\n")
cat("==========================================\n")

sample_counts <- table(
    seu_merged$sample_id
)

print(sample_counts)

write.csv(
    as.data.frame(sample_counts),
    file = file.path(
        output_dir,
        "cells_per_sample.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# 9. Create QC table
# ------------------------------------------------------------

cat("\n")
cat("==========================================\n")
cat("Creating QC metrics table\n")
cat("==========================================\n")

qc_metadata <- seu_merged[[]]

qc_summary <- data.frame(
    cell = rownames(qc_metadata),
    sample_id = qc_metadata$sample_id,
    nCount_RNA = qc_metadata$nCount_RNA,
    nFeature_RNA = qc_metadata$nFeature_RNA,
    percent_mt = qc_metadata$percent.mt
)

write.csv(
    qc_summary,
    file = file.path(
        output_dir,
        "QC_metrics_per_cell.csv"
    ),
    row.names = FALSE
)

cat("QC table saved.\n")


# ------------------------------------------------------------
# 10. QC summary statistics
# ------------------------------------------------------------

cat("\n")
cat("==========================================\n")
cat("QC summary statistics\n")
cat("==========================================\n")

qc_stats <- summary(
    qc_summary[
        ,
        c(
            "nCount_RNA",
            "nFeature_RNA",
            "percent_mt"
        )
    ]
)

print(qc_stats)

write.csv(
    as.data.frame(
        qc_stats
    ),
    file = file.path(
        output_dir,
        "QC_summary_statistics.csv"
    )
)


# ------------------------------------------------------------
# 11. Per-sample QC statistics
# ------------------------------------------------------------

cat("\n")
cat("==========================================\n")
cat("Per-sample QC statistics\n")
cat("==========================================\n")

for (sample_id in samples) {

    cat("\n------------------------------------------\n")
    cat(sample_id, "\n")
    cat("------------------------------------------\n")

    sample_qc <- qc_summary[
        qc_summary$sample_id == sample_id,
        ,
        drop = FALSE
    ]

    print(
        summary(
            sample_qc[
                ,
                c(
                    "nCount_RNA",
                    "nFeature_RNA",
                    "percent_mt"
                )
            ]
        )
    )
}


# ------------------------------------------------------------
# 12. QC plots
# ------------------------------------------------------------

cat("\n")
cat("==========================================\n")
cat("Generating QC plots\n")
cat("==========================================\n")


# ------------------------------------------------------------
# 12.1 nFeature_RNA
# ------------------------------------------------------------

p_nfeature <- VlnPlot(
    seu_merged,
    features = "nFeature_RNA",
    group.by = "sample_id",
    pt.size = 0.1
) +
    ggtitle(
        "Number of detected genes per cell"
    ) +
    theme(
        plot.title = element_text(
            hjust = 0.5
        )
    )

ggsave(
    filename = file.path(
        output_dir,
        "QC_nFeature_RNA_by_sample.png"
    ),
    plot = p_nfeature,
    width = 9,
    height = 6,
    dpi = 300
)


# ------------------------------------------------------------
# 12.2 nCount_RNA
# ------------------------------------------------------------

p_ncount <- VlnPlot(
    seu_merged,
    features = "nCount_RNA",
    group.by = "sample_id",
    pt.size = 0.1
) +
    ggtitle(
        "RNA counts per cell"
    ) +
    theme(
        plot.title = element_text(
            hjust = 0.5
        )
    )

ggsave(
    filename = file.path(
        output_dir,
        "QC_nCount_RNA_by_sample.png"
    ),
    plot = p_ncount,
    width = 9,
    height = 6,
    dpi = 300
)


# ------------------------------------------------------------
# 12.3 Mitochondrial percentage
# ------------------------------------------------------------

p_mt <- VlnPlot(
    seu_merged,
    features = "percent.mt",
    group.by = "sample_id",
    pt.size = 0.1
) +
    ggtitle(
        "Mitochondrial RNA percentage"
    ) +
    theme(
        plot.title = element_text(
            hjust = 0.5
        )
    )

ggsave(
    filename = file.path(
        output_dir,
        "QC_percent_mt_by_sample.png"
    ),
    plot = p_mt,
    width = 9,
    height = 6,
    dpi = 300
)


# ------------------------------------------------------------
# 12.4 nCount vs nFeature
# ------------------------------------------------------------

p_count_feature <- FeatureScatter(
    seu_merged,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
)

ggsave(
    filename = file.path(
        output_dir,
        "QC_nCount_vs_nFeature.png"
    ),
    plot = p_count_feature,
    width = 8,
    height = 6,
    dpi = 300
)


# ------------------------------------------------------------
# 12.5 nCount vs mitochondrial percentage
# ------------------------------------------------------------

p_count_mt <- FeatureScatter(
    seu_merged,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
)

ggsave(
    filename = file.path(
        output_dir,
        "QC_nCount_vs_percent_mt.png"
    ),
    plot = p_count_mt,
    width = 8,
    height = 6,
    dpi = 300
)


# ------------------------------------------------------------
# 13. Save merged raw Seurat object
# ------------------------------------------------------------

cat("\n")
cat("==========================================\n")
cat("Saving merged Seurat object\n")
cat("==========================================\n")

merged_file <- file.path(
    output_dir,
    "GBM_3samples_Seurat_raw.rds"
)

saveRDS(
    seu_merged,
    file = merged_file
)

cat(
    "Saved: ",
    merged_file,
    "\n",
    sep = ""
)


# ------------------------------------------------------------
# 14. Final summary
# ------------------------------------------------------------

cat("\n")
cat("==========================================\n")
cat("QC IMPORT COMPLETED SUCCESSFULLY\n")
cat("==========================================\n")

cat(
    "Number of samples: ",
    length(samples),
    "\n",
    sep = ""
)

cat(
    "Total genes: ",
    nrow(seu_merged),
    "\n",
    sep = ""
)

cat(
    "Total cells: ",
    ncol(seu_merged),
    "\n",
    sep = ""
)

cat("\n")
cat("No cells were filtered.\n")

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat("\n")
cat("Next step:\n")
cat(
    "Inspect QC distributions and choose ",
    "data-driven filtering thresholds.\n",
    sep = ""
)

cat("\n==========================================\n")
cat("DONE\n")
cat("==========================================\n")
