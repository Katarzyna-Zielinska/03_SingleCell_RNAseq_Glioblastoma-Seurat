#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

cat("\n")
cat("============================================================\n")
cat("CNV REFERENCE POPULATION VALIDATION\n")
cat("============================================================\n\n")

cat("Seurat version: ", as.character(packageVersion("Seurat")), "\n")
cat("SeuratObject version: ", as.character(packageVersion("SeuratObject")), "\n\n")

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

seurat_file <- "results/seurat_annotation/GBM_3samples_Seurat_annotated.rds"

counts_file <- "results/cnv_input/CNV_raw_counts_GENCODE_v47.rds"

metadata_file <- "results/cnv_input/CNV_cell_metadata.csv"

output_dir <- "results/cnv_input"

# ------------------------------------------------------------
# LOAD SEURAT
# ------------------------------------------------------------

cat("Loading Seurat object...\n")

obj <- readRDS(seurat_file)

cat("Object loaded successfully.\n\n")

cat("Genes: ", nrow(obj), "\n")
cat("Cells: ", ncol(obj), "\n\n")

# ------------------------------------------------------------
# LOAD COUNTS
# ------------------------------------------------------------

cat("Loading raw counts...\n")

counts <- readRDS(counts_file)

cat("Raw counts loaded.\n")

cat("Genes: ", nrow(counts), "\n")
cat("Cells: ", ncol(counts), "\n\n")

# ------------------------------------------------------------
# LOAD METADATA
# ------------------------------------------------------------

cat("Loading CNV metadata...\n")

meta <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Metadata loaded.\n\n")

cat("Metadata rows: ", nrow(meta), "\n")
cat("Metadata columns: ", ncol(meta), "\n\n")

# ------------------------------------------------------------
# IDENTIFY CELL ID COLUMN
# ------------------------------------------------------------

# The CSV was created from Seurat metadata.
# Therefore the cell IDs are stored in the first column,
# which may have an empty column name after read.csv().

cell_id_column <- 1

cell_ids <- meta[[cell_id_column]]

if (is.null(cell_ids)) {
  stop("ERROR: Could not identify cell IDs in metadata.")
}

cell_ids <- as.character(cell_ids)

cat("Metadata cell ID column: first column\n")
cat("First metadata cell IDs:\n")
print(head(cell_ids))

cat("\n")

# Rename first column to "cell"
colnames(meta)[cell_id_column] <- "cell"

# ------------------------------------------------------------
# BASIC CELL ID CHECKS
# ------------------------------------------------------------

count_cells <- colnames(counts)

if (is.null(count_cells)) {
  stop("ERROR: Counts matrix has no column names.")
}

count_cells <- as.character(count_cells)

cat("Cells in counts:   ", length(count_cells), "\n")
cat("Cells in metadata: ", length(cell_ids), "\n\n")

# Check duplicates

if (anyDuplicated(count_cells) > 0) {
  stop("ERROR: Duplicate cell IDs found in counts matrix.")
}

if (anyDuplicated(cell_ids) > 0) {
  stop("ERROR: Duplicate cell IDs found in metadata.")
}

# Check missing / extra cells

missing_in_metadata <- setdiff(count_cells, cell_ids)

extra_in_metadata <- setdiff(cell_ids, count_cells)

if (length(missing_in_metadata) > 0) {

  cat("Cells present in counts but missing from metadata:\n")
  print(head(missing_in_metadata, 20))

  stop("ERROR: Counts contains cells missing from metadata.")
}

if (length(extra_in_metadata) > 0) {

  cat("Cells present in metadata but absent from counts:\n")
  print(head(extra_in_metadata, 20))

  stop("ERROR: Metadata contains cells absent from counts.")
}

cat("Cell identifiers: PASSED\n")

# ------------------------------------------------------------
# REORDER METADATA
# ------------------------------------------------------------

meta <- meta[
  match(count_cells, meta$cell),
  ,
  drop = FALSE
]

rownames(meta) <- meta$cell

# Final alignment check

if (!identical(colnames(counts), rownames(meta))) {
  stop("ERROR: Cell alignment failed after reordering.")
}

cat("Cell order: PASSED\n\n")

# ------------------------------------------------------------
# REQUIRED METADATA COLUMNS
# ------------------------------------------------------------

required_columns <- c(
  "sample_id",
  "cell_type"
)

missing_columns <- setdiff(
  required_columns,
  colnames(meta)
)

if (length(missing_columns) > 0) {

  cat("Missing metadata columns:\n")
  print(missing_columns)

  stop("ERROR: Required metadata columns are missing.")
}

cat("Required metadata columns: PASSED\n\n")

# ------------------------------------------------------------
# SAMPLE SUMMARY
# ------------------------------------------------------------

cat("============================================================\n")
cat("SAMPLE SUMMARY\n")
cat("============================================================\n\n")

sample_counts <- table(meta$sample_id)

print(sample_counts)

# ------------------------------------------------------------
# CELL TYPE SUMMARY
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CELL TYPE SUMMARY\n")
cat("============================================================\n\n")

celltype_counts <- table(meta$cell_type)

print(celltype_counts)

# ------------------------------------------------------------
# SAMPLE x CELL TYPE
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SAMPLE x CELL TYPE\n")
cat("============================================================\n\n")

sample_celltype <- table(
  meta$sample_id,
  meta$cell_type
)

print(sample_celltype)

# ------------------------------------------------------------
# REFERENCE POPULATIONS
# ------------------------------------------------------------

reference_types <- c(
  "Astrocyte-like",
  "Oligodendrocyte",
  "Microglia-Macrophage",
  "Endothelial",
  "Mast-cell"
)

cat("\n============================================================\n")
cat("POTENTIAL CNV REFERENCE POPULATIONS\n")
cat("============================================================\n\n")

reference_present <- reference_types[
  reference_types %in% unique(meta$cell_type)
]

reference_missing <- reference_types[
  !reference_types %in% unique(meta$cell_type)
]

cat("Reference populations present:\n\n")

for (ct in reference_present) {

  n <- sum(meta$cell_type == ct)

  cat(
    sprintf(
      "%-25s %6d cells\n",
      ct,
      n
    )
  )
}

if (length(reference_missing) > 0) {

  cat("\nReference populations missing:\n\n")

  for (ct in reference_missing) {

    cat(
      sprintf(
        "%-25s MISSING\n",
        ct
      )
    )
  }
}

# ------------------------------------------------------------
# REFERENCE BY SAMPLE
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("REFERENCE POPULATIONS BY SAMPLE\n")
cat("============================================================\n\n")

reference_meta <- meta[
  meta$cell_type %in% reference_present,
  ,
  drop = FALSE
]

reference_by_sample <- table(
  reference_meta$sample_id,
  reference_meta$cell_type
)

print(reference_by_sample)

# ------------------------------------------------------------
# REFERENCE PERCENTAGE
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("REFERENCE PERCENTAGE WITHIN SAMPLE\n")
cat("============================================================\n\n")

sample_totals <- table(meta$sample_id)

reference_percent <- sweep(
  reference_by_sample,
  1,
  sample_totals[rownames(reference_by_sample)],
  FUN = "/"
) * 100

reference_percent <- round(
  reference_percent,
  3
)

print(reference_percent)

# ------------------------------------------------------------
# REFERENCE CELL IDs
# ------------------------------------------------------------

reference_cells <- meta$cell[
  meta$cell_type %in% reference_present
]

cat("\nTotal reference cells: ", length(reference_cells), "\n")

# ------------------------------------------------------------
# CHECK REFERENCE COVERAGE
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("REFERENCE COVERAGE\n")
cat("============================================================\n\n")

for (sample in names(sample_counts)) {

  sample_total <- sum(
    meta$sample_id == sample
  )

  sample_reference <- sum(
    meta$sample_id == sample &
      meta$cell_type %in% reference_present
  )

  percentage <- 100 *
    sample_reference /
    sample_total

  cat(
    sprintf(
      "%s: %d / %d cells (%.2f%% reference)\n",
      sample,
      sample_reference,
      sample_total,
      percentage
    )
  )
}

# ------------------------------------------------------------
# REFERENCE QUALITY FLAGS
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("REFERENCE QUALITY FLAGS\n")
cat("============================================================\n\n")

for (sample in names(sample_counts)) {

  sample_reference <- sum(
    meta$sample_id == sample &
      meta$cell_type %in% reference_present
  )

  if (sample_reference == 0) {

    cat(
      sample,
      ": WARNING - no reference cells\n"
    )

  } else if (sample_reference < 50) {

    cat(
      sample,
      ": WARNING - fewer than 50 reference cells\n"
    )

  } else {

    cat(
      sample,
      ": reference population available\n"
    )
  }
}

# ------------------------------------------------------------
# SAVE REFERENCE CELL METADATA
# ------------------------------------------------------------

reference_metadata <- meta[
  meta$cell_type %in% reference_present,
  ,
  drop = FALSE
]

write.csv(
  reference_metadata,
  file = file.path(
    output_dir,
    "CNV_reference_cells_metadata.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# SAVE REFERENCE CELL IDS
# ------------------------------------------------------------

write.table(
  reference_cells,
  file = file.path(
    output_dir,
    "CNV_reference_cell_ids.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# ------------------------------------------------------------
# SAVE REFERENCE BY SAMPLE
# ------------------------------------------------------------

write.csv(
  as.data.frame.matrix(reference_by_sample),
  file = file.path(
    output_dir,
    "CNV_reference_by_sample.csv"
  ),
  row.names = TRUE
)

# ------------------------------------------------------------
# SAVE REFERENCE PERCENTAGES
# ------------------------------------------------------------

write.csv(
  as.data.frame.matrix(reference_percent),
  file = file.path(
    output_dir,
    "CNV_reference_percent_by_sample.csv"
  ),
  row.names = TRUE
)

# ------------------------------------------------------------
# SUMMARY FILE
# ------------------------------------------------------------

summary_file <- file.path(
  output_dir,
  "CNV_reference_validation_summary.txt"
)

sink(summary_file)

cat("CNV REFERENCE POPULATION VALIDATION\n")
cat("==================================\n\n")

cat("Total genes:", nrow(counts), "\n")
cat("Total cells:", ncol(counts), "\n\n")

cat("Cells per sample:\n")
print(sample_counts)

cat("\nCells per cell type:\n")
print(celltype_counts)

cat("\nReference populations:\n")
print(reference_present)

cat("\nReference cells per sample:\n")
print(reference_by_sample)

cat("\nReference percentage per sample:\n")
print(reference_percent)

sink()

# ------------------------------------------------------------
# FINAL CHECK
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("REFERENCE VALIDATION COMPLETED\n")
cat("============================================================\n\n")

cat("Cell alignment: PASSED\n")
cat("Metadata alignment: PASSED\n")
cat("Reference populations identified: ", length(reference_present), "\n")
cat("Total reference cells: ", length(reference_cells), "\n\n")

cat("Saved files:\n")
cat(
  file.path(
    output_dir,
    "CNV_reference_cells_metadata.csv"
  ),
  "\n"
)

cat(
  file.path(
    output_dir,
    "CNV_reference_cell_ids.txt"
  ),
  "\n"
)

cat(
  file.path(
    output_dir,
    "CNV_reference_by_sample.csv"
  ),
  "\n"
)

cat(
  file.path(
    output_dir,
    "CNV_reference_percent_by_sample.csv"
  ),
  "\n"
)

cat(
  summary_file,
  "\n"
)

cat("\nNext step:\n")
cat("Evaluate whether the reference populations are suitable for\n")
cat("sample-aware CNV inference before running the CNV algorithm.\n\n")
