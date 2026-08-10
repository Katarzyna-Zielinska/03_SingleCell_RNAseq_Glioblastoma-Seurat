#!/usr/bin/env Rscript

# ============================================================
# SCRIPT 26
# FULL OBSERVATION VALIDATION FOR inferCNV
# ============================================================

suppressPackageStartupMessages({
  library(infercnv)
  library(Matrix)
})

cat("\n")
cat("============================================================\n")
cat("FULL OBSERVATION VALIDATION FOR inferCNV\n")
cat("============================================================\n\n")

cat("infercnv version: ",
    as.character(packageVersion("infercnv")), "\n\n")

# ============================================================
# PATHS
# ============================================================

project_dir <- "/home/katja/Cancer-Bioinformatics-Portfolio/03_SingleCell-RNAseq-Glioblastoma-Seurat"

infercnv_file <- file.path(
  project_dir,
  "results/infercnv/GBM_3samples_infercnv_full_final.rds"
)

metadata_file <- file.path(
  project_dir,
  "results/cnv_input/CNV_cell_metadata.csv"
)

output_dir <- file.path(
  project_dir,
  "results/infercnv/inspection_full"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Input inferCNV object:\n")
cat(infercnv_file, "\n\n")

cat("Input cell metadata:\n")
cat(metadata_file, "\n\n")

cat("Output directory:\n")
cat(output_dir, "\n\n")


# ============================================================
# CHECK INPUT FILES
# ============================================================

if (!file.exists(infercnv_file)) {
  stop("ERROR: inferCNV object not found:\n", infercnv_file)
}

if (!file.exists(metadata_file)) {
  stop("ERROR: metadata file not found:\n", metadata_file)
}


# ============================================================
# LOAD INFERCNV OBJECT
# ============================================================

cat("Loading FULL inferCNV object...\n")

infercnv_obj <- readRDS(infercnv_file)

cat("Object loaded successfully.\n\n")

if (!inherits(infercnv_obj, "infercnv")) {
  stop("ERROR: Loaded object is not an infercnv object.")
}


# ============================================================
# LOAD METADATA
# ============================================================

cat("Loading original cell metadata...\n")

meta <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Metadata rows:", nrow(meta), "\n")
cat("Metadata columns:", ncol(meta), "\n\n")

cat("Metadata column names:\n")
print(colnames(meta))
cat("\n")


# ============================================================
# IDENTIFY CELL ID COLUMN
# ============================================================

# The CSV generated from Seurat contains cell IDs
# in the first column, whose column name may be empty.

cell_id_col <- 1

meta_cell_ids <- as.character(meta[[cell_id_col]])

cat("Metadata cell ID column index:", cell_id_col, "\n")
cat("Metadata cell ID column name:",
    ifelse(
      is.na(colnames(meta)[cell_id_col]) ||
        colnames(meta)[cell_id_col] == "",
      "<empty>",
      colnames(meta)[cell_id_col]
    ),
    "\n\n")

cat("First metadata cell IDs:\n")
print(head(meta_cell_ids))
cat("\n")


# ============================================================
# BASIC INFERCNV DIMENSIONS
# ============================================================

expr_mat <- infercnv_obj@expr.data
count_mat <- infercnv_obj@count.data

infercnv_cells <- colnames(expr_mat)

if (is.null(infercnv_cells)) {
  stop("ERROR: inferCNV expression matrix has no cell identifiers.")
}

cat("============================================================\n")
cat("INFERCNV OBJECT DIMENSIONS\n")
cat("============================================================\n\n")

cat("Expression matrix:\n")
cat("Genes:", nrow(expr_mat), "\n")
cat("Cells:", ncol(expr_mat), "\n\n")

cat("Count matrix:\n")
cat("Genes:", nrow(count_mat), "\n")
cat("Cells:", ncol(count_mat), "\n\n")

cat("Number of cell identifiers:", length(infercnv_cells), "\n\n")


# ============================================================
# CELL IDENTIFIER VALIDATION
# ============================================================

cat("============================================================\n")
cat("CELL IDENTIFIER VALIDATION\n")
cat("============================================================\n\n")

matched_cells <- infercnv_cells %in% meta_cell_ids

cat(
  "Cells in inferCNV object:",
  length(infercnv_cells),
  "\n"
)

cat(
  "Cells matched to metadata:",
  sum(matched_cells),
  "\n"
)

cat(
  "Cells missing from metadata:",
  sum(!matched_cells),
  "\n\n"
)

if (any(!matched_cells)) {

  cat("First missing cell IDs:\n")
  print(head(infercnv_cells[!matched_cells], 20))
  cat("\n")

  stop("ERROR: Some inferCNV cells are missing from metadata.")
}

cat("Cell identifiers: PASSED\n")


# ============================================================
# CELL ORDER
# ============================================================

metadata_index <- match(infercnv_cells, meta_cell_ids)

if (!all(meta_cell_ids[metadata_index] == infercnv_cells)) {
  stop("ERROR: Cell order mismatch between inferCNV object and metadata.")
}

cat("Cell order: PASSED\n\n")


# ============================================================
# REORDER METADATA TO INFERCNV ORDER
# ============================================================

meta_aligned <- meta[metadata_index, , drop = FALSE]

rownames(meta_aligned) <- infercnv_cells

if (!all(rownames(meta_aligned) == colnames(expr_mat))) {
  stop("ERROR: Metadata alignment failed.")
}

cat("Metadata alignment: PASSED\n\n")


# ============================================================
# REQUIRED METADATA COLUMNS
# ============================================================

required_columns <- c(
  "sample_id",
  "cell_type"
)

missing_columns <- setdiff(
  required_columns,
  colnames(meta_aligned)
)

if (length(missing_columns) > 0) {
  stop(
    "ERROR: Missing required metadata columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

cat("Required metadata columns: PASSED\n\n")


# ============================================================
# IDENTIFY REFERENCE AND OBSERVATION CELLS
# ============================================================

reference_indices <- sort(
  unique(
    unlist(
      infercnv_obj@reference_grouped_cell_indices,
      use.names = FALSE
    )
  )
)

observation_indices <- sort(
  unique(
    unlist(
      infercnv_obj@observation_grouped_cell_indices,
      use.names = FALSE
    )
  )
)

reference_indices <- reference_indices[
  reference_indices >= 1 &
    reference_indices <= ncol(expr_mat)
]

observation_indices <- observation_indices[
  observation_indices >= 1 &
    observation_indices <= ncol(expr_mat)
]

reference_cells <- infercnv_cells[reference_indices]
observation_cells <- infercnv_cells[observation_indices]

cat("============================================================\n")
cat("INFERCNV CELL GROUPS\n")
cat("============================================================\n\n")

cat("Reference cells:", length(reference_cells), "\n")
cat("Observation cells:", length(observation_cells), "\n")
cat("Total cells:", length(infercnv_cells), "\n\n")

if (length(reference_cells) + length(observation_cells) !=
    length(infercnv_cells)) {

  warning(
    "Reference + observation cells do not equal total cells."
  )
}

if (length(intersect(reference_cells, observation_cells)) > 0) {
  stop("ERROR: Some cells are simultaneously reference and observation.")
}

cat("Reference/observation assignment: PASSED\n\n")


# ============================================================
# OBSERVATION METADATA
# ============================================================

obs_meta <- meta_aligned[observation_cells, , drop = FALSE]

ref_meta <- meta_aligned[reference_cells, , drop = FALSE]


# ============================================================
# SAMPLE REPRESENTATION
# ============================================================

cat("============================================================\n")
cat("OBSERVATION CELLS BY SAMPLE\n")
cat("============================================================\n\n")

obs_by_sample <- as.data.frame(
  table(obs_meta$sample_id),
  stringsAsFactors = FALSE
)

colnames(obs_by_sample) <- c(
  "sample_id",
  "observation_cells"
)

print(obs_by_sample)

cat("\n")

obs_sample_total <- as.data.frame(
  table(meta_aligned$sample_id),
  stringsAsFactors = FALSE
)

colnames(obs_sample_total) <- c(
  "sample_id",
  "total_cells"
)

sample_summary <- merge(
  obs_sample_total,
  obs_by_sample,
  by = "sample_id",
  all = TRUE
)

sample_summary$observation_percent <-
  100 *
  sample_summary$observation_cells /
  sample_summary$total_cells

print(sample_summary)

cat("\n")


# ============================================================
# CELL TYPE REPRESENTATION
# ============================================================

cat("============================================================\n")
cat("OBSERVATION CELLS BY CELL TYPE\n")
cat("============================================================\n\n")

obs_by_celltype <- as.data.frame(
  table(obs_meta$cell_type),
  stringsAsFactors = FALSE
)

colnames(obs_by_celltype) <- c(
  "cell_type",
  "observation_cells"
)

obs_by_celltype <- obs_by_celltype[
  order(-obs_by_celltype$observation_cells),
]

print(obs_by_celltype)

cat("\n")


# ============================================================
# SAMPLE x CELL TYPE
# ============================================================

cat("============================================================\n")
cat("OBSERVATION CELLS BY SAMPLE AND CELL TYPE\n")
cat("============================================================\n\n")

sample_celltype <- as.data.frame(
  table(
    sample_id = obs_meta$sample_id,
    cell_type = obs_meta$cell_type
  ),
  stringsAsFactors = FALSE
)

sample_celltype <- sample_celltype[
  sample_celltype$Freq > 0,
]

colnames(sample_celltype)[3] <- "cells"

print(sample_celltype)

cat("\n")


# ============================================================
# CNV SIGNAL
# ============================================================

cat("============================================================\n")
cat("OBSERVATION CNV SIGNAL\n")
cat("============================================================\n\n")

cat("Calculating per-cell mean signal...\n")

obs_expr <- expr_mat[, observation_indices, drop = FALSE]

per_cell_mean <- Matrix::colMeans(obs_expr)

per_cell_median <- apply(
  obs_expr,
  2,
  median
)

per_cell_sd <- apply(
  obs_expr,
  2,
  sd
)

per_cell_min <- apply(
  obs_expr,
  2,
  min
)

per_cell_max <- apply(
  obs_expr,
  2,
  max
)

observation_signal <- data.frame(
  cell_id = observation_cells,
  sample_id = obs_meta$sample_id,
  cell_type = obs_meta$cell_type,
  mean_signal = per_cell_mean,
  median_signal = per_cell_median,
  sd_signal = per_cell_sd,
  min_signal = per_cell_min,
  max_signal = per_cell_max,
  stringsAsFactors = FALSE
)

cat("Signal calculations completed.\n\n")


# ============================================================
# GLOBAL SIGNAL SUMMARY
# ============================================================

cat("============================================================\n")
cat("GLOBAL OBSERVATION CNV SIGNAL\n")
cat("============================================================\n\n")

global_signal <- data.frame(
  metric = c(
    "minimum",
    "first_quartile",
    "median",
    "mean",
    "third_quartile",
    "maximum"
  ),
  value = c(
    min(observation_signal$mean_signal),
    quantile(
      observation_signal$mean_signal,
      0.25
    ),
    median(observation_signal$mean_signal),
    mean(observation_signal$mean_signal),
    quantile(
      observation_signal$mean_signal,
      0.75
    ),
    max(observation_signal$mean_signal)
  )
)

print(global_signal)

cat("\n")


# ============================================================
# SIGNAL BY SAMPLE
# ============================================================

cat("============================================================\n")
cat("CNV SIGNAL BY SAMPLE\n")
cat("============================================================\n\n")

signal_by_sample <- aggregate(
  mean_signal ~ sample_id,
  data = observation_signal,
  FUN = function(x) {
    c(
      cells = length(x),
      mean = mean(x),
      sd = sd(x),
      median = median(x),
      min = min(x),
      max = max(x)
    )
  }
)

signal_by_sample <- do.call(
  data.frame,
  signal_by_sample
)

colnames(signal_by_sample) <- c(
  "sample_id",
  "cells",
  "mean_signal",
  "sd_signal",
  "median_signal",
  "min_signal",
  "max_signal"
)

print(signal_by_sample)

cat("\n")


# ============================================================
# SIGNAL BY CELL TYPE
# ============================================================

cat("============================================================\n")
cat("CNV SIGNAL BY CELL TYPE\n")
cat("============================================================\n\n")

signal_by_celltype <- aggregate(
  mean_signal ~ cell_type,
  data = observation_signal,
  FUN = function(x) {
    c(
      cells = length(x),
      mean = mean(x),
      sd = sd(x),
      median = median(x),
      min = min(x),
      max = max(x)
    )
  }
)

signal_by_celltype <- do.call(
  data.frame,
  signal_by_celltype
)

colnames(signal_by_celltype) <- c(
  "cell_type",
  "cells",
  "mean_signal",
  "sd_signal",
  "median_signal",
  "min_signal",
  "max_signal"
)

signal_by_celltype <- signal_by_celltype[
  order(-signal_by_celltype$cells),
]

print(signal_by_celltype)

cat("\n")


# ============================================================
# HIGH / LOW SIGNAL CELLS
# ============================================================

cat("============================================================\n")
cat("EXTREME SIGNAL CELLS\n")
cat("============================================================\n\n")

signal_q1 <- quantile(
  observation_signal$mean_signal,
  0.01
)

signal_q99 <- quantile(
  observation_signal$mean_signal,
  0.99
)

cat("1st percentile:", signal_q1, "\n")
cat("99th percentile:", signal_q99, "\n\n")

low_signal_cells <- observation_signal[
  observation_signal$mean_signal <= signal_q1,
]

high_signal_cells <- observation_signal[
  observation_signal$mean_signal >= signal_q99,
]

cat(
  "Cells in lowest 1%:",
  nrow(low_signal_cells),
  "\n"
)

cat(
  "Cells in highest 1%:",
  nrow(high_signal_cells),
  "\n\n"
)


# ============================================================
# HIGH SIGNAL CELLS BY SAMPLE
# ============================================================

high_signal_by_sample <- as.data.frame(
  table(high_signal_cells$sample_id),
  stringsAsFactors = FALSE
)

colnames(high_signal_by_sample) <- c(
  "sample_id",
  "high_signal_cells"
)

cat("Highest-signal cells by sample:\n")
print(high_signal_by_sample)

cat("\n")


# ============================================================
# HIGH SIGNAL CELLS BY CELL TYPE
# ============================================================

high_signal_by_celltype <- as.data.frame(
  table(high_signal_cells$cell_type),
  stringsAsFactors = FALSE
)

colnames(high_signal_by_celltype) <- c(
  "cell_type",
  "high_signal_cells"
)

high_signal_by_celltype <- high_signal_by_celltype[
  high_signal_by_celltype$high_signal_cells > 0,
]

cat("Highest-signal cells by cell type:\n")
print(high_signal_by_celltype)

cat("\n")


# ============================================================
# CHECK FOR NUMERIC PROBLEMS
# ============================================================

cat("============================================================\n")
cat("NUMERIC VALIDATION\n")
cat("============================================================\n\n")

na_count <- sum(
  is.na(observation_signal$mean_signal)
)

inf_count <- sum(
  is.infinite(observation_signal$mean_signal)
)

cat("NA mean-signal values:", na_count, "\n")
cat("Infinite mean-signal values:", inf_count, "\n")

if (na_count > 0 || inf_count > 0) {
  stop("ERROR: Invalid numeric CNV signal values detected.")
}

cat("Numeric signal validation: PASSED\n\n")


# ============================================================
# SAVE RESULTS
# ============================================================

write.csv(
  observation_signal,
  file.path(
    output_dir,
    "infercnv_full_observation_per_cell_signal.csv"
  ),
  row.names = FALSE
)

saveRDS(
  observation_signal,
  file.path(
    output_dir,
    "infercnv_full_observation_per_cell_signal.rds"
  )
)

write.csv(
  sample_summary,
  file.path(
    output_dir,
    "infercnv_full_observation_sample_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  obs_by_celltype,
  file.path(
    output_dir,
    "infercnv_full_observation_celltype_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  sample_celltype,
  file.path(
    output_dir,
    "infercnv_full_observation_sample_x_celltype.csv"
  ),
  row.names = FALSE
)

write.csv(
  signal_by_sample,
  file.path(
    output_dir,
    "infercnv_full_observation_signal_by_sample.csv"
  ),
  row.names = FALSE
)

write.csv(
  signal_by_celltype,
  file.path(
    output_dir,
    "infercnv_full_observation_signal_by_celltype.csv"
  ),
  row.names = FALSE
)

write.csv(
  global_signal,
  file.path(
    output_dir,
    "infercnv_full_observation_global_signal_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  high_signal_cells,
  file.path(
    output_dir,
    "infercnv_full_observation_high_signal_cells.csv"
  ),
  row.names = FALSE
)

write.csv(
  low_signal_cells,
  file.path(
    output_dir,
    "infercnv_full_observation_low_signal_cells.csv"
  ),
  row.names = FALSE
)

write.csv(
  high_signal_by_sample,
  file.path(
    output_dir,
    "infercnv_full_observation_high_signal_by_sample.csv"
  ),
  row.names = FALSE
)

write.csv(
  high_signal_by_celltype,
  file.path(
    output_dir,
    "infercnv_full_observation_high_signal_by_celltype.csv"
  ),
  row.names = FALSE
)

write.csv(
  meta_aligned,
  file.path(
    output_dir,
    "infercnv_full_observation_metadata_aligned.csv"
  ),
  row.names = TRUE
)


# ============================================================
# SUMMARY TEXT
# ============================================================

summary_file <- file.path(
  output_dir,
  "infercnv_full_observation_validation_summary.txt"
)

sink(summary_file)

cat("FULL OBSERVATION VALIDATION FOR inferCNV\n")
cat("========================================\n\n")

cat(
  "infercnv version: ",
  as.character(packageVersion("infercnv")),
  "\n\n"
)

cat("Input object:\n")
cat(infercnv_file, "\n\n")

cat("Total cells:", length(infercnv_cells), "\n")
cat("Reference cells:", length(reference_cells), "\n")
cat("Observation cells:", length(observation_cells), "\n\n")

cat("Genes in expression matrix:", nrow(expr_mat), "\n")
cat("Genes in count matrix:", nrow(count_mat), "\n\n")

cat("Observation cells by sample:\n")
print(obs_by_sample)
cat("\n")

cat("Observation cells by cell type:\n")
print(obs_by_celltype)
cat("\n")

cat("Sample x cell type:\n")
print(sample_celltype)
cat("\n")

cat("Global observation CNV signal:\n")
print(global_signal)
cat("\n")

cat("CNV signal by sample:\n")
print(signal_by_sample)
cat("\n")

cat("CNV signal by cell type:\n")
print(signal_by_celltype)
cat("\n")

cat("1st percentile of per-cell mean signal:",
    signal_q1, "\n")

cat("99th percentile of per-cell mean signal:",
    signal_q99, "\n\n")

cat(
  "Cells in lowest 1%:",
  nrow(low_signal_cells),
  "\n"
)

cat(
  "Cells in highest 1%:",
  nrow(high_signal_cells),
  "\n\n"
)

cat("Numeric validation: PASSED\n")
cat("Cell alignment: PASSED\n")
cat("Metadata alignment: PASSED\n")
cat("Reference/observation assignment: PASSED\n\n")

cat("IMPORTANT:\n")
cat("This script validates the observation population and\n")
cat("quantifies CNV signal. It does NOT classify cells as\n")
cat("malignant or non-malignant.\n\n")

cat("END OF SCRIPT 26\n")

sink()


# ============================================================
# FINAL MESSAGE
# ============================================================

cat("============================================================\n")
cat("SCRIPT 26 COMPLETED\n")
cat("============================================================\n\n")

cat("Full observation validation completed successfully.\n\n")

cat("Reference cells:", length(reference_cells), "\n")
cat("Observation cells:", length(observation_cells), "\n")
cat("Total cells:", length(infercnv_cells), "\n\n")

cat("Output directory:\n")
cat(output_dir, "\n\n")

cat("Saved files include:\n")
cat("- infercnv_full_observation_per_cell_signal.csv\n")
cat("- infercnv_full_observation_per_cell_signal.rds\n")
cat("- infercnv_full_observation_sample_summary.csv\n")
cat("- infercnv_full_observation_celltype_summary.csv\n")
cat("- infercnv_full_observation_sample_x_celltype.csv\n")
cat("- infercnv_full_observation_signal_by_sample.csv\n")
cat("- infercnv_full_observation_signal_by_celltype.csv\n")
cat("- infercnv_full_observation_global_signal_summary.csv\n")
cat("- infercnv_full_observation_high_signal_cells.csv\n")
cat("- infercnv_full_observation_high_signal_by_celltype.csv\n")
cat("- infercnv_full_observation_validation_summary.txt\n\n")

cat("IMPORTANT:\n")
cat("Do NOT classify cells as malignant based on this script alone.\n\n")

cat("Script 26 is the final validation step of the main\n")
cat("inferCNV pipeline.\n")
