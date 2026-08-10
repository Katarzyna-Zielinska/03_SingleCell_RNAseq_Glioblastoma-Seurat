#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

cat("\n")
cat("============================================================\n")
cat("CNV REFERENCE SELECTION\n")
cat("============================================================\n\n")

cat("Seurat version: ", as.character(packageVersion("Seurat")), "\n")
cat("SeuratObject version: ", as.character(packageVersion("SeuratObject")), "\n\n")

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

counts_file <- "results/cnv_input/CNV_raw_counts_GENCODE_v47.rds"

metadata_file <- "results/cnv_input/CNV_cell_metadata.csv"

output_dir <- "results/cnv_input"

# ------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------

# Conservative core reference populations.
# These populations were identified as potential
# non-malignant populations during previous validation.

reference_types <- c(
  "Microglia-Macrophage",
  "Oligodendrocyte"
)

# Maximum number of reference cells per sample.
# This prevents one sample from dominating the reference.

max_cells_per_sample <- 100

# Reproducible random selection
set.seed(1234)

# ------------------------------------------------------------
# LOAD COUNTS
# ------------------------------------------------------------

cat("Loading CNV counts...\n")

counts <- readRDS(counts_file)

cat("Counts loaded.\n")
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

# First column contains cell IDs
colnames(meta)[1] <- "cell"

meta$cell <- as.character(meta$cell)

cat("Metadata loaded.\n")
cat("Rows: ", nrow(meta), "\n\n")

# ------------------------------------------------------------
# ALIGN METADATA TO COUNTS
# ------------------------------------------------------------

count_cells <- colnames(counts)

if (is.null(count_cells)) {
  stop("ERROR: Counts matrix has no cell names.")
}

count_cells <- as.character(count_cells)

if (anyDuplicated(count_cells) > 0) {
  stop("ERROR: Duplicate cell IDs in counts.")
}

if (anyDuplicated(meta$cell) > 0) {
  stop("ERROR: Duplicate cell IDs in metadata.")
}

missing_cells <- setdiff(
  count_cells,
  meta$cell
)

extra_cells <- setdiff(
  meta$cell,
  count_cells
)

if (length(missing_cells) > 0) {
  stop("ERROR: Cells missing from metadata.")
}

if (length(extra_cells) > 0) {
  stop("ERROR: Extra cells in metadata.")
}

meta <- meta[
  match(count_cells, meta$cell),
  ,
  drop = FALSE
]

rownames(meta) <- meta$cell

if (!identical(
  colnames(counts),
  rownames(meta)
)) {
  stop("ERROR: Counts and metadata are not aligned.")
}

cat("Cell alignment: PASSED\n\n")

# ------------------------------------------------------------
# CHECK REQUIRED COLUMNS
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
  stop(
    paste(
      "ERROR: Missing metadata columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# AVAILABLE CELL TYPES
# ------------------------------------------------------------

cat("Available cell types:\n\n")
print(table(meta$cell_type))

cat("\n")

# ------------------------------------------------------------
# CHECK REFERENCE TYPES
# ------------------------------------------------------------

available_reference_types <- reference_types[
  reference_types %in% unique(meta$cell_type)
]

missing_reference_types <- reference_types[
  !reference_types %in% unique(meta$cell_type)
]

if (length(missing_reference_types) > 0) {

  cat("WARNING: Missing reference populations:\n")
  print(missing_reference_types)
  cat("\n")
}

if (length(available_reference_types) == 0) {
  stop("ERROR: No reference populations available.")
}

cat("Reference populations selected:\n")
print(available_reference_types)
cat("\n")

# ------------------------------------------------------------
# INITIAL REFERENCE POPULATION
# ------------------------------------------------------------

reference_meta <- meta[
  meta$cell_type %in% available_reference_types,
  ,
  drop = FALSE
]

cat("Initial reference cells: ",
    nrow(reference_meta),
    "\n\n")

# ------------------------------------------------------------
# REFERENCE COUNTS BY SAMPLE
# ------------------------------------------------------------

cat("============================================================\n")
cat("REFERENCE CELLS BEFORE BALANCING\n")
cat("============================================================\n\n")

reference_by_sample <- table(
  reference_meta$sample_id,
  reference_meta$cell_type
)

print(reference_by_sample)

cat("\n")

reference_totals <- table(
  reference_meta$sample_id
)

print(reference_totals)

# ------------------------------------------------------------
# SAMPLE-AWARE BALANCING
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SAMPLE-AWARE REFERENCE SELECTION\n")
cat("============================================================\n\n")

selected_cells <- character(0)

selection_summary <- data.frame(
  sample_id = character(),
  available_reference_cells = integer(),
  selected_reference_cells = integer(),
  stringsAsFactors = FALSE
)

for (sample in names(table(meta$sample_id))) {

  sample_reference <- reference_meta[
    reference_meta$sample_id == sample,
    ,
    drop = FALSE
  ]

  n_available <- nrow(sample_reference)

  n_select <- min(
    n_available,
    max_cells_per_sample
  )

  if (n_available <= max_cells_per_sample) {

    selected <- sample_reference$cell

  } else {

    selected <- sample(
      sample_reference$cell,
      size = max_cells_per_sample,
      replace = FALSE
    )
  }

  selected_cells <- c(
    selected_cells,
    selected
  )

  selection_summary <- rbind(
    selection_summary,
    data.frame(
      sample_id = sample,
      available_reference_cells = n_available,
      selected_reference_cells = length(selected),
      stringsAsFactors = FALSE
    )
  )
}

# ------------------------------------------------------------
# SELECTED REFERENCE METADATA
# ------------------------------------------------------------

selected_meta <- meta[
  meta$cell %in% selected_cells,
  ,
  drop = FALSE
]

# Restore counts order
selected_meta <- selected_meta[
  match(
    selected_cells,
    selected_meta$cell
  ),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# VALIDATION
# ------------------------------------------------------------

if (nrow(selected_meta) != length(selected_cells)) {
  stop("ERROR: Reference cell selection failed.")
}

if (anyDuplicated(selected_cells) > 0) {
  stop("ERROR: Duplicate selected reference cells.")
}

cat("\nSelected reference cells: ",
    nrow(selected_meta),
    "\n\n")

# ------------------------------------------------------------
# FINAL SUMMARY
# ------------------------------------------------------------

cat("============================================================\n")
cat("FINAL REFERENCE SUMMARY\n")
cat("============================================================\n\n")

cat("Reference cell types:\n")
print(table(selected_meta$cell_type))

cat("\nReference cells by sample:\n")
print(table(selected_meta$sample_id))

cat("\nReference cells by sample and cell type:\n")
print(
  table(
    selected_meta$sample_id,
    selected_meta$cell_type
  )
)

# ------------------------------------------------------------
# REFERENCE FRACTION
# ------------------------------------------------------------

cat("\nReference fraction of all cells:\n\n")

sample_totals <- table(meta$sample_id)

selected_totals <- table(
  selected_meta$sample_id
)

reference_fraction <- data.frame(
  sample_id = names(sample_totals),
  total_cells = as.integer(sample_totals),
  selected_reference_cells = as.integer(
    selected_totals[names(sample_totals)]
  )
)

reference_fraction$selected_reference_cells[
  is.na(reference_fraction$selected_reference_cells)
] <- 0

reference_fraction$reference_percent <-
  100 *
  reference_fraction$selected_reference_cells /
  reference_fraction$total_cells

reference_fraction$reference_percent <-
  round(
    reference_fraction$reference_percent,
    3
  )

print(reference_fraction)

# ------------------------------------------------------------
# EXTRACT REFERENCE COUNTS
# ------------------------------------------------------------

cat("\nExtracting reference count matrix...\n")

reference_counts <- counts[
  ,
  selected_cells,
  drop = FALSE
]

# ------------------------------------------------------------
# FINAL MATRIX CHECK
# ------------------------------------------------------------

if (!identical(
  colnames(reference_counts),
  selected_cells
)) {
  stop(
    "ERROR: Reference counts and selected cell IDs do not match."
  )
}

if (!identical(
  colnames(reference_counts),
  selected_meta$cell
)) {
  stop(
    "ERROR: Reference counts and metadata are not aligned."
  )
}

cat("Reference counts dimensions:\n")
cat(
  nrow(reference_counts),
  "genes x",
  ncol(reference_counts),
  "cells\n"
)

cat("Reference matrix alignment: PASSED\n")

# ------------------------------------------------------------
# SAVE REFERENCE COUNTS
# ------------------------------------------------------------

reference_counts_file <- file.path(
  output_dir,
  "CNV_reference_counts_GENCODE_v47.rds"
)

saveRDS(
  reference_counts,
  reference_counts_file
)

# ------------------------------------------------------------
# SAVE REFERENCE METADATA
# ------------------------------------------------------------

reference_metadata_file <- file.path(
  output_dir,
  "CNV_reference_metadata_selected.csv"
)

write.csv(
  selected_meta,
  reference_metadata_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# SAVE REFERENCE CELL IDS
# ------------------------------------------------------------

reference_ids_file <- file.path(
  output_dir,
  "CNV_reference_cell_ids_selected.txt"
)

write.table(
  selected_cells,
  reference_ids_file,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# ------------------------------------------------------------
# SAVE SELECTION SUMMARY
# ------------------------------------------------------------

selection_summary_file <- file.path(
  output_dir,
  "CNV_reference_selection_by_sample.csv"
)

write.csv(
  selection_summary,
  selection_summary_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# SAVE REFERENCE COMPOSITION
# ------------------------------------------------------------

reference_composition <- as.data.frame(
  table(
    selected_meta$sample_id,
    selected_meta$cell_type
  )
)

colnames(reference_composition) <- c(
  "sample_id",
  "cell_type",
  "cells"
)

write.csv(
  reference_composition,
  file.path(
    output_dir,
    "CNV_reference_composition.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# SAVE SUMMARY TEXT
# ------------------------------------------------------------

summary_file <- file.path(
  output_dir,
  "CNV_reference_selection_summary.txt"
)

sink(summary_file)

cat("CNV REFERENCE SELECTION\n")
cat("=======================\n\n")

cat("Reference populations:\n")
print(available_reference_types)

cat("\nMaximum cells per sample:",
    max_cells_per_sample,
    "\n\n")

cat("Selection by sample:\n")
print(selection_summary)

cat("\nFinal reference composition:\n")
print(
  table(
    selected_meta$sample_id,
    selected_meta$cell_type
  )
)

cat("\nReference fractions:\n")
print(reference_fraction)

cat("\nTotal selected reference cells:",
    nrow(selected_meta),
    "\n")

sink()

# ------------------------------------------------------------
# FINAL MESSAGE
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("CNV REFERENCE SELECTION COMPLETED\n")
cat("============================================================\n\n")

cat("Reference populations:\n")
cat(paste(
  available_reference_types,
  collapse = ", "
))
cat("\n\n")

cat("Maximum reference cells per sample: ",
    max_cells_per_sample,
    "\n")

cat("Total selected reference cells: ",
    nrow(selected_meta),
    "\n\n")

cat("Saved files:\n")
cat(reference_counts_file, "\n")
cat(reference_metadata_file, "\n")
cat(reference_ids_file, "\n")
cat(selection_summary_file, "\n")
cat(
  file.path(
    output_dir,
    "CNV_reference_composition.csv"
  ),
  "\n"
)
cat(summary_file, "\n\n")

cat("Next step:\n")
cat("Validate the balanced reference matrix before CNV inference.\n\n")
