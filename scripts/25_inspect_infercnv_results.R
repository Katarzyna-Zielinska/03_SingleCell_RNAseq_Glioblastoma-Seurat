#!/usr/bin/env Rscript

# ============================================================
# 25_inspect_infercnv_results.R
#
# Project 03:
# Single-cell RNA-seq Glioblastoma
#
# Purpose:
#   Quantitative inspection of the FULL inferCNV result.
#
# IMPORTANT:
#   This script does NOT classify cells as malignant.
#   It only describes the inferCNV output and prepares
#   quantitative summaries for downstream validation.
# ============================================================


suppressPackageStartupMessages({
  library(infercnv)
  library(Matrix)
})


cat("\n")
cat("============================================================\n")
cat("FULL INFERCNV RESULT INSPECTION\n")
cat("============================================================\n\n")


# ============================================================
# 1. VERSION
# ============================================================

cat(
  "infercnv version: ",
  as.character(packageVersion("infercnv")),
  "\n\n",
  sep = ""
)


# ============================================================
# 2. PATHS
# ============================================================

project_dir <- getwd()

infercnv_file <- file.path(
  project_dir,
  "results",
  "infercnv",
  "GBM_3samples_infercnv_full_final.rds"
)

metadata_file <- file.path(
  project_dir,
  "results",
  "cnv_input",
  "CNV_cell_metadata.csv"
)

output_dir <- file.path(
  project_dir,
  "results",
  "infercnv",
  "inspection_full"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


cat("Input inferCNV object:\n")
cat(infercnv_file, "\n\n")

cat("Input cell metadata:\n")
cat(metadata_file, "\n\n")

cat("Output directory:\n")
cat(output_dir, "\n\n")


# ============================================================
# 3. CHECK INPUT FILES
# ============================================================

if (!file.exists(infercnv_file)) {
  stop(
    "ERROR: inferCNV object not found:\n",
    infercnv_file
  )
}

if (!file.exists(metadata_file)) {
  stop(
    "ERROR: cell metadata not found:\n",
    metadata_file
  )
}


# ============================================================
# 4. LOAD INFERCNV OBJECT
# ============================================================

cat("Loading FULL inferCNV object...\n")

infercnv_obj <- readRDS(infercnv_file)

cat("Object loaded successfully.\n\n")


cat("Class:\n")
print(class(infercnv_obj))
cat("\n")

cat("Slots:\n")
print(slotNames(infercnv_obj))
cat("\n")


# ============================================================
# 5. EXPRESSION MATRIX
# ============================================================

if (!"expr.data" %in% slotNames(infercnv_obj)) {
  stop(
    "ERROR: inferCNV object does not contain expr.data."
  )
}

expr <- infercnv_obj@expr.data


cat("Expression matrix:\n")
cat("Genes:", nrow(expr), "\n")
cat("Cells:", ncol(expr), "\n\n")

cat("Expression matrix class:\n")
print(class(expr))
cat("\n")


# ============================================================
# 6. COUNT MATRIX
# ============================================================

if ("count.data" %in% slotNames(infercnv_obj)) {

  count_data <- infercnv_obj@count.data

  cat("Count matrix:\n")
  cat("Genes:", nrow(count_data), "\n")
  cat("Cells:", ncol(count_data), "\n\n")

  rm(count_data)

} else {

  cat("Count matrix slot not available.\n\n")

}


# ============================================================
# 7. GENE ORDER
# ============================================================

if ("gene_order" %in% slotNames(infercnv_obj)) {

  gene_order <- infercnv_obj@gene_order

  cat("Gene-order object:\n")
  cat("Rows:", nrow(gene_order), "\n")
  cat("Columns:", ncol(gene_order), "\n\n")

  cat("Column names:\n")
  print(colnames(gene_order))
  cat("\n")

  cat("First 10 rows:\n")
  print(head(gene_order, 10))
  cat("\n")

  if ("chr" %in% colnames(gene_order)) {

    cat("Genes per chromosome:\n")
    print(table(gene_order$chr))
    cat("\n")

  }

} else {

  warning(
    "Gene-order slot not available."
  )

}


# ============================================================
# 8. REFERENCE / OBSERVATION GROUPS
# ============================================================

cat("============================================================\n")
cat("CELL GROUPS\n")
cat("============================================================\n\n")


reference_indices <-
  infercnv_obj@reference_grouped_cell_indices

observation_indices <-
  infercnv_obj@observation_grouped_cell_indices


reference_cells <- unlist(
  reference_indices,
  use.names = FALSE
)

observation_cells <- unlist(
  observation_indices,
  use.names = FALSE
)


cat(
  "Reference cells:",
  length(reference_cells),
  "\n"
)

cat(
  "Observation cells:",
  length(observation_cells),
  "\n"
)

cat(
  "Total cells:",
  ncol(expr),
  "\n\n"
)


if (
  length(reference_cells) +
    length(observation_cells) !=
  ncol(expr)
) {

  warning(
    "Reference + observation cells do not equal total cells."
  )

}


# ============================================================
# 9. CELL IDENTIFIERS
# ============================================================

cell_ids <- colnames(expr)

if (is.null(cell_ids)) {
  stop(
    "ERROR: inferCNV expression matrix has no cell identifiers."
  )
}


cat(
  "Number of cell identifiers:",
  length(cell_ids),
  "\n\n"
)

cat("First 10 cell identifiers:\n")
print(head(cell_ids, 10))
cat("\n")

cat("Last 10 cell identifiers:\n")
print(tail(cell_ids, 10))
cat("\n")


# ============================================================
# 10. LOAD ORIGINAL METADATA
# ============================================================

cat("============================================================\n")
cat("CELL METADATA\n")
cat("============================================================\n\n")


meta <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat(
  "Metadata rows:",
  nrow(meta),
  "\n"
)

cat(
  "Metadata columns:",
  ncol(meta),
  "\n\n"
)


cat("Metadata column names:\n")
print(colnames(meta))
cat("\n")


# First column contains cell IDs.
metadata_cell_ids <- as.character(
  meta[[1]]
)


cat("First metadata cell IDs:\n")
print(head(metadata_cell_ids))
cat("\n")


# ============================================================
# 11. MATCH CELLS
# ============================================================

match_idx <- match(
  cell_ids,
  metadata_cell_ids
)

matched <- !is.na(match_idx)


cat(
  "Cells in inferCNV object:",
  length(cell_ids),
  "\n"
)

cat(
  "Cells matched to metadata:",
  sum(matched),
  "\n"
)

cat(
  "Cells missing from metadata:",
  sum(!matched),
  "\n\n"
)


if (any(!matched)) {

  cat("First missing cells:\n")

  print(
    head(
      cell_ids[!matched],
      20
    )
  )

  cat("\n")

  stop(
    "ERROR: Some inferCNV cells are missing from metadata."
  )

}


meta_infercnv <- meta[
  match_idx,
  ,
  drop = FALSE
]


# ============================================================
# 12. ALIGNMENT CHECK
# ============================================================

if (
  !identical(
    cell_ids,
    as.character(meta_infercnv[[1]])
  )
) {

  stop(
    "ERROR: Cell order mismatch after metadata alignment."
  )

}


cat("Cell identifiers: PASSED\n")
cat("Cell order: PASSED\n\n")


# ============================================================
# 13. REQUIRED METADATA
# ============================================================

required_columns <- c(
  "sample_id",
  "cell_type"
)

missing_columns <- setdiff(
  required_columns,
  colnames(meta_infercnv)
)


if (length(missing_columns) > 0) {

  stop(
    "ERROR: Missing metadata columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )

}


cat("Required metadata columns: PASSED\n\n")


# ============================================================
# 14. ADD INFERCNV GROUP
# ============================================================

infercnv_group <- rep(
  "observation",
  length(cell_ids)
)

infercnv_group[
  reference_cells
] <- "reference"


meta_infercnv$infercnv_group <-
  infercnv_group


cat("inferCNV groups:\n")
print(
  table(
    meta_infercnv$infercnv_group
  )
)
cat("\n")


# ============================================================
# 15. CELLS BY SAMPLE
# ============================================================

sample_summary <- as.data.frame(
  table(
    meta_infercnv$sample_id
  ),
  stringsAsFactors = FALSE
)

colnames(sample_summary) <- c(
  "sample_id",
  "cells"
)


cat("Cells by sample:\n")
print(sample_summary)
cat("\n")


# ============================================================
# 16. SAMPLE x INFERCNV GROUP
# ============================================================

sample_group_summary <- as.data.frame(
  table(
    meta_infercnv$sample_id,
    meta_infercnv$infercnv_group
  ),
  stringsAsFactors = FALSE
)

colnames(sample_group_summary) <- c(
  "sample_id",
  "infercnv_group",
  "cells"
)


cat("Sample x inferCNV group:\n")
print(sample_group_summary)
cat("\n")


# ============================================================
# 17. CELL TYPES
# ============================================================

celltype_summary <- as.data.frame(
  table(
    meta_infercnv$cell_type
  ),
  stringsAsFactors = FALSE
)

colnames(celltype_summary) <- c(
  "cell_type",
  "cells"
)


celltype_summary <- celltype_summary[
  order(
    -celltype_summary$cells
  ),
  ,
  drop = FALSE
]


cat(
  "Cell types represented in FULL inferCNV object:\n"
)

print(celltype_summary)
cat("\n")


# ============================================================
# 18. SAMPLE x CELL TYPE
# ============================================================

sample_celltype_summary <- as.data.frame(
  table(
    meta_infercnv$sample_id,
    meta_infercnv$cell_type
  ),
  stringsAsFactors = FALSE
)

colnames(sample_celltype_summary) <- c(
  "sample_id",
  "cell_type",
  "cells"
)


sample_celltype_summary <-
  sample_celltype_summary[
    sample_celltype_summary$cells > 0,
    ,
    drop = FALSE
  ]


cat("Sample x cell type:\n")
print(sample_celltype_summary)
cat("\n")


# ============================================================
# 19. BASIC EXPRESSION / CNV SIGNAL STATISTICS
# ============================================================

cat("============================================================\n")
cat("QUANTITATIVE CNV SIGNAL\n")
cat("============================================================\n\n")


cat(
  "Calculating per-cell mean CNV signal...\n"
)


cell_mean_signal <- Matrix::colMeans(
  expr
)


cat(
  "Calculating per-cell median CNV signal...\n"
)


cell_median_signal <- apply(
  expr,
  2,
  median
)


cat(
  "Calculating per-cell standard deviation...\n"
)


cell_sd_signal <- apply(
  expr,
  2,
  sd
)


cat(
  "Calculating per-cell minimum...\n"
)


cell_min_signal <- apply(
  expr,
  2,
  min
)


cat(
  "Calculating per-cell maximum...\n"
)


cell_max_signal <- apply(
  expr,
  2,
  max
)


cat("Signal calculations completed.\n\n")


# ============================================================
# 20. GLOBAL SIGNAL SUMMARY
# ============================================================

global_signal_summary <- data.frame(
  metric = c(
    "minimum",
    "first_quartile",
    "median",
    "mean",
    "third_quartile",
    "maximum"
  ),
  value = c(
    min(cell_mean_signal),
    as.numeric(
      quantile(
        cell_mean_signal,
        0.25
      )
    ),
    median(cell_mean_signal),
    mean(cell_mean_signal),
    as.numeric(
      quantile(
        cell_mean_signal,
        0.75
      )
    ),
    max(cell_mean_signal)
  )
)


cat(
  "Overall per-cell mean CNV signal:\n"
)

print(global_signal_summary)
cat("\n")


# ============================================================
# 21. SIGNAL BY INFERCNV GROUP
# ============================================================

group_split <- split(
  seq_along(cell_mean_signal),
  meta_infercnv$infercnv_group
)


group_signal_list <- lapply(
  names(group_split),
  function(group_name) {

    idx <- group_split[[group_name]]

    data.frame(
      group = group_name,
      cells = length(idx),
      mean_signal = mean(
        cell_mean_signal[idx]
      ),
      sd_signal = sd(
        cell_mean_signal[idx]
      ),
      median_signal = median(
        cell_mean_signal[idx]
      ),
      min_signal = min(
        cell_min_signal[idx]
      ),
      max_signal = max(
        cell_max_signal[idx]
      ),
      stringsAsFactors = FALSE
    )

  }
)


group_signal_summary <- do.call(
  rbind,
  group_signal_list
)

rownames(group_signal_summary) <- NULL


cat("CNV signal by inferCNV group:\n")
print(group_signal_summary)
cat("\n")


# ============================================================
# 22. SIGNAL BY SAMPLE
# ============================================================

sample_split <- split(
  seq_along(cell_mean_signal),
  meta_infercnv$sample_id
)


sample_signal_list <- lapply(
  names(sample_split),
  function(sample_name) {

    idx <- sample_split[[sample_name]]

    data.frame(
      sample_id = sample_name,
      cells = length(idx),
      mean_signal = mean(
        cell_mean_signal[idx]
      ),
      sd_signal = sd(
        cell_mean_signal[idx]
      ),
      median_signal = median(
        cell_median_signal[idx]
      ),
      min_signal = min(
        cell_min_signal[idx]
      ),
      max_signal = max(
        cell_max_signal[idx]
      ),
      stringsAsFactors = FALSE
    )

  }
)


sample_signal_summary <- do.call(
  rbind,
  sample_signal_list
)

rownames(sample_signal_summary) <- NULL


cat("CNV signal by sample:\n")
print(sample_signal_summary)
cat("\n")


# ============================================================
# 23. SIGNAL BY CELL TYPE
# ============================================================

celltype_split <- split(
  seq_along(cell_mean_signal),
  meta_infercnv$cell_type
)


celltype_signal_list <- lapply(
  names(celltype_split),
  function(celltype_name) {

    idx <- celltype_split[[celltype_name]]

    data.frame(
      cell_type = celltype_name,
      cells = length(idx),
      mean_signal = mean(
        cell_mean_signal[idx]
      ),
      sd_signal = sd(
        cell_mean_signal[idx]
      ),
      median_signal = median(
        cell_median_signal[idx]
      ),
      min_signal = min(
        cell_min_signal[idx]
      ),
      max_signal = max(
        cell_max_signal[idx]
      ),
      stringsAsFactors = FALSE
    )

  }
)


celltype_signal_summary <- do.call(
  rbind,
  celltype_signal_list
)

rownames(celltype_signal_summary) <- NULL


celltype_signal_summary <-
  celltype_signal_summary[
    order(
      -celltype_signal_summary$cells
    ),
    ,
    drop = FALSE
  ]


cat("CNV signal by cell type:\n")
print(celltype_signal_summary)
cat("\n")


# ============================================================
# 24. PER-CELL QUANTITATIVE TABLE
# ============================================================

cell_signal_table <- data.frame(
  cell_id = cell_ids,
  sample_id = meta_infercnv$sample_id,
  cell_type = meta_infercnv$cell_type,
  infercnv_group = meta_infercnv$infercnv_group,
  mean_signal = as.numeric(
    cell_mean_signal
  ),
  median_signal = as.numeric(
    cell_median_signal
  ),
  sd_signal = as.numeric(
    cell_sd_signal
  ),
  min_signal = as.numeric(
    cell_min_signal
  ),
  max_signal = as.numeric(
    cell_max_signal
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 25. SAVE PER-CELL TABLE
# ============================================================

write.csv(
  cell_signal_table,
  file.path(
    output_dir,
    "infercnv_full_per_cell_signal.csv"
  ),
  row.names = FALSE
)


saveRDS(
  cell_signal_table,
  file.path(
    output_dir,
    "infercnv_full_per_cell_signal.rds"
  )
)


# ============================================================
# 26. SAVE SUMMARY TABLES
# ============================================================

write.csv(
  sample_summary,
  file.path(
    output_dir,
    "infercnv_full_cells_by_sample.csv"
  ),
  row.names = FALSE
)


write.csv(
  sample_group_summary,
  file.path(
    output_dir,
    "infercnv_full_sample_x_group.csv"
  ),
  row.names = FALSE
)


write.csv(
  celltype_summary,
  file.path(
    output_dir,
    "infercnv_full_cells_by_celltype.csv"
  ),
  row.names = FALSE
)


write.csv(
  sample_celltype_summary,
  file.path(
    output_dir,
    "infercnv_full_sample_x_celltype.csv"
  ),
  row.names = FALSE
)


write.csv(
  global_signal_summary,
  file.path(
    output_dir,
    "infercnv_full_global_signal_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  group_signal_summary,
  file.path(
    output_dir,
    "infercnv_full_signal_by_group.csv"
  ),
  row.names = FALSE
)


write.csv(
  sample_signal_summary,
  file.path(
    output_dir,
    "infercnv_full_signal_by_sample.csv"
  ),
  row.names = FALSE
)


write.csv(
  celltype_signal_summary,
  file.path(
    output_dir,
    "infercnv_full_signal_by_celltype.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 27. SAVE ALIGNED METADATA
# ============================================================

meta_output <- meta_infercnv

meta_output$mean_signal <-
  as.numeric(cell_mean_signal)

meta_output$median_signal <-
  as.numeric(cell_median_signal)

meta_output$sd_signal <-
  as.numeric(cell_sd_signal)

meta_output$min_signal <-
  as.numeric(cell_min_signal)

meta_output$max_signal <-
  as.numeric(cell_max_signal)


write.csv(
  meta_output,
  file.path(
    output_dir,
    "infercnv_full_metadata_aligned.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 28. SUMMARY TEXT FILE
# ============================================================

summary_file <- file.path(
  output_dir,
  "infercnv_full_inspection_summary.txt"
)


sink(summary_file)


cat(
  "FULL INFERCNV RESULT INSPECTION\n"
)

cat(
  "================================\n\n"
)


cat(
  "infercnv version: ",
  as.character(packageVersion("infercnv")),
  "\n\n",
  sep = ""
)


cat(
  "Input object:\n"
)

cat(
  infercnv_file,
  "\n\n"
)


cat(
  "Expression matrix:\n"
)

cat(
  "Genes:",
  nrow(expr),
  "\n"
)

cat(
  "Cells:",
  ncol(expr),
  "\n\n"
)


cat(
  "Reference cells:",
  length(reference_cells),
  "\n"
)

cat(
  "Observation cells:",
  length(observation_cells),
  "\n\n"
)


cat(
  "Cells by sample:\n"
)

print(sample_summary)


cat(
  "\nSample x inferCNV group:\n"
)

print(sample_group_summary)


cat(
  "\nCell types:\n"
)

print(celltype_summary)


cat(
  "\nOverall per-cell mean CNV signal:\n"
)

print(global_signal_summary)


cat(
  "\nCNV signal by inferCNV group:\n"
)

print(group_signal_summary)


cat(
  "\nCNV signal by sample:\n"
)

print(sample_signal_summary)


cat(
  "\nCNV signal by cell type:\n"
)

print(celltype_signal_summary)


cat(
  "\nIMPORTANT:\n"
)

cat(
  "These results are descriptive quantitative summaries of the\n",
  "FULL inferCNV output. They do NOT classify cells as malignant.\n",
  "No malignant-cell threshold is applied in this script.\n",
  "The next step is validation of the full observation population.\n",
  sep = ""
)


sink()


# ============================================================
# 29. FINAL MESSAGE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SCRIPT 25 COMPLETED\n")
cat("============================================================\n\n")

cat("FULL inferCNV object inspected successfully.\n\n")

cat("Reference cells:    ",
    length(reference_cells),
    "\n",
    sep = "")

cat("Observation cells:  ",
    length(observation_cells),
    "\n",
    sep = "")

cat("Total cells:        ",
    ncol(expr),
    "\n",
    sep = "")

cat("Genes in expression matrix: ",
    nrow(expr),
    "\n\n",
    sep = "")

cat("Output directory:\n")
cat(output_dir, "\n\n")

cat("Saved files include:\n")
cat("- infercnv_full_per_cell_signal.csv\n")
cat("- infercnv_full_per_cell_signal.rds\n")
cat("- infercnv_full_cells_by_sample.csv\n")
cat("- infercnv_full_sample_x_group.csv\n")
cat("- infercnv_full_cells_by_celltype.csv\n")
cat("- infercnv_full_sample_x_celltype.csv\n")
cat("- infercnv_full_global_signal_summary.csv\n")
cat("- infercnv_full_signal_by_group.csv\n")
cat("- infercnv_full_signal_by_sample.csv\n")
cat("- infercnv_full_signal_by_celltype.csv\n")
cat("- infercnv_full_metadata_aligned.csv\n")
cat("- infercnv_full_inspection_summary.txt\n\n")

cat("IMPORTANT:\n")
cat(
  "Do NOT classify cells as malignant based on this script alone.\n",
  sep = ""
)

cat("\n")
cat("Next step: FULL observation validation (script 26).\n\n")
