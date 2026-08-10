#!/usr/bin/env Rscript

# ============================================================
# INFERCNV CNV INFERENCE
# Project 03 - Single-cell RNA-seq Glioblastoma
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(infercnv)
})

options(scipen = 100)

cat("\n")
cat("============================================================\n")
cat("INFERCNV CNV INFERENCE\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. VERSIONS
# ------------------------------------------------------------

cat(
  "Seurat version: ",
  as.character(packageVersion("Seurat")),
  "\n"
)

cat(
  "SeuratObject version: ",
  as.character(packageVersion("SeuratObject")),
  "\n"
)

cat(
  "infercnv version: ",
  as.character(packageVersion("infercnv")),
  "\n\n"
)


# ------------------------------------------------------------
# 2. PATHS
# ------------------------------------------------------------

project_dir <- "/home/katja/Cancer-Bioinformatics-Portfolio/03_SingleCell-RNAseq-Glioblastoma-Seurat"

counts_file <- file.path(
  project_dir,
  "results/cnv_input/CNV_raw_counts_GENCODE_v47.rds"
)

gene_order_csv <- file.path(
  project_dir,
  "results/cnv_input/CNV_gene_order_GENCODE_v47.csv"
)

metadata_file <- file.path(
  project_dir,
  "results/cnv_input/CNV_cell_metadata.csv"
)

reference_ids_file <- file.path(
  project_dir,
  "results/cnv_input/CNV_reference_cell_ids_selected.txt"
)

output_dir <- file.path(
  project_dir,
  "results/infercnv"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Input files:\n")
cat(counts_file, "\n")
cat(gene_order_csv, "\n")
cat(metadata_file, "\n")
cat(reference_ids_file, "\n\n")

cat("Output directory:\n")
cat(output_dir, "\n\n")


# ------------------------------------------------------------
# 3. CHECK INPUT FILES
# ------------------------------------------------------------

required_files <- c(
  counts_file,
  gene_order_csv,
  metadata_file,
  reference_ids_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  cat("ERROR: Missing input files:\n")

  for (f in missing_files) {
    cat("-", f, "\n")
  }

  stop("Required input files are missing.")
}


# ------------------------------------------------------------
# 4. LOAD RAW COUNTS
# ------------------------------------------------------------

cat("Loading raw counts...\n")

counts <- readRDS(counts_file)

cat("Counts loaded.\n")
cat("Genes:", nrow(counts), "\n")
cat("Cells:", ncol(counts), "\n\n")

if (!inherits(counts, "Matrix")) {
  counts <- as(counts, "dgCMatrix")
}

cat("Counts class:\n")
print(class(counts))
cat("\n")


# ------------------------------------------------------------
# 5. LOAD METADATA
# ------------------------------------------------------------

cat("Loading cell metadata...\n")

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "Metadata rows:",
  nrow(metadata),
  "\n"
)

cat(
  "Metadata columns:",
  ncol(metadata),
  "\n\n"
)

cat("Metadata column names:\n")
print(colnames(metadata))

cat("\n")


# ------------------------------------------------------------
# 6. IDENTIFY CELL ID COLUMN
# ------------------------------------------------------------

cell_id_col <- 1

cell_ids <- as.character(
  metadata[[cell_id_col]]
)

if (
  any(is.na(cell_ids)) ||
  any(cell_ids == "")
) {
  stop("ERROR: Invalid cell IDs in metadata.")
}

cat(
  "Metadata cell ID column:",
  colnames(metadata)[cell_id_col],
  "\n\n"
)

cat("First metadata cell IDs:\n")
print(head(cell_ids))

cat("\n")


# ------------------------------------------------------------
# 7. CELL ALIGNMENT
# ------------------------------------------------------------

count_cell_ids <- colnames(counts)

if (is.null(count_cell_ids)) {
  stop("ERROR: Counts matrix has no cell identifiers.")
}

cat("Checking cell identifiers...\n")

if (!setequal(count_cell_ids, cell_ids)) {

  missing_in_metadata <- setdiff(
    count_cell_ids,
    cell_ids
  )

  missing_in_counts <- setdiff(
    cell_ids,
    count_cell_ids
  )

  cat(
    "Cells missing from metadata:",
    length(missing_in_metadata),
    "\n"
  )

  cat(
    "Cells missing from counts:",
    length(missing_in_counts),
    "\n"
  )

  stop("ERROR: Cell identifiers do not match.")
}


if (!identical(count_cell_ids, cell_ids)) {

  cat(
    "Cell order differs.\n"
  )

  cat(
    "Reordering metadata to match counts...\n"
  )

  metadata <- metadata[
    match(
      count_cell_ids,
      cell_ids
    ),
    ,
    drop = FALSE
  ]

  cell_ids <- as.character(
    metadata[[cell_id_col]]
  )
}


if (!identical(count_cell_ids, cell_ids)) {
  stop(
    "ERROR: Cell order could not be aligned."
  )
}

cat("Cell identifiers: PASSED\n")
cat("Cell order: PASSED\n\n")


# ------------------------------------------------------------
# 8. LOAD GENE ORDER
# ------------------------------------------------------------

cat("Loading gene order...\n")

gene_order <- read.csv(
  gene_order_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "Gene order rows:",
  nrow(gene_order),
  "\n"
)

cat("Gene order columns:\n")
print(colnames(gene_order))

cat("\n")


required_gene_columns <- c(
  "chromosome",
  "start",
  "end",
  "seurat_gene"
)

missing_gene_columns <- setdiff(
  required_gene_columns,
  colnames(gene_order)
)

if (length(missing_gene_columns) > 0) {

  stop(
    paste(
      "ERROR: Missing gene-order columns:",
      paste(
        missing_gene_columns,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 9. PREPARE INFERCNV GENE ORDER
# ------------------------------------------------------------

cat(
  "Preparing inferCNV gene-order table...\n"
)

gene_order_infercnv <- data.frame(
  chr = as.character(
    gene_order$chromosome
  ),
  start = as.numeric(
    gene_order$start
  ),
  stop = as.numeric(
    gene_order$end
  ),
  row.names = as.character(
    gene_order$seurat_gene
  ),
  stringsAsFactors = FALSE
)


duplicated_genes <- duplicated(
  rownames(gene_order_infercnv)
)

if (any(duplicated_genes)) {

  cat(
    "Removing duplicated genes:",
    sum(duplicated_genes),
    "\n"
  )

  gene_order_infercnv <- gene_order_infercnv[
    !duplicated_genes,
    ,
    drop = FALSE
  ]
}


standard_chr <- paste0(
  "chr",
  c(
    1:22,
    "X",
    "Y"
  )
)

keep_chr <- gene_order_infercnv$chr %in% standard_chr

gene_order_infercnv <- gene_order_infercnv[
  keep_chr,
  ,
  drop = FALSE
]


chr_levels <- paste0(
  "chr",
  c(
    1:22,
    "X",
    "Y"
  )
)

gene_order_infercnv$chr <- factor(
  gene_order_infercnv$chr,
  levels = chr_levels
)

gene_order_infercnv <- gene_order_infercnv[
  order(
    gene_order_infercnv$chr,
    gene_order_infercnv$start
  ),
  ,
  drop = FALSE
]

gene_order_infercnv$chr <- as.character(
  gene_order_infercnv$chr
)

cat(
  "Genes in prepared gene-order table:",
  nrow(gene_order_infercnv),
  "\n\n"
)


# ------------------------------------------------------------
# 10. MATCH COUNTS TO GENE ORDER
# ------------------------------------------------------------

count_genes <- rownames(counts)

if (is.null(count_genes)) {
  stop(
    "ERROR: Counts matrix has no gene identifiers."
  )
}

common_genes <- intersect(
  rownames(gene_order_infercnv),
  count_genes
)

cat(
  "Genes in counts:",
  length(count_genes),
  "\n"
)

cat(
  "Genes in gene-order:",
  nrow(gene_order_infercnv),
  "\n"
)

cat(
  "Genes shared:",
  length(common_genes),
  "\n\n"
)


if (length(common_genes) < 1000) {

  stop(
    "ERROR: Too few genes shared between counts and gene order."
  )
}


gene_order_infercnv <- gene_order_infercnv[
  common_genes,
  ,
  drop = FALSE
]

counts <- counts[
  rownames(gene_order_infercnv),
  ,
  drop = FALSE
]

cat(
  "Final counts matrix:",
  nrow(counts),
  "genes x",
  ncol(counts),
  "cells\n"
)

cat(
  "Final gene-order table:",
  nrow(gene_order_infercnv),
  "genes\n\n"
)


if (
  !identical(
    rownames(counts),
    rownames(gene_order_infercnv)
  )
) {

  stop(
    "ERROR: Gene order mismatch."
  )
}

cat(
  "Gene order alignment: PASSED\n\n"
)


# ------------------------------------------------------------
# 11. WRITE INFERCNV GENE ORDER
# ------------------------------------------------------------

infercnv_gene_order_file <- file.path(
  output_dir,
  "infercnv_gene_order.txt"
)

write.table(
  gene_order_infercnv,
  file = infercnv_gene_order_file,
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = TRUE
)

cat(
  "Gene-order file saved:\n"
)

cat(
  infercnv_gene_order_file,
  "\n\n"
)


# ------------------------------------------------------------
# 12. LOAD REFERENCE CELL IDS
# ------------------------------------------------------------

cat(
  "Loading reference cell IDs...\n"
)

reference_cells <- readLines(
  reference_ids_file
)

reference_cells <- trimws(
  reference_cells
)

reference_cells <- reference_cells[
  reference_cells != ""
]

reference_cells <- unique(
  reference_cells
)

cat(
  "Reference cells:",
  length(reference_cells),
  "\n"
)


missing_reference <- setdiff(
  reference_cells,
  colnames(counts)
)

if (length(missing_reference) > 0) {

  cat(
    "Reference cells missing:",
    length(missing_reference),
    "\n"
  )

  stop(
    "ERROR: Some reference cells are missing from counts."
  )
}

cat(
  "Reference cell IDs: PASSED\n\n"
)


# ------------------------------------------------------------
# 13. CREATE CELL ANNOTATION
# ------------------------------------------------------------

cat(
  "Creating inferCNV cell annotations...\n"
)

annotation <- data.frame(
  cell_id = colnames(counts),
  group = ifelse(
    colnames(counts) %in% reference_cells,
    "normal_reference",
    "observation"
  ),
  stringsAsFactors = FALSE
)

cat(
  "Annotation groups:\n"
)

print(
  table(annotation$group)
)

cat("\n")


annotation_file <- file.path(
  output_dir,
  "infercnv_cell_annotations.txt"
)

write.table(
  annotation,
  file = annotation_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat(
  "Annotation file saved:\n"
)

cat(
  annotation_file,
  "\n\n"
)


# ------------------------------------------------------------
# 14. FINAL CELL COUNTS
# ------------------------------------------------------------

reference_count <- sum(
  annotation$group == "normal_reference"
)

observation_count <- sum(
  annotation$group == "observation"
)

cat(
  "Reference cells:",
  reference_count,
  "\n"
)

cat(
  "Observation cells:",
  observation_count,
  "\n"
)

cat(
  "Total cells:",
  reference_count +
    observation_count,
  "\n\n"
)


if (
  reference_count !=
    length(reference_cells)
) {

  stop(
    "ERROR: Reference cell count mismatch."
  )
}


if (
  reference_count +
    observation_count !=
    ncol(counts)
) {

  stop(
    "ERROR: Cell count mismatch."
  )
}


# ------------------------------------------------------------
# 15. SAMPLE REPRESENTATION
# ------------------------------------------------------------

if (
  !"sample_id" %in%
  colnames(metadata)
) {

  stop(
    "ERROR: sample_id column not found."
  )
}


cat(
  "Cells per sample:\n"
)

print(
  table(metadata$sample_id)
)

cat("\n")


reference_metadata <- metadata[
  metadata[[cell_id_col]]
    %in%
    reference_cells,
  ,
  drop = FALSE
]

cat(
  "Reference cells per sample:\n"
)

print(
  table(
    reference_metadata$sample_id
  )
)

cat("\n")


observation_metadata <- metadata[
  !(
    metadata[[cell_id_col]]
      %in%
      reference_cells
  ),
  ,
  drop = FALSE
]

cat(
  "Observation cells per sample:\n"
)

print(
  table(
    observation_metadata$sample_id
  )
)

cat("\n")


# ------------------------------------------------------------
# 16. CREATE INFERCNV OBJECT
# ------------------------------------------------------------

cat(
  "Creating inferCNV object...\n\n"
)

cat(
  "IMPORTANT:\n"
)

cat(
  "max_cells_per_group = NULL\n"
)

cat(
  "Full observation population retained.\n\n"
)


infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = counts,
  gene_order_file = infercnv_gene_order_file,
  annotations_file = annotation_file,
  ref_group_names = "normal_reference",
  delim = "\t",
  max_cells_per_group = NULL,
  min_max_counts_per_cell = c(
    100,
    Inf
  ),
  chr_exclude = c(
    "chrX",
    "chrY",
    "chrM"
  )
)


cat("\n")
cat(
  "inferCNV object created successfully.\n\n"
)

cat(
  "Genes:",
  nrow(infercnv_obj@expr.data),
  "\n"
)

cat(
  "Cells:",
  ncol(infercnv_obj@expr.data),
  "\n\n"
)


# ------------------------------------------------------------
# 17. VERIFY CELL RETENTION
# ------------------------------------------------------------

infercnv_cells <- colnames(
  infercnv_obj@expr.data
)

cat(
  "Cell retention check:\n"
)

cat(
  "Original cells:",
  ncol(counts),
  "\n"
)

cat(
  "inferCNV cells:",
  length(infercnv_cells),
  "\n"
)

missing_after_creation <- setdiff(
  colnames(counts),
  infercnv_cells
)

cat(
  "Cells missing after CreateInfercnvObject:",
  length(missing_after_creation),
  "\n\n"
)


if (
  length(missing_after_creation) > 0
) {

  cat(
    "WARNING: Some cells were removed during object creation.\n"
  )

} else {

  cat(
    "All input cells retained: PASSED\n\n"
  )
}


# ------------------------------------------------------------
# 18. RUN INFERCNV
# ------------------------------------------------------------

cat(
  "============================================================\n"
)

cat(
  "STARTING INFERCNV::RUN()\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "Analysis mode: samples\n"
)

cat(
  "Window length: 101\n"
)

cat(
  "Cutoff: 0.1\n"
)

cat(
  "Denoising: TRUE\n"
)

cat(
  "HMM: FALSE\n"
)

cat(
  "Reference clustering: DISABLED\n"
)

cat(
  "num_ref_groups: NULL\n"
)

cat(
  "Threads: 4\n\n"
)


result <- infercnv::run(

  infercnv_obj,

  cutoff = 0.1,

  min_cells_per_gene = 3,

  out_dir = output_dir,

  window_length = 101,

  smooth_method = "pyramidinal",

  num_ref_groups = NULL,

  ref_subtract_use_mean_bounds = TRUE,

  cluster_by_groups = TRUE,

  cluster_references = FALSE,

  k_obs_groups = 1,

  hclust_method = "ward.D2",

  max_centered_threshold = 3,

  scale_data = FALSE,

  HMM = FALSE,

  denoise = TRUE,

  noise_filter = NA,

  sim_method = "meanvar",

  sim_foreground = FALSE,

  reassignCNVs = TRUE,

  analysis_mode = "samples",

  tumor_subcluster_partition_method = "leiden",

  tumor_subcluster_pval = 0.1,

  k_nn = 20,

  leiden_method = "PCA",

  leiden_function = "modularity",

  leiden_resolution = "auto",

  z_score_filter = 0.8,

  sd_amplifier = 1.5,

  outlier_method_bound = "average_bound",

  prune_outliers = FALSE,

  mask_nonDE_genes = FALSE,

  test.use = "wilcoxon",

  require_DE_all_normals = "any",

  hspike_aggregate_normals = FALSE,

  no_plot = FALSE,

  no_prelim_plot = FALSE,

  write_expr_matrix = FALSE,

  write_phylo = FALSE,

  output_format = "png",

  plot_chr_scale = FALSE,

  useRaster = TRUE,

  num_threads = 4,

  plot_steps = FALSE,

  inspect_subclusters = TRUE,

  resume_mode = FALSE,

  png_res = 150,

  plot_probabilities = TRUE,

  save_rds = TRUE,

  save_final_rds = TRUE,

  diagnostics = TRUE,

  remove_genes_at_chr_ends = FALSE,

  up_to_step = 100
)


# ------------------------------------------------------------
# 19. SAVE FINAL OBJECT
# ------------------------------------------------------------

cat("\n")

cat(
  "============================================================\n"
)

cat(
  "INFERCNV RUN COMPLETED\n"
)

cat(
  "============================================================\n\n"
)


final_rds <- file.path(
  output_dir,
  "GBM_3samples_infercnv_full_final.rds"
)

saveRDS(
  result,
  file = final_rds
)

cat(
  "Final inferCNV object saved:\n"
)

cat(
  final_rds,
  "\n\n"
)


# ------------------------------------------------------------
# 20. FINAL SUMMARY
# ------------------------------------------------------------

cat(
  "FINAL SUMMARY\n"
)

cat(
  "-------------\n"
)

cat(
  "Input genes:",
  nrow(counts),
  "\n"
)

cat(
  "Input cells:",
  ncol(counts),
  "\n"
)

cat(
  "Reference cells:",
  reference_count,
  "\n"
)

cat(
  "Observation cells:",
  observation_count,
  "\n"
)

cat(
  "Final inferCNV genes:",
  nrow(result@expr.data),
  "\n"
)

cat(
  "Final inferCNV cells:",
  ncol(result@expr.data),
  "\n"
)

cat("\n")

cat(
  "Output directory:\n"
)

cat(
  output_dir,
  "\n\n"
)

cat(
  "Next step:\n"
)

cat(
  "Inspect the full inferCNV result before any malignant-cell classification.\n"
)

cat("\n")

cat(
  "============================================================\n"
)

cat(
  "DONE\n"
)

cat(
  "============================================================\n"
)
