#!/usr/bin/env Rscript

# ============================================================
# SAMPLE COMPOSITION ANALYSIS
# Glioblastoma scRNA-seq
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

cat("\n")
cat("============================================================\n")
cat("SAMPLE COMPOSITION ANALYSIS\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# PROJECT PATHS
# ------------------------------------------------------------

project_dir <- "/home/katja/Cancer-Bioinformatics-Portfolio/03_SingleCell-RNAseq-Glioblastoma-Seurat"

input_file <- file.path(
  project_dir,
  "results",
  "seurat_annotation",
  "GBM_3samples_Seurat_annotated.rds"
)

output_dir <- file.path(
  project_dir,
  "results",
  "sample_composition"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# PACKAGE VERSIONS
# ------------------------------------------------------------

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
# LOAD SEURAT OBJECT
# ------------------------------------------------------------

cat("Loading Seurat object:\n")
cat(input_file, "\n\n")

if (!file.exists(input_file)) {
  stop(
    paste(
      "Input file does not exist:",
      input_file
    )
  )
}

obj <- readRDS(input_file)

cat("Object loaded successfully.\n\n")

cat("Genes: ", nrow(obj), "\n", sep = "")
cat("Cells: ", ncol(obj), "\n\n", sep = "")

# ------------------------------------------------------------
# CHECK METADATA
# ------------------------------------------------------------

required_columns <- c(
  "sample_id",
  "cell_type"
)

missing_columns <- setdiff(
  required_columns,
  colnames(obj@meta.data)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "Missing metadata columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# SAMPLE COUNTS
# ------------------------------------------------------------

cat("============================================================\n")
cat("TOTAL CELLS PER SAMPLE\n")
cat("============================================================\n\n")

sample_counts <- as.data.frame(
  table(obj$sample_id)
)

colnames(sample_counts) <- c(
  "sample_id",
  "cells"
)

print(sample_counts, row.names = FALSE)

write.csv(
  sample_counts,
  file.path(
    output_dir,
    "sample_cell_counts.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# CELL TYPE COUNTS
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CELL TYPE COUNTS BY SAMPLE\n")
cat("============================================================\n\n")

count_table <- table(
  obj$sample_id,
  obj$cell_type
)

print(count_table)

count_df <- as.data.frame(count_table)

colnames(count_df) <- c(
  "sample_id",
  "cell_type",
  "cells"
)

write.csv(
  count_df,
  file.path(
    output_dir,
    "sample_celltype_counts.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# PERCENTAGE COMPOSITION
# ------------------------------------------------------------

composition_df <- count_df %>%
  group_by(sample_id) %>%
  mutate(
    total_cells = sum(cells),
    percent = 100 * cells / total_cells
  ) %>%
  ungroup()

write.csv(
  composition_df,
  file.path(
    output_dir,
    "sample_celltype_composition_percent.csv"
  ),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("CELL TYPE COMPOSITION (%)\n")
cat("============================================================\n\n")

composition_print <- composition_df %>%
  select(
    sample_id,
    cell_type,
    percent
  ) %>%
  arrange(
    sample_id,
    desc(percent)
  )

print(
  composition_print,
  row.names = FALSE
)

# ------------------------------------------------------------
# DOMINANT CELL TYPE
# ------------------------------------------------------------

dominant_df <- composition_df %>%
  group_by(sample_id) %>%
  slice_max(
    order_by = percent,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(
    sample_id,
    dominant_cell_type = cell_type,
    dominant_percent = percent,
    total_cells
  )

cat("\n============================================================\n")
cat("DOMINANT CELL TYPE PER SAMPLE\n")
cat("============================================================\n\n")

print(
  dominant_df,
  row.names = FALSE
)

write.csv(
  dominant_df,
  file.path(
    output_dir,
    "dominant_cell_type_per_sample.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# SHANNON DIVERSITY
# ------------------------------------------------------------

shannon_df <- composition_df %>%
  group_by(sample_id) %>%
  summarise(
    shannon_diversity = -sum(
      ifelse(
        percent > 0,
        (percent / 100) *
          log(percent / 100),
        0
      )
    ),
    total_cells = first(total_cells),
    .groups = "drop"
  )

cat("\n============================================================\n")
cat("SHANNON DIVERSITY\n")
cat("============================================================\n\n")

print(
  shannon_df,
  row.names = FALSE
)

write.csv(
  shannon_df,
  file.path(
    output_dir,
    "sample_shannon_diversity.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# COUNT BARPLOT
# ------------------------------------------------------------

cat("\nGenerating cell count barplot...\n")

p_count <- ggplot(
  count_df,
  aes(
    x = sample_id,
    y = cells,
    fill = cell_type
  )
) +
  geom_col() +
  labs(
    title = "Cellular Composition by Sample",
    x = "Sample",
    y = "Number of cells",
    fill = "Cell type"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "Cell_Composition_Counts.png"
  ),
  plot = p_count,
  width = 12,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# PERCENT BARPLOT
# ------------------------------------------------------------

cat("Generating percentage composition plot...\n")

p_percent <- ggplot(
  composition_df,
  aes(
    x = sample_id,
    y = percent,
    fill = cell_type
  )
) +
  geom_col() +
  labs(
    title = "Relative Cellular Composition by Sample",
    x = "Sample",
    y = "Cells (%)",
    fill = "Cell type"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "Cell_Composition_Percent.png"
  ),
  plot = p_percent,
  width = 12,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# COMPOSITION HEATMAP
# ------------------------------------------------------------

cat("Generating composition heatmap...\n")

heatmap_df <- composition_df %>%
  select(
    sample_id,
    cell_type,
    percent
  ) %>%
  pivot_wider(
    names_from = cell_type,
    values_from = percent,
    values_fill = 0
  )

heatmap_matrix <- as.data.frame(
  heatmap_df
)

rownames(heatmap_matrix) <- heatmap_matrix$sample_id

heatmap_matrix$sample_id <- NULL

heatmap_matrix <- as.matrix(
  heatmap_matrix
)

png(
  filename = file.path(
    output_dir,
    "Cell_Composition_Heatmap.png"
  ),
  width = 3200,
  height = 2200,
  res = 300
)

pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.1f",
  main = "Cellular Composition (%)"
)

dev.off()

# ------------------------------------------------------------
# FINAL SUMMARY
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("SAMPLE COMPOSITION ANALYSIS COMPLETED\n")
cat("============================================================\n\n")

cat(
  "Output directory:\n",
  output_dir,
  "\n\n",
  sep = ""
)

cat("Saved files:\n")

output_files <- list.files(
  output_dir,
  full.names = FALSE
)

for (f in output_files) {
  cat(" - ", f, "\n", sep = "")
}

cat("\n")
cat("Important:\n")
cat(
  "The three samples should be interpreted separately when\n",
  "evaluating cellular composition because strong sample-specific\n",
  "differences are present in the current annotation.\n",
  sep = ""
)

cat("\n")
cat("Next step:\n")
cat(
  "Inspect sample composition before malignant/non-malignant\n",
  "and CNV/state analysis.\n",
  sep = ""
)

cat("\n")
