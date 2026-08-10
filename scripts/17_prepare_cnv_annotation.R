#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

cat("\n")
cat("============================================================\n")
cat("PREPARING GENE ANNOTATION FOR CNV ANALYSIS\n")
cat("GENCODE v47 / GRCh38\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

project_dir <- normalizePath(".", mustWork = TRUE)

seurat_file <- file.path(
  project_dir,
  "results/seurat_annotation/GBM_3samples_Seurat_annotated.rds"
)

gtf_file <- file.path(
  project_dir,
  "data/reference/annotation/gencode.v47.primary_assembly.annotation.gtf"
)

output_dir <- file.path(
  project_dir,
  "results/cnv_annotation"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. VERSIONS
# ------------------------------------------------------------

cat(
  "Seurat version:",
  as.character(packageVersion("Seurat")),
  "\n"
)

cat(
  "SeuratObject version:",
  as.character(packageVersion("SeuratObject")),
  "\n\n"
)

# ------------------------------------------------------------
# 3. INPUT CHECK
# ------------------------------------------------------------

if (!file.exists(seurat_file)) {
  stop(
    "Seurat object not found:\n",
    seurat_file
  )
}

if (!file.exists(gtf_file)) {
  stop(
    "GTF file not found:\n",
    gtf_file
  )
}

cat("Seurat object:\n")
cat(seurat_file, "\n\n")

cat("GTF annotation:\n")
cat(gtf_file, "\n\n")

# ------------------------------------------------------------
# 4. LOAD SEURAT
# ------------------------------------------------------------

cat("Loading Seurat object...\n")

obj <- readRDS(seurat_file)

cat("Object loaded successfully.\n\n")

seurat_genes <- rownames(obj)

cat(
  "Genes in Seurat:",
  length(seurat_genes),
  "\n"
)

cat(
  "Cells in Seurat:",
  ncol(obj),
  "\n\n"
)

# ------------------------------------------------------------
# 5. GENE ID / SYMBOL SUMMARY
# ------------------------------------------------------------

is_ensg <- grepl(
  "^ENSG[0-9]+$",
  seurat_genes
)

cat(
  "ENSG genes:",
  sum(is_ensg),
  "\n"
)

cat(
  "Gene symbols / other:",
  sum(!is_ensg),
  "\n\n"
)

# ------------------------------------------------------------
# 6. READ GTF
# ------------------------------------------------------------

cat("Reading GENCODE GTF...\n")
cat("This may take some time.\n\n")

gtf <- read.delim(
  gtf_file,
  header = FALSE,
  sep = "\t",
  comment.char = "#",
  quote = "",
  stringsAsFactors = FALSE
)

colnames(gtf) <- c(
  "chromosome",
  "source",
  "feature",
  "start",
  "end",
  "score",
  "strand",
  "frame",
  "attribute"
)

cat(
  "GTF rows:",
  nrow(gtf),
  "\n\n"
)

# ------------------------------------------------------------
# 7. KEEP GENE FEATURES
# ------------------------------------------------------------

gene_gtf <- gtf[
  gtf$feature == "gene",
]

cat(
  "Gene records:",
  nrow(gene_gtf),
  "\n\n"
)

# ------------------------------------------------------------
# 8. EXTRACT GENE_ID
# ------------------------------------------------------------

gene_gtf$gene_id <- sub(
  '.*gene_id "([^"]+)".*',
  "\\1",
  gene_gtf$attribute
)

# Remove version suffix
# ENSG00000123456.7 -> ENSG00000123456

gene_gtf$gene_id <- sub(
  "\\..*$",
  "",
  gene_gtf$gene_id
)

# ------------------------------------------------------------
# 9. EXTRACT GENE_NAME
# ------------------------------------------------------------

gene_gtf$gene_name <- sub(
  '.*gene_name "([^"]+)".*',
  "\\1",
  gene_gtf$attribute
)

# ------------------------------------------------------------
# 10. EXTRACT GENE_TYPE
# ------------------------------------------------------------

gene_gtf$gene_type <- sub(
  '.*gene_type "([^"]+)".*',
  "\\1",
  gene_gtf$attribute
)

# ------------------------------------------------------------
# 11. CLEAN GTF TABLE
# ------------------------------------------------------------

gene_annotation <- gene_gtf[
  ,
  c(
    "chromosome",
    "start",
    "end",
    "strand",
    "gene_id",
    "gene_name",
    "gene_type"
  )
]

# ------------------------------------------------------------
# 12. REMOVE DUPLICATE GENE IDS
# ------------------------------------------------------------

gene_annotation <- gene_annotation[
  !duplicated(gene_annotation$gene_id),
]

cat(
  "Unique GTF gene IDs:",
  nrow(gene_annotation),
  "\n\n"
)

# ------------------------------------------------------------
# 13. MATCH BY ENSEMBL GENE ID
# ------------------------------------------------------------

ensg_genes <- seurat_genes[
  is_ensg
]

ensg_match <- match(
  ensg_genes,
  gene_annotation$gene_id
)

matched_ensg <- !is.na(ensg_match)

cat(
  "ENSG genes matched:",
  sum(matched_ensg),
  "/",
  length(ensg_genes),
  "\n"
)

# ------------------------------------------------------------
# 14. MATCH NON-ENSG GENES BY GENE SYMBOL
# ------------------------------------------------------------

symbol_genes <- seurat_genes[
  !is_ensg
]

symbol_match <- match(
  symbol_genes,
  gene_annotation$gene_name
)

matched_symbol <- !is.na(symbol_match)

cat(
  "Gene symbols matched:",
  sum(matched_symbol),
  "/",
  length(symbol_genes),
  "\n\n"
)

# ------------------------------------------------------------
# 15. BUILD MATCHED TABLE
# ------------------------------------------------------------

# ENSG matches

ensg_table <- gene_annotation[
  ensg_match[matched_ensg],
  ,
  drop = FALSE
]

ensg_table$seurat_gene <- ensg_genes[
  matched_ensg
]

ensg_table$match_type <- "gene_id"

# Gene symbol matches

symbol_table <- gene_annotation[
  symbol_match[matched_symbol],
  ,
  drop = FALSE
]

symbol_table$seurat_gene <- symbol_genes[
  matched_symbol
]

symbol_table$match_type <- "gene_name"

# ------------------------------------------------------------
# 16. COMBINE
# ------------------------------------------------------------

cnv_annotation <- rbind(
  ensg_table,
  symbol_table
)

# ------------------------------------------------------------
# 17. REMOVE DUPLICATE SEURAT GENES
# ------------------------------------------------------------

cnv_annotation <- cnv_annotation[
  !duplicated(cnv_annotation$seurat_gene),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 18. STANDARD CHROMOSOMES
# ------------------------------------------------------------

standard_chromosomes <- c(
  paste0("chr", 1:22),
  "chrX",
  "chrY"
)

cnv_annotation <- cnv_annotation[
  cnv_annotation$chromosome %in%
    standard_chromosomes,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 19. CHROMOSOME ORDER
# ------------------------------------------------------------

chromosome_order <- c(
  paste0("chr", 1:22),
  "chrX",
  "chrY"
)

cnv_annotation$chromosome <- factor(
  cnv_annotation$chromosome,
  levels = chromosome_order,
  ordered = TRUE
)

# ------------------------------------------------------------
# 20. SORT GENES BY GENOMIC POSITION
# ------------------------------------------------------------

cnv_annotation <- cnv_annotation[
  order(
    cnv_annotation$chromosome,
    cnv_annotation$start
  ),
  ,
  drop = FALSE
]

cnv_annotation$chromosome <- as.character(
  cnv_annotation$chromosome
)

# ------------------------------------------------------------
# 21. SUMMARY
# ------------------------------------------------------------

total_genes <- length(seurat_genes)

matched_genes <- length(
  unique(cnv_annotation$seurat_gene)
)

unmatched_genes <- total_genes - matched_genes

match_percentage <- (
  matched_genes /
    total_genes
) * 100

cat("\n")
cat("============================================================\n")
cat("MATCHING SUMMARY\n")
cat("============================================================\n\n")

cat(
  "Total Seurat genes:",
  total_genes,
  "\n"
)

cat(
  "Matched genes:",
  matched_genes,
  "\n"
)

cat(
  "Unmatched genes:",
  unmatched_genes,
  "\n"
)

cat(
  "Match percentage:",
  round(match_percentage, 2),
  "%\n\n"
)

# ------------------------------------------------------------
# 22. MATCH TYPE SUMMARY
# ------------------------------------------------------------

cat("Match type:\n\n")

print(
  table(cnv_annotation$match_type)
)

cat("\n")

# ------------------------------------------------------------
# 23. CHROMOSOME SUMMARY
# ------------------------------------------------------------

cat("Genes per chromosome:\n\n")

print(
  table(cnv_annotation$chromosome)
)

cat("\n")

# ------------------------------------------------------------
# 24. DUPLICATE CHECK
# ------------------------------------------------------------

duplicate_seurat_genes <- duplicated(
  cnv_annotation$seurat_gene
)

duplicate_gene_ids <- duplicated(
  cnv_annotation$gene_id
)

cat(
  "Duplicated Seurat genes:",
  sum(duplicate_seurat_genes),
  "\n"
)

cat(
  "Duplicated GENCODE gene IDs:",
  sum(duplicate_gene_ids),
  "\n\n"
)

# ------------------------------------------------------------
# 25. SAVE MAIN ANNOTATION
# ------------------------------------------------------------

output_csv <- file.path(
  output_dir,
  "GENCODE_v47_CNV_gene_annotation.csv"
)

write.csv(
  cnv_annotation,
  output_csv,
  row.names = FALSE
)

# ------------------------------------------------------------
# 26. SAVE MATCH TABLE
# ------------------------------------------------------------

match_table <- cnv_annotation[
  ,
  c(
    "seurat_gene",
    "gene_id",
    "gene_name",
    "gene_type",
    "chromosome",
    "start",
    "end",
    "match_type"
  )
]

match_table_file <- file.path(
  output_dir,
  "Seurat_genes_matched_to_GENCODE.csv"
)

write.csv(
  match_table,
  match_table_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 27. SAVE UNMATCHED GENES
# ------------------------------------------------------------

matched_gene_names <- unique(
  cnv_annotation$seurat_gene
)

unmatched <- seurat_genes[
  !seurat_genes %in%
    matched_gene_names
]

unmatched_file <- file.path(
  output_dir,
  "Seurat_unmatched_genes.csv"
)

write.csv(
  data.frame(
    seurat_gene = unmatched
  ),
  unmatched_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 28. SAVE SUMMARY FILE
# ------------------------------------------------------------

summary_file <- file.path(
  output_dir,
  "CNV_annotation_summary.txt"
)

summary_lines <- c(
  "CNV GENE ANNOTATION SUMMARY",
  "============================",
  "",
  "Reference: GENCODE v47",
  "Genome: GRCh38",
  "Ensembl: 113",
  "",
  paste(
    "Total Seurat genes:",
    total_genes
  ),
  paste(
    "Matched genes:",
    matched_genes
  ),
  paste(
    "Unmatched genes:",
    unmatched_genes
  ),
  paste(
    "Match percentage:",
    round(match_percentage, 2),
    "%"
  ),
  "",
  "Match type:",
  paste(
    capture.output(
      print(table(cnv_annotation$match_type))
    ),
    collapse = "\n"
  ),
  "",
  paste(
    "Standard chromosome genes:",
    nrow(cnv_annotation)
  ),
  ""
)

writeLines(
  summary_lines,
  summary_file
)

# ------------------------------------------------------------
# 29. SHOW EXAMPLES
# ------------------------------------------------------------

cat("First 10 rows:\n\n")

print(
  head(
    cnv_annotation,
    10
  )
)

cat("\n")

cat("Example symbol-based matches:\n\n")

print(
  head(
    cnv_annotation[
      cnv_annotation$match_type == "gene_name",
      ],
    10
  )
)

cat("\n")

cat("Example ENSG-based matches:\n\n")

print(
  head(
    cnv_annotation[
      cnv_annotation$match_type == "gene_id",
      ],
    10
  )
)

# ------------------------------------------------------------
# 30. FINAL
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("CNV ANNOTATION PREPARATION COMPLETED\n")
cat("============================================================\n\n")

cat("Saved files:\n\n")

cat(
  output_csv,
  "\n"
)

cat(
  match_table_file,
  "\n"
)

cat(
  unmatched_file,
  "\n"
)

cat(
  summary_file,
  "\n\n"
)

cat("Next step:\n")
cat(
  "Inspect the new matching percentage before CNV inference.\n"
)

cat("\n")
