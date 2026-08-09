Single-Cell RNA-seq Analysis of Glioblastoma (GBM)

A complete single-cell RNA-seq (scRNA-seq) analysis workflow for glioblastoma (GBM), including raw sequencing data quality control, alignment and quantification with STARsolo, single-cell analysis and annotation using Seurat, and copy-number variation (CNV) analysis using inferCNV.

The project was designed as an end-to-end bioinformatics workflow, starting from raw FASTQ files and ending with an annotated single-cell dataset and quantitative CNV analysis.

Project Overview

This project presents a complete single-cell RNA-seq analysis pipeline performed on glioblastoma samples.

The workflow combines command-line bioinformatics tools with R/Bioconductor and Seurat for downstream single-cell analysis.

The main objectives were to:

perform quality control of raw sequencing data,
process single-cell RNA-seq reads,
generate gene-by-cell expression matrices,
perform single-cell quality control,
filter low-quality cells,
normalize the expression data,
perform dimensionality reduction,
identify cell clusters,
identify marker genes,
annotate major cell populations,
investigate sample and cell-type composition,
prepare genomic annotation for CNV analysis,
select reference cell populations,
perform CNV inference using inferCNV,
validate the resulting CNV data and observation cells.
Analysis Pipeline
                  Raw FASTQ files
                         │
                         ▼
                      FastQC
                         │
                         ▼
                     MultiQC
                         │
                         ▼
              Reference preparation
                         │
                         ▼
                      STARsolo
                         │
                         ▼
             Gene × Cell count matrix
                         │
                         ▼
                  Seurat object
                         │
                         ▼
                Single-cell QC
                         │
                         ▼
                  Cell filtering
                         │
                         ▼
                  Normalization
                         │
                         ▼
                       PCA
                         │
                         ▼
                  Clustering
                         │
                         ▼
                       UMAP
                         │
                         ▼
                  Marker analysis
                         │
                         ▼
              Cell-type annotation
                         │
                         ▼
             Sample composition
                         │
                         ▼
              ┌─────────────────────┐
              │      CNV analysis   │
              └─────────────────────┘
                         │
                         ▼
              GENCODE gene annotation
                         │
                         ▼
                CNV input preparation
                         │
                         ▼
              Reference population
                     validation
                         │
                         ▼
              Balanced reference set
                         │
                         ▼
                     inferCNV
                         │
                         ▼
              CNV result inspection
                         │
                         ▼
              Observation validation
Biological Question

What cellular populations are present in glioblastoma samples, how are they distributed across samples, and what large-scale copy-number variation patterns can be detected at the single-cell level?

Dataset
Information	Value
Organism	Homo sapiens
Disease	Glioblastoma (GBM)
Data type	Single-cell RNA-seq
Samples	3
Sample IDs	SRR10353960, SRR10353961, SRR10353962
Total cells	13,632
Reference genome	GRCh38
Gene annotation	GENCODE v47
Ensembl release	113

The dataset consists of three glioblastoma single-cell RNA-seq samples.

Bioinformatics Workflow
1. Raw Data Quality Control

Raw sequencing data were evaluated using:

FastQC
MultiQC

The quality-control stage was used to assess sequencing quality before downstream processing.

2. Reference Preparation

The workflow includes preparation of:

reference genome,
GENCODE gene annotation,
STAR genome index,
single-cell barcode whitelist.
3. Single-cell Alignment and Quantification

Single-cell reads were processed using STARsolo.

STARsolo was used to:

identify cell barcodes,
process UMI information,
align reads to the reference genome,
generate gene-level count matrices.

The resulting expression data were used as input for downstream Seurat analysis.

Seurat Analysis

The downstream single-cell analysis was performed using Seurat.

Quality Control

Quality-control metrics were evaluated for individual cells before filtering.

The workflow included:

number of detected genes,
number of RNA counts,
mitochondrial RNA percentage,
identification of low-quality cells.
Cell Filtering

Low-quality cells were removed according to the established QC criteria.

The resulting filtered dataset contained:

13,632 cells

Normalization

Gene expression data were normalized using the Seurat workflow.

The normalized expression matrix was subsequently used for dimensionality reduction and clustering.

Principal Component Analysis

Principal Component Analysis (PCA) was performed to reduce the dimensionality of the expression data and identify the major sources of transcriptional variation.

Clustering and UMAP

Cells were clustered according to their transcriptional profiles.

UMAP was used to visualize the resulting cellular structure in a low-dimensional space.

Marker Gene Analysis

Differentially expressed marker genes were identified for individual clusters.

Marker analysis was used to characterize the transcriptional profiles of the identified cell populations.

Cell-Type Annotation

Cell populations were annotated using marker-based analysis.

The final annotated dataset contained the following cell populations:

Astrocyte-like
Cycling-G2M
Cycling-Mixed
Cycling-S-phase
Endothelial
Highly-Cycling
Hypoxic-Mesenchymal-like
Mast-cell
Mesenchymal-like
Microglia-Macrophage
Neural-like
Neural-like-mixed
Oligodendrocyte
OPC-like
OPC-Neural-like
Vascular-Mesenchymal-like

The annotation therefore captures multiple major cellular compartments of the glioblastoma microenvironment, including neural, glial, vascular, immune, mesenchymal and cycling populations.

Sample Composition

The three samples contained:

Sample	Cells
SRR10353960	4,791
SRR10353961	4,652
SRR10353962	4,189
Total	13,632

Cell-type composition was examined across samples to evaluate differences in cellular representation.

CNV Analysis

A separate CNV analysis workflow was developed to investigate large-scale copy-number variation from single-cell RNA-seq expression data.

Because CNV inference depends strongly on genomic gene ordering and reference populations, the workflow includes several validation steps before inferCNV analysis.

GENCODE Gene Annotation

The Seurat gene identifiers were matched against:

GENCODE v47 — GRCh38

A total of:

38,606

genes were present in the Seurat object.

After matching against GENCODE:

37,344 genes

were successfully annotated.

The final matching rate was:

96.73%

The annotated genes included protein-coding genes, lncRNAs and additional gene types.

CNV Input Preparation

The CNV workflow generated:

normalized expression matrix,
raw count matrix,
genomic gene-order information,
cell metadata,
reference-cell metadata.

The final CNV input contained:

37,344 genes × 13,632 cells

The genes were ordered according to their genomic coordinates.

The gene-order validation confirmed correct ordering within chromosomes.

CNV Reference Populations

Potential non-malignant reference populations were evaluated before CNV inference.

The candidate populations included:

Astrocyte-like
Oligodendrocyte
Microglia-Macrophage
Endothelial
Mast-cell

For the final balanced reference set, two populations were selected:

Microglia-Macrophage
Oligodendrocyte

A total of 300 reference cells were selected:

Sample	Reference cells
SRR10353960	100
SRR10353961	100
SRR10353962	100

This ensured balanced reference representation across all three samples.

inferCNV Analysis

CNV inference was performed using inferCNV 1.18.1.

The final inferCNV analysis included:

13,632 cells

including:

300 reference cells
13,332 observation cells

The analysis was performed using genomic gene ordering based on GENCODE v47.

The inferCNV workflow included:

genomic gene ordering,
reference population definition,
sequencing-depth normalization,
expression transformation,
genomic smoothing,
CNV signal estimation,
observation-cell analysis,
CNV result inspection.
CNV Validation

The resulting inferCNV object was quantitatively inspected and validated.

The final inspected inferCNV object contained:

6,528 genes × 600 cells

with:

300 reference cells
300 observation cells

The validation confirmed representation of all three samples in the analyzed observation population.

The full observation-cell validation additionally confirmed that cells from all annotated populations were represented in the inferCNV analysis, although the proportions differed between cell types.

Results
UMAP

UMAP visualization demonstrates the transcriptional organization of the glioblastoma single-cell dataset and the separation of major cellular populations.

Cell-Type Annotation

The annotated UMAP highlights the cellular heterogeneity of the glioblastoma samples and identifies neural, glial, immune, vascular, mesenchymal and proliferative populations.

Sample Composition

The distribution of cell populations across the three samples demonstrates substantial differences in cellular composition between individual GBM samples.

CNV Heatmap

inferCNV analysis was used to visualize genomic CNV patterns across reference and observation cells.

The analysis provides a genome-wide view of relative CNV signal across ordered genomic regions.

Important Interpretation Note

The inferCNV analysis in this project was used to estimate and inspect CNV patterns in the single-cell dataset.

The project does not automatically classify individual cells as malignant or non-malignant based solely on the inferCNV output.

CNV patterns from scRNA-seq should be interpreted together with:

cell-type annotation,
marker gene expression,
sample context,
reference population quality,
genomic structure of the inferred signal.

Therefore, the current results should be considered a CNV inference and validation workflow, rather than a definitive malignant-cell classification.

Main Findings

The analysis demonstrates substantial cellular heterogeneity within the GBM dataset.

The final annotation identified multiple cellular populations, including:

Neural-like cells
Mesenchymal-like cells
Hypoxic-Mesenchymal-like cells
OPC-like cells
Microglia-Macrophage cells
Astrocyte-like cells
Oligodendrocytes
Endothelial cells
Mast cells
Cycling populations

The three samples also showed markedly different cell-type compositions.

The CNV workflow successfully integrated:

GENCODE v47 genomic annotation,
raw single-cell counts,
balanced reference populations,
sample-aware cell metadata,
inferCNV-based CNV estimation.
Technologies
R
RStudio
Seurat 5.5.1
SeuratObject 5.4.0
inferCNV 1.18.1
Bioconductor
STAR
STARsolo
FastQC
MultiQC
GENCODE
GRCh38
Linux
WSL2
Bash
Git & GitHub
Tools Used
Tool	Purpose
FastQC	Sequencing quality control
MultiQC	Aggregated QC reporting
STARsolo	Single-cell alignment and quantification
Seurat	Single-cell RNA-seq analysis
GENCODE	Gene annotation and genomic coordinates
inferCNV	Single-cell CNV inference
R	Statistical analysis and visualization
Bash	Workflow execution and data processing
Linux / WSL2	Computational environment
Repository Structure
03_SingleCell-RNAseq-Glioblastoma-Seurat/
│
├── data/
│   ├── counts/
│   ├── metadata/
│   ├── raw_fastq/
│   └── reference/
│
├── docs/
│
├── figures/
│
├── results/
│   ├── fastqc/
│   ├── multiqc/
│   ├── starsolo/
│   ├── seurat_qc/
│   ├── seurat_filter/
│   ├── seurat_normalization/
│   ├── seurat_pca/
│   ├── seurat_clustering/
│   ├── seurat_markers/
│   ├── seurat_marker_validation/
│   ├── seurat_annotation/
│   ├── sample_composition/
│   ├── cnv_annotation/
│   ├── cnv_input/
│   └── infercnv/
│
├── scripts/
│   ├── 01_download_fastq.sh
│   ├── 02_fastqc.sh
│   ├── 03_multiqc.sh
│   ├── 04_download_reference.sh
│   ├── 05_download_star_index.sh
│   ├── 06_download_whitelist.sh
│   ├── 07_starsolo.sh
│   ├── 08_seurat_qc.R
│   ├── 09_seurat_filter.R
│   ├── 10_seurat_normalization.R
│   ├── 11_seurat_pca.R
│   ├── 12_seurat_clustering_umap.R
│   ├── 13_seurat_markers.R
│   ├── 14_seurat_marker_validation.R
│   ├── 15_seurat_annotation.R
│   ├── 16_sample_composition.R
│   ├── 17_prepare_cnv_annotation.R
│   ├── 18_prepare_cnv_input.R
│   ├── 19_validate_cnv_input.R
│   ├── 20_prepare_cnv_counts.R
│   ├── 21_validate_cnv_reference.R
│   ├── 22_select_cnv_reference.R
│   ├── 23_validate_selected_cnv_reference.R
│   ├── 24_run_infercnv.R
│   ├── 25_inspect_infercnv_results.R
│   └── 26_validate_infercnv_observations.R
│
├── .gitignore
├── LICENSE
└── README.md
Reproducibility

The project is organized as a sequential workflow consisting of numbered scripts.

The scripts can be executed in order:

01 → 02 → 03 → ... → 26

The numbered structure makes it possible to follow the complete analysis from raw sequencing data through single-cell annotation and CNV inference.

Skills Demonstrated

This project demonstrates practical experience with:

Single-cell RNA-seq analysis
Raw sequencing data QC
FASTQ processing
STARsolo
Gene-by-cell count matrices
Seurat
Single-cell quality control
Data normalization
PCA
UMAP
Clustering
Marker gene analysis
Cell-type annotation
Sample composition analysis
Genomic annotation
GENCODE
GRCh38
CNV input preparation
Reference population selection
inferCNV
CNV validation
R programming
Bash scripting
Linux / WSL2
Reproducible bioinformatics workflows
Git & GitHub project organization
Future Improvements

Possible extensions of this project include:

More detailed malignant-cell identification using integrated CNV and expression evidence
Subclonal CNV analysis
Cell-level CNV scoring
Comparison of CNV patterns between annotated cellular populations
Integration of CNV profiles with GBM molecular subtypes
Additional validation using external genomic datasets
Automated workflow implementation using Snakemake or Nextflow
Containerization using Docker
Interactive visualization of single-cell CNV profiles
Project Status

Completed

The project currently includes:

raw-data QC,
STARsolo processing,
Seurat-based single-cell analysis,
cell-type annotation,
sample composition analysis,
CNV preparation and validation,
inferCNV analysis,
inferCNV result inspection,
observation-cell validation.

The final analysis stage is script 26.

Author

Katarzyna Zielińska
