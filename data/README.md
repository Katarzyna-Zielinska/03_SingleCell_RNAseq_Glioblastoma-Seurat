# Data

The raw sequencing data and reference files are not included in this repository because the files are too large for GitHub.

The single-cell RNA-seq data used in this project consist of three samples:

- SRR10353960
- SRR10353961
- SRR10353962

The original sequencing data can be downloaded from the NCBI Sequence Read Archive (SRA):

- SRR10353960
- SRR10353961
- SRR10353962

The project uses the following reference resources:

- Human GRCh38 reference genome
- GENCODE v47 gene annotation
- 10x Genomics 737K-august-2016 cell barcode whitelist
- STAR genome index generated from the GRCh38 reference genome and GENCODE v47 annotation

The raw FASTQ files, SRA files, reference genome, gene annotation, and STAR genome index are not stored in this repository.

The scripts required to download or prepare the reference resources and process the sequencing data are provided in the `scripts/` directory.
