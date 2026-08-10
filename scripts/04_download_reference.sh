#!/usr/bin/env bash

# ============================================================
# Project:
# Single-cell RNA-seq Analysis of Human Glioblastoma
#
# Dataset:
# GSE139448 (PRJNA579895)
#
# Script:
# 04_download_reference.sh
#
# Author:
# Katarzyna Zielinska
#
# Description:
# Download GRCh38 genome and GENCODE annotation.
# ============================================================

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

GENOME_DIR="${PROJECT_DIR}/data/reference/genome"
ANNOTATION_DIR="${PROJECT_DIR}/data/reference/annotation"

mkdir -p "${GENOME_DIR}"
mkdir -p "${ANNOTATION_DIR}"

echo "Downloading GRCh38 genome..."

wget -c \
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz \
-P "${GENOME_DIR}"

echo "Downloading GENCODE annotation..."

wget -c \
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz \
-P "${ANNOTATION_DIR}"

echo "Decompressing genome..."

gunzip -kf \
"${GENOME_DIR}/GRCh38.primary_assembly.genome.fa.gz"

echo "Decompressing annotation..."

gunzip -kf \
"${ANNOTATION_DIR}/gencode.v47.primary_assembly.annotation.gtf.gz"

echo
echo "Reference downloaded successfully."
