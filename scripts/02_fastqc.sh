#!/usr/bin/env bash

# ============================================================
# Project:
# Single-cell RNA-seq Analysis of Human Glioblastoma
#
# Dataset:
# GSE139448
#
# Script:
# 02_fastqc.sh
#
# Author:
# Katarzyna Zielinska
#
# Description:
# Quality control of raw FASTQ files using FastQC.
# ============================================================

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

RAW_DIR="${PROJECT_DIR}/data/raw_fastq"
RESULTS_DIR="${PROJECT_DIR}/results/fastqc"

mkdir -p "${RESULTS_DIR}"

echo "=========================================="
echo "Running FastQC..."
echo "=========================================="

fastqc \
    --threads "$(nproc)" \
    --outdir "${RESULTS_DIR}" \
    "${RAW_DIR}"/*.fastq.gz

echo
echo "=========================================="
echo "FastQC completed successfully."
echo "Results saved to:"
echo "${RESULTS_DIR}"
echo "=========================================="
