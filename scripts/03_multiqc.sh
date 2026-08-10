#!/usr/bin/env bash

# ============================================================
# Project:
# Single-cell RNA-seq Analysis of Human Glioblastoma
#
# Dataset:
# GSE139448 (PRJNA579895)
#
# Script:
# 03_multiqc.sh
#
# Author:
# Katarzyna Zielinska
#
# Description:
# Aggregate FastQC reports into a single MultiQC report.
# ============================================================

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

FASTQC_DIR="${PROJECT_DIR}/results/fastqc"
OUT_DIR="${PROJECT_DIR}/results/multiqc"

mkdir -p "${OUT_DIR}"

echo "=========================================="
echo "Running MultiQC..."
echo "=========================================="

multiqc \
    "${FASTQC_DIR}" \
    --outdir "${OUT_DIR}" \
    --force

echo
echo "=========================================="
echo "MultiQC completed successfully."
echo "Results saved to:"
echo "${OUT_DIR}"
echo "=========================================="
