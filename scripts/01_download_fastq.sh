#!/usr/bin/env bash

# ============================================================
# Project:
# Single-cell RNA-seq Analysis of Human Glioblastoma
#
# Dataset:
# GSE139448 (PRJNA579895)
#
# Script:
# 01_download_fastq.sh
#
# Author:
# Katarzyna Zielinska
#
# Description:
# Download raw FASTQ files from the SRA database.
# ============================================================

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_DIR="${PROJECT_DIR}/data/raw_fastq"

mkdir -p "${RAW_DIR}"
cd "${RAW_DIR}"

RUNS=(
    SRR10353960
    SRR10353961
    SRR10353962
)

for RUN in "${RUNS[@]}"; do
    echo "========================================"
    echo "Downloading ${RUN}"
    echo "========================================"

    if [[ ! -d "${RUN}" ]]; then
        prefetch "${RUN}"
    else
        echo "${RUN} already downloaded."
    fi

    if [[ ! -f "${RUN}_1.fastq.gz" || ! -f "${RUN}_2.fastq.gz" ]]; then
        fasterq-dump "${RUN}" --split-files --threads "$(nproc)"
        pigz -f -p "$(nproc)" "${RUN}_1.fastq"
        pigz -f -p "$(nproc)" "${RUN}_2.fastq"
    else
        echo "FASTQ files already exist."
    fi
done

echo
echo "Download completed successfully."
