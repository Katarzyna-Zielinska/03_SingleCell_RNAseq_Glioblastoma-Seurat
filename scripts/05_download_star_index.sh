#!/usr/bin/env bash

# ============================================================
# Project:
# Single-cell RNA-seq Analysis of Human Glioblastoma
#
# Dataset:
# GSE139448 (PRJNA579895)
#
# Script:
# 05_download_star_index.sh
#
# Author:
# Katarzyna Zielinska
#
# Description:
# Download pre-built STAR 2.7.11b genome index
# (GRCh38 Cell Ranger 2024-A compatible with STARsolo)
# ============================================================

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

DOWNLOAD_DIR="${PROJECT_DIR}/data/reference/downloads"
INDEX_DIR="${PROJECT_DIR}/data/reference/star_index"

mkdir -p "${DOWNLOAD_DIR}"
mkdir -p "${INDEX_DIR}"

ARCHIVE="GRCh38_star_2_7_11b.tar.gz"

URL="https://zenodo.org/records/11181586/files/${ARCHIVE}?download=1"

cd "${DOWNLOAD_DIR}"

echo "=========================================="
echo "Downloading STAR index..."
echo "=========================================="

wget -c -O "${ARCHIVE}" "${URL}"

echo
echo "=========================================="
echo "Extracting STAR index..."
echo "=========================================="

tar -xzf "${ARCHIVE}" -C "${INDEX_DIR}"

echo
echo "=========================================="
echo "Extraction completed."
echo "=========================================="

echo
echo "Files in STAR index directory:"

ls "${INDEX_DIR}"
