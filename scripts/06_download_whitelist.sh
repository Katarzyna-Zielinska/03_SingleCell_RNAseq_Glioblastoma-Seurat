#!/usr/bin/env bash

# ============================================================
# Project:
# Single-cell RNA-seq Analysis of Human Glioblastoma
#
# Dataset:
# GSE139448 (PRJNA579895)
#
# Script:
# 06_download_whitelist.sh
#
# Author:
# Katarzyna Zielinska
#
# Description:
# Download the official 10x Genomics v2 barcode whitelist
# for STARsolo.
# ============================================================

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

WHITELIST_DIR="${PROJECT_DIR}/data/reference/whitelists"

mkdir -p "${WHITELIST_DIR}"

URL="https://raw.githubusercontent.com/10XGenomics/supernova/master/tenkit/lib/python/tenkit/barcodes/737K-august-2016.txt"

OUTPUT="${WHITELIST_DIR}/737K-august-2016.txt"

echo "=========================================="
echo "Downloading 10x Genomics v2 whitelist..."
echo "=========================================="

wget -c \
    -O "${OUTPUT}" \
    "${URL}"

echo
echo "=========================================="
echo "Download completed successfully."
echo "=========================================="

echo
echo "Whitelist location:"
echo "${OUTPUT}"

echo
echo "File information:"
ls -lh "${OUTPUT}"
