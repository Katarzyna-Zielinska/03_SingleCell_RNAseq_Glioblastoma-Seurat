#!/usr/bin/env bash

# ============================================================
# Project:
# Single-cell RNA-seq Analysis of Human Glioblastoma
#
# Dataset:
# GSE139448 (PRJNA579895)
#
# Script:
# 07_starsolo.sh
#
# Author:
# Katarzyna Zielinska
#
# Description:
# Quantification of single-cell RNA-seq data using STARsolo
# (STAR 2.7.11b, 10x Chromium v2).
# ============================================================

set -euo pipefail

##############################################################
# Directories
##############################################################

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

RAW_DIR="${PROJECT_DIR}/data/raw_fastq"

GENOME_DIR="${PROJECT_DIR}/data/reference/star_index"

WHITELIST="${PROJECT_DIR}/data/reference/whitelists/737K-august-2016.txt"

RESULTS_DIR="${PROJECT_DIR}/results/starsolo"

LOG_DIR="${PROJECT_DIR}/logs"

TMP_ROOT="${HOME}/starsolo_tmp"

THREADS=$(nproc)

SAMPLES=(
SRR10353960
SRR10353961
SRR10353962
)

mkdir -p "${RESULTS_DIR}"
mkdir -p "${LOG_DIR}"
mkdir -p "${TMP_ROOT}"

##############################################################
# Cleanup
##############################################################

cleanup(){

echo
echo "Interrupted."
echo "Temporary files kept in:"
echo "${TMP_ROOT}"

}

trap cleanup INT TERM

##############################################################
# Check requirements
##############################################################

check_requirements(){

echo "=========================================="
echo "Checking requirements..."
echo "=========================================="

command -v STAR >/dev/null || {
echo "STAR not found."
exit 1
}

command -v zcat >/dev/null || {
echo "zcat not found."
exit 1
}

[[ -f "${WHITELIST}" ]] || {
echo "Whitelist missing."
exit 1
}

[[ -f "${GENOME_DIR}/Genome" ]] || {
echo "STAR index missing."
exit 1
}

for SAMPLE in "${SAMPLES[@]}"
do

[[ -f "${RAW_DIR}/${SAMPLE}_1.fastq.gz" ]] || {
echo "${SAMPLE}_1.fastq.gz missing."
exit 1
}

[[ -f "${RAW_DIR}/${SAMPLE}_2.fastq.gz" ]] || {
echo "${SAMPLE}_2.fastq.gz missing."
exit 1
}

done

echo
echo "STAR version:"
STAR --version

echo
echo "Requirements OK."
echo

}

##############################################################
# Run one sample
##############################################################

run_sample(){

SAMPLE=$1

OUT_PREFIX="${RESULTS_DIR}/${SAMPLE}_"

TMP_DIR="${TMP_ROOT}/${SAMPLE}"

LOG_FILE="${LOG_DIR}/${SAMPLE}.log"

if [[ -d "${OUT_PREFIX}Solo.out/Gene" ]]
then

echo
echo "${SAMPLE} already completed."
echo "Skipping."

return

fi

rm -rf "${TMP_DIR}"

echo
echo "=========================================="
echo "Sample:"
echo "${SAMPLE}"
echo "Started:"
date
echo "=========================================="

START=$(date +%s)

STAR \
--runThreadN "${THREADS}" \
--genomeDir "${GENOME_DIR}" \
--outTmpDir "${TMP_DIR}" \
--readFilesIn \
"${RAW_DIR}/${SAMPLE}_2.fastq.gz" \
"${RAW_DIR}/${SAMPLE}_1.fastq.gz" \
--readFilesCommand zcat \
--soloType CB_UMI_Simple \
--soloCBstart 1 \
--soloCBlen 16 \
--soloUMIstart 17 \
--soloUMIlen 10 \
--soloCBwhitelist "${WHITELIST}" \
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
--soloUMIfiltering MultiGeneUMI_CR \
--soloUMIdedup 1MM_CR \
--soloCellFilter EmptyDrops_CR \
--soloFeatures Gene GeneFull Velocyto \
--soloMultiMappers EM \
--clipAdapterType CellRanger4 \
--outFilterScoreMin 30 \
--outSAMtype None \
--outFileNamePrefix "${OUT_PREFIX}" \
2>&1 | tee "${LOG_FILE}"

END=$(date +%s)

echo
echo "Finished:"
date

echo "Elapsed: $(( (END-START)/60 )) minutes"

rm -rf "${TMP_DIR}"

}

##############################################################
# Main
##############################################################

START_ALL=$(date +%s)

check_requirements

echo
echo "=========================================="
echo "Running STARsolo"
echo "=========================================="

for SAMPLE in "${SAMPLES[@]}"
do

run_sample "${SAMPLE}"

done

END_ALL=$(date +%s)

echo
echo "=========================================="
echo "STARsolo finished."
echo "=========================================="

echo
echo "Results:"
echo "${RESULTS_DIR}"

echo
echo "Logs:"
echo "${LOG_DIR}"

echo
echo "Total time:"
echo "$(( (END_ALL-START_ALL)/60 )) minutes"

echo
echo "Done."
