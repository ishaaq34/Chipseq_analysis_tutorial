#!/bin/bash
set -euo pipefail

if [ ! -f "sample_id.txt" ]; then
    echo "Error: sample_id.txt not found!"
    exit 1
fi

echo "Creating file lists from sample_id.txt..."

# All BAM files
BAM_FILES=$(while read sample; do echo "encode_bam/${sample}.bam"; done < sample_id.txt | tr '\n' ' ')

# Labels (basenames)
LABELS=$(cat sample_id.txt | tr '\n' ' ')

# IP samples only (excluding Input)
IP_FILES=$(grep -v "Input" sample_id.txt | while read sample; do echo "encode_bam/${sample}.bam"; done | tr '\n' ' ')
IP_LABELS=$(grep -v "Input" sample_id.txt | tr '\n' ' ')

echo "All MultiBamSummary Command with dynamic lists:"
echo "multiBamSummary bins -b $IP_FILES --labels $IP_LABELS -o deeptools_qc/matrix.npz"
