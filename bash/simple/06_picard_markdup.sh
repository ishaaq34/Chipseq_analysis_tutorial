#!/bin/bash
set -euo pipefail

mkdir -p picard_markdup picard_markdup_metrics

# Run Picard MarkDuplicates (Keeping duplicates for QC)
picard MarkDuplicates \
  I=picard_rg_bam/H3K27me3_IP_rep1.RG.bam \
  O=picard_markdup/H3K27me3_IP_rep1.marked.bam \
  M=picard_markdup_metrics/H3K27me3_IP_rep1.metrics.txt \
  REMOVE_DUPLICATES=false

# Index
samtools index picard_markdup/H3K27me3_IP_rep1.marked.bam
