#!/bin/bash
set -euo pipefail

mkdir -p picard_dedup_bam picard_dedup_metrics

# Run Picard to REMOVE duplicates
picard MarkDuplicates \
  I=picard_rg_bam/H3K27me3_IP_rep1.RG.bam \
  O=picard_dedup_bam/H3K27me3_IP_rep1.dedup.bam \
  M=picard_dedup_metrics/H3K27me3_IP_rep1.metrics.txt \
  REMOVE_DUPLICATES=true

# Index
samtools index picard_dedup_bam/H3K27me3_IP_rep1.dedup.bam
