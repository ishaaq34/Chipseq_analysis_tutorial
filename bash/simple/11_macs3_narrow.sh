#!/bin/bash
set -euo pipefail

mkdir -p macs3_results

# H3K9ac Replicate 1 (Narrow Peak)
macs3 callpeak \
  -t encode_bam/H3K9ac_ENCFF534IPX.bam \
  -c encode_bam/Input_ENCFF110SOB.bam \
  -f BAM \
  -g hs \
  -n H3K9ac_ENCFF534IPX \
  -p 0.01 \
  --keep-dup all \
  --outdir macs3_results

# H3K9ac Replicate 2
macs3 callpeak \
  -t encode_bam/H3K9ac_ENCFF193NPE.bam \
  -c encode_bam/Input_ENCFF919XCV.bam \
  -f BAM \
  -g hs \
  -n H3K9ac_ENCFF193NPE \
  -p 0.01 \
  --keep-dup all \
  --outdir macs3_results

# CEBPA Replicates (Narrow Peak) can be added similarly
