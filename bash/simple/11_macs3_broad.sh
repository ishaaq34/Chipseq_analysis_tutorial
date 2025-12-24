#!/bin/bash
set -euo pipefail

mkdir -p macs3_results

# H3K27me3 Replicate 1 (Broad Peak)
macs3 callpeak \
  -t encode_bam/H3K27me3_ENCFF532DQH.bam \
  -c encode_bam/Input_ENCFF110SOB.bam \
  -f BAM \
  -g hs \
  -n H3K27me3_ENCFF532DQH \
  --broad -p 0.01\
  --keep-dup all \
  --outdir macs3_results

# H3K27me3 Replicate 2
macs3 callpeak \
  -t encode_bam/H3K27me3_ENCFF164ALR.bam \
  -c encode_bam/Input_ENCFF919XCV.bam \
  -f BAM \
  -g hs \
  -n H3K27me3_ENCFF164ALR \
   --broad -p 0.01\
  --keep-dup all \
  --outdir macs3_results
