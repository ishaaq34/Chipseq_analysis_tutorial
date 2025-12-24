#!/bin/bash
set -euo pipefail

mkdir -p picard_rg_bam

# Example for H3K27me3_IP_rep1
picard AddOrReplaceReadGroups \
  I=bowalign_filtered/H3K27me3_IP_rep1.filtered.bam \
  O=picard_rg_bam/H3K27me3_IP_rep1.RG.bam \
  RGID=H3K27me3_IP_rep1 \
  RGSM=H3K9ac

# Note: Adjust RGSM to your biological sample name
# RGID should be unique per library/replicate
