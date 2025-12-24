#!/bin/bash
set -euo pipefail

mkdir -p deeptools_qc

# Excluding Input samples for correlation/PCA of IP enrichment
multiBamSummary bins \
  -b encode_bam/ceb_ENCFF327JFG.bam encode_bam/ceb_ENCFF744SVA.bam \
     encode_bam/H3K27me3_ENCFF164ALR.bam encode_bam/H3K27me3_ENCFF532DQH.bam \
     encode_bam/H3K9ac_ENCFF193NPE.bam encode_bam/H3K9ac_ENCFF534IPX.bam \
  --numberOfProcessors 4 \
  -o deeptools_qc/matrix.npz \
  --outRawCounts deeptools_qc/matrix.tab
