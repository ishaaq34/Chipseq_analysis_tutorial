#!/bin/bash
set -euo pipefail

mkdir -p deeptools_qc

plotCoverage \
  -b encode_bam/ceb_ENCFF327JFG.bam encode_bam/H3K9ac_ENCFF193NPE.bam encode_bam/Input_ENCFF110SOB.bam \
  -o deeptools_qc/coverage_histogram.pdf \
  --plotFileFormat pdf \
  --dpi 600 \
  --smartLabels \
  --numberOfSamples 1000000 \
  --ignoreDuplicates \
  --minMappingQuality 30 \
  --outRawCounts deeptools_qc/coverage_counts.txt
