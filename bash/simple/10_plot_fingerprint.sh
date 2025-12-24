#!/bin/bash
set -euo pipefail

mkdir -p deeptools_qc

plotFingerprint \
  -b encode_bam/ceb_ENCFF327JFG.bam encode_bam/H3K9ac_ENCFF193NPE.bam encode_bam/Input_ENCFF110SOB.bam \
  --skipZeros \
  --numberOfSamples 50000 \
  -T "Fingerprints of all BAM samples" \
  --plotFile deeptools_qc/fingerprints.pdf \
  --plotFileFormat pdf \
  --dpi 600 \
  2>&1 | tee deeptools_qc/plotFingerprint.log
