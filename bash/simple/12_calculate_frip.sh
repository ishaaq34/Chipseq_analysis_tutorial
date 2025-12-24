#!/bin/bash
set -euo pipefail

# Complete FRiP calculation for ceb_ENCFF327JFG
mkdir -p frip_analysis

# Count total reads
TOTAL_READS=$(samtools view -c encode_bam/ceb_ENCFF327JFG.bam)

# Merge overlapping peaks
sort -k1,1 -k2,2n macs3_results/ceb_ENCFF327JFG_peaks.narrowPeak | \
  bedtools merge -i stdin > frip_analysis/ceb_ENCFF327JFG.peaks.merged.bed

# Count reads in peaks
READS_IN_PEAKS=$(samtools view -b encode_bam/ceb_ENCFF327JFG.bam | \
  bedtools intersect -u -a stdin -b frip_analysis/ceb_ENCFF327JFG.peaks.merged.bed | \
  samtools view -c)

# Calculate FRiP
FRIP=$(awk -v n="${READS_IN_PEAKS}" -v d="${TOTAL_READS}" 'BEGIN {printf "%.5f", n/d}')

# Report results
echo "Sample: ceb_ENCFF327JFG"
echo "Total mapped reads : ${TOTAL_READS}"
echo "Reads in peaks     : ${READS_IN_PEAKS}"
echo "FRiP               : ${FRIP}"
