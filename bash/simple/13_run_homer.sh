#!/bin/bash
set -euo pipefail

# Filter peaks for IDR <= 0.05
awk '$12 >= 1.3 {print $1,$2,$3}' OFS="\t" \
  idr/ceb_idr_peaks.txt > idr/ceb_idr_passed.bed

# Run HOMER
findMotifsGenome.pl \
  idr/ceb_idr_passed.bed \
  GRCh38.primary_assembly.genome.fa \
  idr/cebpa_motifs/ \
  -size 200 \
  -mask \
  -p 8
