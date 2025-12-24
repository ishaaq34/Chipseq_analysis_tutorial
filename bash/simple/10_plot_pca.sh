#!/bin/bash
set -euo pipefail

mkdir -p deeptools_qc

plotPCA \
  -in deeptools_qc/matrix.npz \
  -o deeptools_qc/pca.pdf \
  -T "PCA of samples based on binned genome coverage" \
  --transpose \
  --plotWidth 14 \
  --plotHeight 12 \
  --plotFileFormat pdf \
  --dpi 600 \
  --outFileNameData deeptools_qc/pca.tab \
  --markers 's' 's' 'D' 'D' 'o' 'o' \
  --colors  '#1b9e77' '#66c2a5' \
            '#d95f02' '#fc8d62' \
            '#7570b3' '#8da0cb'
