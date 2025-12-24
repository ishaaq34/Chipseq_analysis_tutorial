#!/bin/bash
set -euo pipefail

mkdir -p deeptools_qc

plotCorrelation \
  -in deeptools_qc/matrix.npz \
  --corMethod spearman \
  --skipZeros \
  --whatToPlot heatmap \
  --colorMap RdYlBu \
  --plotNumbers \
  --plotTitle "Spearman correlation of binned genome coverage" \
  --dpi 600 \
  -o deeptools_qc/spearman_corr_plot.pdf \
  --outFileCorMatrix deeptools_qc/spearman_corr_plot.tab
