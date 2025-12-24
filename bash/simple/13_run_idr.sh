#!/bin/bash
set -euo pipefail

mkdir -p idr

# Run IDR on CEBPA replicates
idr --samples \
  macs3_results/ceb_ENCFF327JFG_peaks.narrowPeak \
  macs3_results/ceb_ENCFF744SVA_peaks.narrowPeak \
  --input-file-type narrowPeak \
  --rank signal.value \
  --output-file idr/ceb_idr_peaks.txt \
  --plot \
  --log-output-file idr/ceb_idr.log
