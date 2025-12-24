#!/bin/bash
set -euo pipefail

mkdir -p spp_qc

# Locate run_spp.R
# If installed via conda, it might be in $CONDA_PREFIX/bin/run_spp.R or similar
RUN_SPP=$(which run_spp.R || echo "")

if [ -z "$RUN_SPP" ]; then
    # Fallback to checking typical location or asking user
    if [ -f "/opt/anaconda3/envs/chip/bin/run_spp.R" ]; then
        RUN_SPP="/opt/anaconda3/envs/chip/bin/run_spp.R"
    else
        echo "Error: run_spp.R not found in PATH. Please set path manually."
        exit 1
    fi
fi

echo "Using run_spp.R at: $RUN_SPP"

Rscript "$RUN_SPP" \
      -c=encode_bam/H3K9ac_ENCFF193NPE.bam \
      -savp=spp_qc/H3K9ac_ENCFF193NPE_spp.qc.pdf \
      -out=spp_qc/H3K9ac_ENCFF193NPE_spp.qc.txt
