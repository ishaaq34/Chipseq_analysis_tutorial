#!/bin/bash
set -euo pipefail

RAW_DIR="fastq_raw"
mkdir -p "$RAW_DIR"

if [ ! -f "srr_list.txt" ]; then
    echo "Error: srr_list.txt not found!"
    exit 1
fi

while read -r acc; do
  echo "Downloading accession: $acc"

  fastq-dl \
    --accession "$acc" \
    --provider SRA \
    --cpus 1 \
    --outdir "$RAW_DIR"

  echo "Finished downloading: $acc"
done < srr_list.txt
