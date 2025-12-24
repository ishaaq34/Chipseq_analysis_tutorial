#!/bin/bash
set -euo pipefail

mkdir -p fastq_raw

if ! command -v parallel &> /dev/null; then
    echo "Error: GNU parallel is not installed."
    exit 1
fi

if [ ! -f "srr_list.txt" ]; then
    echo "Error: srr_list.txt not found!"
    exit 1
fi

parallel -j 4 \
  'echo "Starting download: {}" &&
   fastq-dl --accession {} --provider SRA --cpus 1 --outdir fastq_raw &&
   echo "Finished download: {}"' \
  :::: srr_list.txt
