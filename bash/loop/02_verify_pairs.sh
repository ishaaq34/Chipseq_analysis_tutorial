#!/bin/bash
set -euo pipefail

# Checks if sample_id.txt exists
if [ ! -f "sample_id.txt" ]; then
    echo "Error: sample_id.txt not found!"
    exit 1
fi

while read -r sample; do
  echo "sample_id: $sample"

  fq1="${sample}_R1.fastq.gz"
  fq2="${sample}_R2.fastq.gz"

  echo "paired end: $fq1 : $fq2"
done < sample_id.txt
