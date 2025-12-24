#!/bin/bash
set -euo pipefail

mkdir -p samtools_dedup_bam

if [ ! -f "sample_id.txt" ]; then
    echo "Error: sample_id.txt not found!"
    exit 1
fi

while read -r sample; do
  echo "Processing $sample..."

  # All steps piped together - no temp files needed!
  # Input uses bowalign_filtered folder from previous step
  samtools collate -u -O "bowalign_filtered/${sample}.filtered.bam" | \
    samtools fixmate -m -u - - | \
    samtools sort -u - | \
    samtools markdup -r - "samtools_dedup_bam/${sample}.dedup.bam"

  # Index the final output
  samtools index "samtools_dedup_bam/${sample}.dedup.bam"

  echo "Finished $sample"
done < sample_id.txt
