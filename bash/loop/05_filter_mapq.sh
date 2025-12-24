#!/bin/bash
set -euo pipefail

mkdir -p bowalign_filtered

while read -r sample; do
  echo "Filtering multi-mappers for: $sample"
  samtools view -b -q 30 "bowalign/${sample}.sorted.bam" > "bowalign_filtered/${sample}.filtered.bam"
  samtools index "bowalign_filtered/${sample}.filtered.bam"
done < sample_id.txt
