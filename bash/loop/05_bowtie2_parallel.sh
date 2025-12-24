#!/bin/bash
set -euo pipefail

mkdir -p bowalign bowalign_log

if ! command -v parallel &> /dev/null; then
    echo "Error: GNU parallel is not installed."
    exit 1
fi

parallel -j 2 '
  bowtie2 \
    -x genome_index/ce_index \
    -U fastq_cleaned/{}.clean.fastq.gz \
    -p 4 \
    --no-unal \
    2> bowalign_log/{}.log \
  | samtools sort \
      -@ 2 \
      -m 1G \
      -o bowalign/{}.sorted.bam

  samtools index bowalign/{}.sorted.bam
' :::: sample_id.txt
