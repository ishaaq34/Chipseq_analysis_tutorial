#!/bin/bash
set -euo pipefail

RAW_DIR="fastq_raw"
OUT_DIR="bowtie_align"

mkdir -p "$OUT_DIR"

if [ ! -f "sample_id.txt" ]; then
    echo "Error: sample_id.txt not found!"
    exit 1
fi

while read -r sample; do
  echo "sample_id: $sample"

  fq1="${RAW_DIR}/${sample}_R1.fastq.gz"
  fq2="${RAW_DIR}/${sample}_R2.fastq.gz"
  bam="${OUT_DIR}/${sample}.sorted.bam"

  echo "inputs:"
  echo "  $fq1"
  echo "  $fq2"

  echo "output:"
  echo "  $bam"

  echo "bowtie2 command: bowtie2 -x hg38_index -1 $fq1 -2 $fq2 | samtools sort -o $bam"
  echo
done < sample_id.txt
