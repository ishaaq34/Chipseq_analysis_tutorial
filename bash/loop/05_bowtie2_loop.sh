#!/bin/bash
set -euo pipefail

mkdir -p bowalign bowalign_log

if [ ! -f "sample_id.txt" ]; then
    echo "Error: sample_id.txt not found!"
    exit 1
fi

while read -r sample; do
  echo "Aligning $sample"

  # Single-end Loop
  bowtie2 -x genome_index/ce_index \
    -U "fastq_cleaned/${sample}.clean.fastq.gz" \
    -p 6 --no-unal \
    2> "bowalign_log/${sample}.bowtie2.log" \
    | samtools sort -@ 6 -o "bowalign/${sample}.sorted.bam"

  samtools index "bowalign/${sample}.sorted.bam"

  # Paired-end Loop (Uncomment to use)
  # bowtie2 -x genome_index/ce_index \
  #   -1 "fastq_cleaned/${sample}_R1.clean.fastq.gz" \
  #   -2 "fastq_cleaned/${sample}_R2.clean.fastq.gz" \
  #   -p 6 --no-unal \
  #   2> "bowalign_log/${sample}.log" | \
  # samtools sort -@ 6 -o "bowalign/${sample}.sorted.bam"
  # samtools index "bowalign/${sample}.sorted.bam"

done < sample_id.txt
