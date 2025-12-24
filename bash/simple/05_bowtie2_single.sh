#!/bin/bash
set -euo pipefail

mkdir -p bowalign

# Single-end example
bowtie2 -x genome_index/ce_index \
  -U fastq_cleaned/H3K27me3_IP_rep1.clean.fastq.gz \
  -p 6 --no-unal \
  2> bowalign/H3K27me3_IP_rep1.log | \
  samtools sort -@ 6 -o bowalign/H3K27me3_IP_rep1.sorted.bam

samtools index bowalign/H3K27me3_IP_rep1.sorted.bam

# Paired-end example (commented out)
# bowtie2 -x genome_index/ce_index \
#   -1 fastq_cleaned/Sample_R1.clean.fastq.gz \
#   -2 fastq_cleaned/Sample_R2.clean.fastq.gz \
#   -p 6 --no-unal \
#   2> bowalign/Sample.log | samtools sort -@ 6 -o bowalign/Sample.sorted.bam
# samtools index bowalign/Sample.sorted.bam
