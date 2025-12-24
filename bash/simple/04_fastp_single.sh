#!/bin/bash
set -euo pipefail

# Ensure directories exist
mkdir -p fastq_cleaned

# Single-end cleaning
echo "Running fastp for single-end..."
fastp -i fastq_raw/H3K27me3_IP_rep1.fastq.gz -o fastq_cleaned/H3K27me3_IP_rep1.clean.fastq.gz

# Note for paired-end:
# fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz
