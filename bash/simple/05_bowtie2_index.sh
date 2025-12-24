#!/bin/bash
set -euo pipefail

mkdir -p genome_index

# Assuming genome file exists
# bowtie2-build genome_index/ce.fa genome_index/ce_index
echo "Please download genome first. Example command:"
echo "bowtie2-build genome_index/ce.fa genome_index/ce_index"
