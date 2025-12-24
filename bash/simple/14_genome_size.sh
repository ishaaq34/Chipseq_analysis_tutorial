#!/bin/bash
set -euo pipefail

# Method 1: faSize (General)
samtools faidx genome.fasta chr11 chr12 > chr11_chr12.fasta
faSize chr11_chr12.fasta > chr11_chr12.faSize.txt
awk '{nonN = $2 - $5; sum += nonN} END {print sum}' chr11_chr12.faSize.txt

# Method 2: khmer (Unique - Recommended)
# Requires khmer package
unique-kmers.py -k 21 chr11_chr12.fasta
