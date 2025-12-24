#!/bin/bash
set -euo pipefail

# Method 1: Step-by-step
# ls *.fastq.gz > samples.txt
# sed 's/.fastq.gz//' samples.txt > sample_id.txt

# Method 2: One-line command (Recommended)
# For Single-End:
# ls *.fastq.gz | sed 's/.fastq.gz//' > sample_id.txt

# For Paired-End (looking for _R1):
echo "Generating sample_id.txt from *_R1.fastq.gz files..."
ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//' > sample_id.txt

echo "Sample list created:"
cat sample_id.txt
