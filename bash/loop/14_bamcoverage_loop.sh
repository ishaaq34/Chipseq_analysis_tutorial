#!/bin/bash
set -euo pipefail

mkdir -p bigwigs

if [ ! -f "sample_id.txt" ]; then
    echo "Warning: sample_id.txt not found. Please create it first."
    exit 1
fi

# Loop through each sample ID in the text file
while read id; do
  
  echo "Generating BigWig for: $id"
  
  bamCoverage \
    -b "encode_bam/${id}.bam" \
    -o "bigwigs/${id}.bw" \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 220798375 \
    --smoothLength 30 \
    --numberOfProcessors 4 \
    --ignoreDuplicates
    
done < sample_id.txt
