#!/bin/bash
set -euo pipefail

cd fastq_raw/ || { echo "Directory fastq_raw/ not found"; exit 1; }

# Rename SRR7297996 (H3K27me3 IP rep1)
if [ -f "SRR7297996.fastq.gz" ]; then
    mv SRR7297996.fastq.gz H3K27me3_IP_rep1.fastq.gz
fi

# Rename SRR7297997 (H3K27me3 IP rep2)
if [ -f "SRR7297997.fastq.gz" ]; then
    mv SRR7297997.fastq.gz H3K27me3_IP_rep2.fastq.gz
fi

# Rename SRR7298011 (Input rep1)
if [ -f "SRR7298011.fastq.gz" ]; then
    mv SRR7298011.fastq.gz Input_rep1.fastq.gz
fi

# Rename SRR7298012 (Input rep2)
if [ -f "SRR7298012.fastq.gz" ]; then
    mv SRR7298012.fastq.gz Input_rep2.fastq.gz
fi

echo "Renaming complete!"
