#!/bin/bash
set -euo pipefail

mkdir -p QC_results

# Example for one sample.
# For loop version, wrap this in a while loop reading sample_id.txt

SAMPLE="H3K27me3_IP_rep1"
INPUT_BAM="picard_dedup_bam/${SAMPLE}.dedup.bam"

if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM $INPUT_BAM not found."
    exit 1
fi

echo "Calculating PBC metrics for $SAMPLE..."

# 1. Create a "5-prime" BED file
# (Convert BAM to simple coordinates, keeping only the start position of each read)
bedtools bamtobed -i "$INPUT_BAM" \
  | awk 'BEGIN{OFS="\t"} ($6=="+"){print $1,$2,$2+1} ($6=="-"){print $1,$3-1,$3}' \
  | sort -k1,1 -k2,2n \
  > "QC_results/${SAMPLE}.read5.bed"

# 2. Compute NRF, PBC1, PBC2
# (Count how many times each position appears)
uniq -c "QC_results/${SAMPLE}.read5.bed" \
  | awk '{c=$1; total+=c; uniq++; if(c==1) single++; if(c==2) double++;} \
    END{ if(total==0){print "NRF=NA\tPBC1=NA\tPBC2=NA"; exit} \
    NRF=uniq/total; \
    PBC1=single/uniq; \
    PBC2=(double? single/double:"Inf"); \
    printf("NRF=%.3f\tPBC1=%.3f\tPBC2=%s\n", NRF, PBC1, PBC2); }' \
  > "QC_results/${SAMPLE}.pbc.txt"

cat "QC_results/${SAMPLE}.pbc.txt"
