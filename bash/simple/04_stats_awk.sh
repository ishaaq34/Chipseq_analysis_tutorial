#!/bin/bash
set -euo pipefail

# Check for gunzip or zcat or gzcat
if command -v gunzip >/dev/null; then
    ZCAT="gunzip -c"
elif command -v zcat >/dev/null; then
    ZCAT="zcat"
elif command -v gzcat >/dev/null; then
    ZCAT="gzcat" # macOS
else
    echo "Error: No gzip decompressor found (gunzip, zcat, gzcat)."
    exit 1
fi

FILE="fastq_raw/H3K27me3_IP_rep1.fastq.gz"

if [ ! -f "$FILE" ]; then
    echo "Error: $FILE not found."
    exit 1
fi

echo "Counting reads for $FILE..."
$ZCAT "$FILE" | wc -l | awk '{print $1/4 " reads"}'

echo "Counting bases (coverage) for $FILE..."
$ZCAT "$FILE" | awk 'NR%4==2 {b+=length($0)} END{print b " bases"}'
