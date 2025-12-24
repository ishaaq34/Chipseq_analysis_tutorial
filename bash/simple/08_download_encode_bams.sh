#!/bin/bash
set -euo pipefail

mkdir -p encode_bam
cd encode_bam

# Download ENCODE BAM files
# Note: Requires wget. On macOS, install with 'brew install wget'
cat <<EOF | xargs -n 2 -P 4 wget -c -O
ceb_ENCFF327JFG.bam        https://www.encodeproject.org/files/ENCFF327JFG/@@download/ENCFF327JFG.bam
ceb_ENCFF744SVA.bam        https://www.encodeproject.org/files/ENCFF744SVA/@@download/ENCFF744SVA.bam
H3K27me3_ENCFF164ALR.bam   https://www.encodeproject.org/files/ENCFF164ALR/@@download/ENCFF164ALR.bam
H3K27me3_ENCFF532DQH.bam   https://www.encodeproject.org/files/ENCFF532DQH/@@download/ENCFF532DQH.bam
H3K9ac_ENCFF193NPE.bam     https://www.encodeproject.org/files/ENCFF193NPE/@@download/ENCFF193NPE.bam
H3K9ac_ENCFF534IPX.bam     https://www.encodeproject.org/files/ENCFF534IPX/@@download/ENCFF534IPX.bam
Input_ENCFF110SOB.bam      https://www.encodeproject.org/files/ENCFF110SOB/@@download/ENCFF110SOB.bam
Input_ENCFF919XCV.bam      https://www.encodeproject.org/files/ENCFF919XCV/@@download/ENCFF919XCV.bam
EOF

echo "Download complete."
