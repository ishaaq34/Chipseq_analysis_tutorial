#!/bin/bash
set -euo pipefail

# Extract only chr11 and chr12 for testing
samtools view -b -h ceb_ENCFF327JFG.bam chr11 chr12 | samtools sort -o ceb_ENCFF327JFG.chr11_12.bam
samtools index ceb_ENCFF327JFG.chr11_12.bam
