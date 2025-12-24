#!/bin/bash
set -euo pipefail

# Output format: count MAPQ_score
samtools view ceb_ENCFF327JFG.bam | awk '{print $5}' | sort -n | uniq -c
