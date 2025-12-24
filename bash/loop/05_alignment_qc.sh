#!/bin/bash
set -euo pipefail

mkdir -p bowalign_qc

while read -r sample; do
  echo "Running QC for: $sample"
  
  # Generate flagstat summary
  samtools flagstat "bowalign/${sample}.sorted.bam" > "bowalign_qc/${sample}.flagstat.txt"
  
  # Generate comprehensive stats (for MultiQC)
  samtools stats "bowalign/${sample}.sorted.bam" > "bowalign_qc/${sample}.stats.txt"
  
done < sample_id.txt

# Generate MultiQC report with all QC results
echo "Generating MultiQC report..."
if command -v multiqc &> /dev/null; then
    multiqc bowalign_qc/ -o bowalign_qc/ -n alignment_qc_report
    echo "MultiQC report saved: bowalign_qc/alignment_qc_report.html"
else
    echo "MultiQC not found, skipping report generation."
fi
