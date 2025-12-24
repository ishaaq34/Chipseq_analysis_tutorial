#!/bin/bash
set -euo pipefail

echo "======================================"
echo "DRY RUN: BAM Deduplication and QC Pipeline"
echo "======================================"
echo "This will show what commands WOULD be executed"
echo ""

# Check if sample_id.txt exists
if [[ ! -f "sample_id.txt" ]]; then
  echo "ERROR: sample_id.txt not found!"
  exit 1
fi

echo "Step 1: Directory Creation"
echo "------------------------------------------------------"
echo "[DRY RUN] mkdir -p samtools_dedup_bam"
echo "[DRY RUN] mkdir -p samtools_dedup_metrics"
echo "[DRY RUN] mkdir -p deeptools_qc"
echo ""

echo "Step 2: BAM Deduplication Loop"
echo "------------------------------------------------------"

while read -r sample; do
  echo "Sample: $sample"
  echo ""
  
  # Show the piped command
  echo "[DRY RUN] samtools collate -u -O \"bowalign_filtered/${sample}.filtered.bam\" | \\"
  echo "            samtools fixmate -m -u - - | \\"
  echo "            samtools sort -u - | \\"
  echo "            samtools markdup -r - \"samtools_dedup_bam/${sample}.dedup.bam\""
  echo ""
  
  # Show indexing command
  echo "[DRY RUN] samtools index \"samtools_dedup_bam/${sample}.dedup.bam\""
  echo ""
  
  # Show flagstat command
  echo "[DRY RUN] samtools flagstat \"samtools_dedup_bam/${sample}.dedup.bam\" > \"samtools_dedup_metrics/${sample}.dedup.flagstat.txt\""
  echo ""
  
  echo "  Would create: samtools_dedup_bam/${sample}.dedup.bam"
  echo "  Would create: samtools_dedup_bam/${sample}.dedup.bam.bai"
  echo "  Would create: samtools_dedup_metrics/${sample}.dedup.flagstat.txt"
  echo ""
  echo "---"
  echo ""

done < sample_id.txt

echo "Step 3: Collecting BAM file paths"
echo "------------------------------------------------------"

# Simulate collecting BAM files
BAM_COUNT=$(wc -l < sample_id.txt)
echo "[DRY RUN] Would collect ${BAM_COUNT} deduplicated BAM files from samtools_dedup_bam/*.dedup.bam"
echo ""

echo "Example BAM list:"
while read -r sample; do
  echo "  - samtools_dedup_bam/${sample}.dedup.bam"
done < sample_id.txt
echo ""

echo "Step 4: plotFingerprint Command"
echo "------------------------------------------------------"
echo "[DRY RUN] plotFingerprint \\"
echo "  -b samtools_dedup_bam/*.dedup.bam \\"
echo "  --skipZeros \\"
echo "  --numberOfSamples 50000 \\"
echo "  -T \"Fingerprints of all BAM samples\" \\"
echo "  --plotFile deeptools_qc/fingerprints.pdf \\"
echo "  --plotFileFormat pdf \\"
echo "  --dpi 600 \\"
echo "  2>&1 | tee deeptools_qc/plotFingerprint.log"
echo ""
echo "  Would create: deeptools_qc/fingerprints.pdf"
echo "  Would create: deeptools_qc/plotFingerprint.log"
echo ""

echo "Step 5: plotCoverage Command"
echo "------------------------------------------------------"
echo "[DRY RUN] plotCoverage \\"
echo "  -b samtools_dedup_bam/*.dedup.bam \\"
echo "  -o deeptools_qc/coverage_histogram.pdf \\"
echo "  --plotFileFormat pdf \\"
echo "  --dpi 600 \\"
echo "  --smartLabels \\"
echo "  --numberOfSamples 1000000 \\"
echo "  --ignoreDuplicates \\"
echo "  --minMappingQuality 30 \\"
echo "  --outRawCounts deeptools_qc/coverage_counts.txt \\"
echo "  2>&1 | tee deeptools_qc/plotCoverage.log"
echo ""
echo "  Would create: deeptools_qc/coverage_histogram.pdf"
echo "  Would create: deeptools_qc/coverage_counts.txt"
echo "  Would create: deeptools_qc/plotCoverage.log"
echo ""

echo "======================================"
echo "DRY RUN Complete!"
echo "======================================"
echo ""
echo "Summary of what WOULD be created:"
echo ""
echo "Directories:"
echo "  - samtools_dedup_bam/"
echo "  - samtools_dedup_metrics/"
echo "  - deeptools_qc/"
echo ""
echo "Per-sample files (${BAM_COUNT} samples):"
echo "  - ${BAM_COUNT} × .dedup.bam files"
echo "  - ${BAM_COUNT} × .dedup.bam.bai index files"
echo "  - ${BAM_COUNT} × .dedup.flagstat.txt reports"
echo ""
echo "QC output files:"
echo "  - deeptools_qc/fingerprints.pdf"
echo "  - deeptools_qc/plotFingerprint.log"
echo "  - deeptools_qc/coverage_histogram.pdf"
echo "  - deeptools_qc/coverage_counts.txt"
echo "  - deeptools_qc/plotCoverage.log"
echo ""
echo "To run the actual pipeline, execute:"
echo "  ./bam_dedup_and_qc.sh"
echo ""
