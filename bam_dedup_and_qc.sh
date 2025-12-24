#!/bin/bash
set -euo pipefail

# Configuration
THREADS=8  # Adjust based on your CPU cores

echo "======================================"
echo "BAM Deduplication and QC Pipeline"
echo "======================================"
echo "Using ${THREADS} threads"
echo ""

# Create output directories
mkdir -p samtools_dedup_bam
mkdir -p samtools_dedup_metrics
mkdir -p deeptools_qc

# Step 1: Deduplicate all BAMs and generate flagstat
echo ""
echo "Step 1: Deduplicating BAMs and generating flagstat..."
echo "------------------------------------------------------"

while read -r sample; do
  echo "Processing $sample..."

  # All steps piped together - no temp files needed!
  # Using THREADS for samtools commands
  samtools collate -@ ${THREADS} -u -O "bowalign_filtered/${sample}.filtered.bam" | \
    samtools fixmate -@ ${THREADS} -m -u - - | \
    samtools sort -@ ${THREADS} -u - | \
    samtools markdup -@ ${THREADS} -r - "samtools_dedup_bam/${sample}.dedup.bam"

  # Index the final output
  samtools index -@ ${THREADS} "samtools_dedup_bam/${sample}.dedup.bam"

  # Generate flagstat for deduplicated BAM
  samtools flagstat -@ ${THREADS} "samtools_dedup_bam/${sample}.dedup.bam" > "samtools_dedup_metrics/${sample}.dedup.flagstat.txt"
  
  echo "  ✓ Deduplicated: samtools_dedup_bam/${sample}.dedup.bam"
  echo "  ✓ Flagstat: samtools_dedup_metrics/${sample}.dedup.flagstat.txt"
  echo ""

done < sample_id.txt

echo "Deduplication complete!"
echo ""

# Step 2: Collect all deduplicated BAM files
echo "Step 2: Preparing for deepTools QC..."
echo "------------------------------------------------------"

# Create array of all deduplicated BAM files
DEDUP_BAMS=(samtools_dedup_bam/*.dedup.bam)

echo "Found ${#DEDUP_BAMS[@]} deduplicated BAM files"
echo ""

# Step 3: Run plotFingerprint
echo "Step 3: Running plotFingerprint..."
echo "------------------------------------------------------"

plotFingerprint \
  -b "${DEDUP_BAMS[@]}" \
  -p ${THREADS} \
  --skipZeros \
  --numberOfSamples 50000 \
  -T "Fingerprints of all BAM samples" \
  --plotFile deeptools_qc/fingerprints.pdf \
  --plotFileFormat pdf \
  --dpi 600 \
  2>&1 | tee deeptools_qc/plotFingerprint.log

echo "  ✓ Fingerprint plot: deeptools_qc/fingerprints.pdf"
echo ""

# Step 4: Run plotCoverage
echo "Step 4: Running plotCoverage..."
echo "------------------------------------------------------"

plotCoverage \
  -b "${DEDUP_BAMS[@]}" \
  -p ${THREADS} \
  -o deeptools_qc/coverage_histogram.pdf \
  --plotFileFormat pdf \
  --dpi 600 \
  --smartLabels \
  --numberOfSamples 1000000 \
  --ignoreDuplicates \
  --minMappingQuality 30 \
  --outRawCounts deeptools_qc/coverage_counts.txt \
  2>&1 | tee deeptools_qc/plotCoverage.log

echo "  ✓ Coverage histogram: deeptools_qc/coverage_histogram.pdf"
echo "  ✓ Coverage counts: deeptools_qc/coverage_counts.txt"
echo ""

# Summary
echo "======================================"
echo "Pipeline Complete!"
echo "======================================"
echo ""
echo "Outputs:"
echo "  - Deduplicated BAMs: samtools_dedup_bam/"
echo "  - Flagstat reports: samtools_dedup_metrics/"
echo "  - deepTools QC plots: deeptools_qc/"
echo ""
echo "Next steps:"
echo "  1. Review flagstat: cat samtools_dedup_metrics/*.flagstat.txt"
echo "  2. Open plots: open deeptools_qc/fingerprints.pdf"
echo "  3. Check coverage: cat deeptools_qc/coverage_counts.txt"
echo ""
