#!/bin/bash
set -euo pipefail

# BigWig Averaging Script
# Averages IP and Input replicates using bigwigCompare --operation mean

BW="bigwig_smoothlength"
mkdir -p averaged_bigwigs

echo "Averaging BigWig replicates..."

# Input replicates
echo "Averaging Input replicates..."
bigwigCompare \
  -b1 "$BW/Input_ENCFF110SOB.RPGC.bw" \
  -b2 "$BW/Input_ENCFF919XCV.RPGC.bw" \
  --operation mean \
  -o averaged_bigwigs/Input_mean.RPGC.bw

# H3K9ac replicates
echo "Averaging H3K9ac replicates..."
bigwigCompare \
  -b1 "$BW/H3K9ac_ENCFF193NPE.RPGC.bw" \
  -b2 "$BW/H3K9ac_ENCFF534IPX.RPGC.bw" \
  --operation mean \
  -o averaged_bigwigs/H3K9ac_mean.RPGC.bw

# H3K27me3 replicates
echo "Averaging H3K27me3 replicates..."
bigwigCompare \
  -b1 "$BW/H3K27me3_ENCFF532DQH.RPGC.bw" \
  -b2 "$BW/H3K27me3_ENCFF164ALR.RPGC.bw" \
  --operation mean \
  -o averaged_bigwigs/H3K27me3_mean.RPGC.bw

# CEBPA replicates
echo "Averaging CEBPA replicates..."
bigwigCompare \
  -b1 "$BW/ceb_ENCFF327JFG.RPGC.bw" \
  -b2 "$BW/ceb_ENCFF744SVA.RPGC.bw" \
  --operation mean \
  -o averaged_bigwigs/ceb_mean.RPGC.bw

echo "✅ All replicates averaged successfully!"
echo "Output directory: averaged_bigwigs/"
echo ""
echo "Generated files:"
ls -lh averaged_bigwigs/
