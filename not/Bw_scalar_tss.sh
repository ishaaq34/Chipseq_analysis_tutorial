
#!/bin/bash
set -euo pipefail


mkdir -p bigwig_smoothlength


# Loop through each sample ID in the text file
cat sample_id.txt | while read id; do
  
  echo "Generating BigWig for: $id"
  
  bamCoverage \
    -b /Volumes/CrucialX6/bws/bam/${id}.bam \                    # Input BAM file (constructed from ID)
    -o bigwig_smoothlenght/${id}.RPGC.bw \                        # Output BigWig file
    --binSize 10 \                               # Resolution (buckets of 10bp)
    --normalizeUsing RPGC \                      # Normalization method (1x coverage)
    --effectiveGenomeSize 2701495761  \            # Mappable size of chr11+chr12
    --smoothLength 30 \                          # Smooth signal to reduce noise
    --numberOfProcessors 8 \                     # Use 4 CPUs for speed

done


BW="bigwig_smoothlength"
REGIONS_GENES="genes.bed"
REGIONS_TSS="TSS.bed"
THREADS=8
BIN_GENE=100
BIN_TSS=50

# Create separate output directories
mkdir -p results_scale_regions
mkdir -p results_tss

# CEBPA (two replicates)
computeMatrix scale-regions \
  -R "$REGIONS_GENES" \
  -S \
    "$BW/ceb_ENCFF327JFG.RPGC.bw" \
    "$BW/ceb_ENCFF744SVA.RPGC.bw" \
  --regionBodyLength 5000 \
  -b 3000 -a 3000 \
  --binSize "$BIN_GENE" \
  --numberOfProcessors "$THREADS" \
  -o results_scale_regions/ceb_combined.mat.gz

plotProfile \
  -m results_scale_regions/ceb_combined.mat.gz \
  --perGroup \
  -out results_scale_regions/ceb_combined_profile.pdf


# H3K27me3
computeMatrix scale-regions \
  -R "$REGIONS_GENES" \
  -S \
    "$BW/H3K27me3_ENCFF532DQH.RPGC.bw" \
    "$BW/H3K27me3_ENCFF164ALR.RPGC.bw" \
  --regionBodyLength 5000 \
  -b 3000 -a 3000 \
  --binSize "$BIN_GENE" \
  --numberOfProcessors "$THREADS" \
  -o results_scale_regions/H3K27me3_combined.mat.gz



plotProfile \
  -m results_scale_regions/H3K27me3_combined.mat.gz \
  --perGroup \
  -out results_scale_regions/H3K27me3_combined.pdf



# H3K9ac
computeMatrix scale-regions \
  -R "$REGIONS_GENES" \
  -S \
    "$BW/H3K9ac_ENCFF193NPE.RPGC.bw" \
    "$BW/H3K9ac_ENCFF534IPX.RPGC.bw" \
  --regionBodyLength 5000 \
  -b 3000 -a 3000 \
  --binSize "$BIN_GENE" \
  --numberOfProcessors "$THREADS" \
  -o results_scale_regions/H3K9ac_combined.mat.gz


plotProfile \
  -m results_scale_regions/H3K9ac_combined.mat.gz \
  --perGroup \
  -out results_scale_regions/H3K9ac_combined.pdf

# Input
computeMatrix scale-regions \
  -R "$REGIONS_GENES" \
  -S \
    "$BW/Input_ENCFF110SOB.RPGC.bw" \
    "$BW/Input_ENCFF919XCV.RPGC.bw" \
  --regionBodyLength 5000 \
  -b 3000 -a 3000 \
  --binSize "$BIN_GENE" \
  --numberOfProcessors "$THREADS" \
  -o results_scale_regions/Input_combined.mat.gz

  
plotProfile \
  -m results_scale_regions/Input_combined.mat.gz \
  --perGroup \
  -out results_scale_regions/Input_combined.pdf
  
##############################
# TSS Analysis (reference-point mode)
##############################

computeMatrix reference-point \
  --referencePoint TSS \
  -R "$REGIONS_TSS" \
  -S \
    "$BW/ceb_ENCFF327JFG.RPGC.bw" \
    "$BW/ceb_ENCFF744SVA.RPGC.bw" \
  -b 3000 -a 3000 \
  --binSize "$BIN_TSS" \
  --numberOfProcessors "$THREADS" \
  -o results_tss/ceb_TSS.mat.gz

plotProfile \
  -m results_tss/ceb_TSS.mat.gz \
  --perGroup \
  -out results_tss/ceb_TSS_profile.pdf



computeMatrix reference-point \
  --referencePoint TSS \
  -R "$REGIONS_TSS" \
  -S \
    "$BW/H3K27me3_ENCFF532DQH.RPGC.bw" \
    "$BW/H3K27me3_ENCFF164ALR.RPGC.bw" \
  -b 3000 -a 3000 \
  --binSize "$BIN_TSS" \
  --numberOfProcessors "$THREADS" \
  -o results_tss/H3K27me3_TSS.mat.gz

plotProfile \
  -m results_tss/H3K27me3_TSS.mat.gz \
  --perGroup \
  -out results_tss/H3K27me3_TSS_profile.pdf


computeMatrix reference-point \
  --referencePoint TSS \
  -R "$REGIONS_TSS" \
  -S \
    "$BW/H3K9ac_ENCFF193NPE.RPGC.bw" \
    "$BW/H3K9ac_ENCFF534IPX.RPGC.bw" \
  -b 3000 -a 3000 \
  --binSize "$BIN_TSS" \
  --numberOfProcessors "$THREADS" \
  -o results_tss/H3K9ac_TSS.mat.gz

plotProfile \
  -m results_tss/H3K9ac_TSS.mat.gz \
  --perGroup \
  -out results_tss/H3K9ac_TSS_profile.pdf


computeMatrix reference-point \
  --referencePoint TSS \
  -R "$REGIONS_TSS" \
  -S \
    "$BW/Input_ENCFF110SOB.RPGC.bw" \
    "$BW/Input_ENCFF919XCV.RPGC.bw" \
  -b 3000 -a 3000 \
  --binSize "$BIN_TSS" \
  --numberOfProcessors "$THREADS" \
  -o results_tss/Input_TSS.mat.gz

plotProfile \
  -m results_tss/Input_TSS.mat.gz \
  --perGroup \
  -out results_tss/Input_TSS_profile.pdf









  
