
mkdir -p bw_plot

# Calculate signals around TSS

computeMatrix reference-point \
    --referencePoint TSS \
    -b 4000 -a 4000 \
    -R tss.bed \
    -S a.ind_bw/*.bw \
    --skipZeros \
    -o bw_plot/matrix_all_bw_TSS.gz \
    --binSize 1000 \
    --numberOfProcessors 4

# Calculate signals across entire gene bodies

computeMatrix scale-regions \
    -R genes.bed \
    -S a.ind_bw/*.bw \
    --regionBodyLength 5000 \
    -b 1000 -a 1000 \
    --binSize 1000 \
    --skipZeros \
    -o bw_plot/matrix_all_bw_scalar.gz \
    --numberOfProcessors 4

plotHeatmap \
  -m bw_plot/matrix_all_bw_TSS.gz \
  -out bw_plot/heatmap_TSS.pdf \
  --colorMap jet \
  --missingDataColor "#FFF6EB" \
  --refPointLabel "TSS" \
  --dpi 600
  
plotHeatmap \
  -m bw_plot/matrix_all_bw_scalar.gz  \
  -out bw_plot/matrix_all_bw_scalar.pdf \
  --colorMap jet \
  --missingDataColor "#FFF6EB" \
  --refPointLabel "TSS" \
  --dpi 600
  
plotProfile -m bw_plot/matrix_all_bw_scalar.gz \
      --perGroup \
      --kmeans 2 \
      -plotType heatmap \
      -out bw_plot/profile_heatmap.pdf
