#!/bin/bash

# Image Migration Script
# Downloads all GitHub asset images to local storage

set -e

echo "Starting image migration..."

# Create images directory
mkdir -p src/images
cd src/images

echo "Downloading images..."

# File 11a
echo "  [1/17] peak_types_narrow_broad.png"
curl -sL "https://github.com/user-attachments/assets/ea7d96e9-d8eb-4e4c-822c-0e7efe82821f" -o peak_types_narrow_broad.png

# File 07
echo "  [2/17] library_complexity_plot.png"
curl -sL "https://github.com/user-attachments/assets/730d3662-6973-4a4e-a9f9-fe565583ed32" -o library_complexity_plot.png

# File 09
echo "  [3/17] cross_correlation_concept.png"
curl -sL "https://github.com/user-attachments/assets/cfc18aa5-2c0d-4772-ad2a-f4ae8c0f578d" -o cross_correlation_concept.png

echo "  [4/17] cross_correlation_plot.png"
curl -sL "https://github.com/user-attachments/assets/696f703d-a182-466f-aeb3-bf19a07b0272" -o cross_correlation_plot.png

echo "  [5/17] fragment_length_estimate.png"
curl -sL "https://github.com/user-attachments/assets/5b03189d-8fb7-4f4c-8145-5b5930921e05" -o fragment_length_estimate.png

# File 10
echo "  [6/17] fingerprint_plot.png"
curl -sL "https://github.com/user-attachments/assets/db1e7f33-6775-4705-b6b3-62a4c0bf7405" -o fingerprint_plot.png

echo "  [7/17] coverage_histogram_full.png"
curl -sL "https://github.com/user-attachments/assets/26404eb5-7041-4674-bf94-d44d0d9edc8b" -o coverage_histogram_full.png

echo "  [8/17] coverage_histogram_zoomed.png"
curl -sL "https://github.com/user-attachments/assets/b8e24d14-d41b-4c55-a268-9982b49026c5" -o coverage_histogram_zoomed.png

echo "  [9/17] correlation_heatmap.png"
curl -sL "https://github.com/user-attachments/assets/652741bb-8b63-4fc2-88d8-eca5508e5938" -o correlation_heatmap.png

echo "  [10/17] pca_plot.png"
curl -sL "https://github.com/user-attachments/assets/d2d785d1-4da3-4cc4-a3e0-1ee612ad3d38" -o pca_plot.png

# File 13 visualization
echo "  [11/17] heatmap_viz_1.png"
curl -sL "https://github.com/user-attachments/assets/294305c8-a488-45b8-bf90-d1c0a0d6a1be" -o heatmap_viz_1.png

echo "  [12/17] heatmap_viz_2.png"
curl -sL "https://github.com/user-attachments/assets/0c10ed9f-6ff9-43a1-82d6-0c029973c56f" -o heatmap_viz_2.png

# File 13 IDR/motifs
echo "  [13/17] motif_1.png"
curl -sL "https://github.com/user-attachments/assets/90e32a53-16c6-423e-b42f-af30f9f7ba7e" -o motif_1.png

echo "  [14/17] motif_2.png"
curl -sL "https://github.com/user-attachments/assets/26f276dc-2525-4775-b155-4ef31b4d01ae" -o motif_2.png

echo "  [15/17] motif_3.png"
curl -sL "https://github.com/user-attachments/assets/21224773-01fd-4d61-8b91-fb92ae71f0f5" -o motif_3.png

echo "  [16/17] idr_plot.png"
curl -sL "https://github.com/user-attachments/assets/1e10e020-1138-46d2-b046-de1a656dd508" -o idr_plot.png

echo "  [17/17] homer_results.png"
curl -sL "https://github.com/user-attachments/assets/5a40f7b3-1a11-4cc7-8405-7001aafdff35" -o homer_results.png

cd ../..

echo ""
echo "✅ Downloaded 17 images to src/images/"
echo ""
echo "Verifying downloads..."
ls -lh src/images/*.png | wc -l
echo ""
echo "Next steps:"
echo "1. Run update_image_references.sh to update markdown files"
echo "2. Commit changes: git add src/images/ src/*.md && git commit -m 'Migrate images to local storage'"
