# Image Migration Guide

**Date:** 2025-12-22  
**Purpose:** Migrate all tutorial images from GitHub assets to local storage

---

## Images Found (17 total)

### File 07 (1 image)

- Line 68: Library complexity plot

### File 09 (3 images)

- Line 22: Cross-correlation concept
- Line 87: Cross-correlation plot
- Line 106: Fragment length estimation

### File 10 (5 images)

- Line 202: Fingerprint plot
- Line 222: Coverage histogram A
- Line 231: Coverage histogram B (zoomed)
- Line 245: Correlation heatmap
- Line 260: PCA plot

### File 11a (1 image)

- Line 27: Peak types table (narrow/broad/exceptions)

### File 13 visualization (2 images)

- Line 156: Heatmap 1
- Line 167: Heatmap 2

### File 13 IDR (5 images)

- Line 142: Motif 1
- Line 152: Motif 2
- Line 159: Motif 3
- Line 168: IDR plot
- Line 264: HOMER results

---

## Migration Steps

### 1. Download Images

For each image URL, download and rename:

```bash
cd /Users/rajaishaqnabikhan/Desktop/codes_bioinformatics/chipseqtut_versions/Chipseq_analysis_tutorial/src/images

# File 11a
curl -L "https://github.com/user-attachments/assets/ea7d96e9-d8eb-4e4c-822c-0e7efe82821f" -o peak_types_narrow_broad.png

# File 07
curl -L "https://github.com/user-attachments/assets/730d3662-6973-4a4e-a9f9-fe565583ed32" -o library_complexity_plot.png

# File 09
curl -L "https://github.com/user-attachments/assets/cfc18aa5-2c0d-4772-ad2a-f4ae8c0f578d" -o cross_correlation_concept.png
curl -L "https://github.com/user-attachments/assets/696f703d-a182-466f-aeb3-bf19a07b0272" -o cross_correlation_plot.png
curl -L "https://github.com/user-attachments/assets/5b03189d-8fb7-4f4c-8145-5b5930921e05" -o fragment_length_estimate.png

# File 10
curl -L "https://github.com/user-attachments/assets/db1e7f33-6775-4705-b6b3-62a4c0bf7405" -o fingerprint_plot.png
curl -L "https://github.com/user-attachments/assets/26404eb5-7041-4674-bf94-d44d0d9edc8b" -o coverage_histogram_full.png
curl -L "https://github.com/user-attachments/assets/b8e24d14-d41b-4c55-a268-9982b49026c5" -o coverage_histogram_zoomed.png
curl -L "https://github.com/user-attachments/assets/652741bb-8b63-4fc2-88d8-eca5508e5938" -o correlation_heatmap.png
curl -L "https://github.com/user-attachments/assets/d2d785d1-4da3-4cc4-a3e0-1ee612ad3d38" -o pca_plot.png

# File 13 visualization
curl -L "https://github.com/user-attachments/assets/294305c8-a488-45b8-bf90-d1c0a0d6a1be" -o heatmap_viz_1.png
curl -L "https://github.com/user-attachments/assets/0c10ed9f-6ff9-43a1-82d6-0c029973c56f" -o heatmap_viz_2.png

# File 13 IDR/motifs
curl -L "https://github.com/user-attachments/assets/90e32a53-16c6-423e-b42f-af30f9f7ba7e" -o motif_1.png
curl -L "https://github.com/user-attachments/assets/26f276dc-2525-4775-b155-4ef31b4d01ae" -o motif_2.png
curl -L "https://github.com/user-attachments/assets/21224773-01fd-4d61-8b91-fb92ae71f0f5" -o motif_3.png
curl -L "https://github.com/user-attachments/assets/1e10e020-1138-46d2-b046-de1a656dd508" -o idr_plot.png
curl -L "https://github.com/user-attachments/assets/5a40f7b3-1a11-4cc7-8405-7001aafdff35" -o homer_results.png
```

### 2. Update Markdown References

Replace HTML img tags with markdown:

**File 11a line 27:**

```markdown
# Before
<img width="661" height="371" alt="Screenshot 2025-12-21 at 12 36 14 PM" src="https://github.com/user-attachments/assets/ea7d96e9-d8eb-4e4c-822c-0e7efe82821f" />

# After
![Peak types: narrow vs broad marks](./images/peak_types_narrow_broad.png)
```

Similar pattern for all 17 images.

### 3. Git Commit

```bash
cd /Users/rajaishaqnabikhan/Desktop/codes_bioinformatics/chipseqtut_versions/Chipseq_analysis_tutorial

# Add images
git add src/images/*.png

# Add updated markdown files
git add src/07_library_complexity.md
git add src/09_strand_cross_correlation.md
git add src/10_bam_summary_fingerprint.md
git add src/11a_macs3_peak_calling.md
git add src/13_visualization_heatmaps.md
git add src/13_idr_consensus_motifs_rk_corrected.md

# Commit
git commit -m "Migrate tutorial images to local storage

- Created src/images/ directory
- Downloaded 17 images from GitHub assets
- Updated all markdown files to use local image paths
- Improves portability and eliminates external dependencies"
```

---

## Verification

After migration, verify:

```bash
# Check all images exist
ls -lh src/images/

# Test markdown rendering
# Open files in markdown preview or mdbook build
```

---

## Benefits

✅ **Portable** - Works offline  
✅ **Version controlled** - Images tracked with code  
✅ **Reliable** - No external URL dependencies  
✅ **Standard** - Pure markdown, no HTML  
