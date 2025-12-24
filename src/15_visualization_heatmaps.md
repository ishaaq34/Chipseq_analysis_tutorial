# Visualization (The Camera Angles)

`computeMatrix` `plotProfile` `plotHeatmap` `TSS-enrichment` `gene-body-profile` `visualization`

## 1. Basic Concept (Camera Modes)

We want to visualize the "Average" pattern of our protein across all genes. To do this, we need to choose our **Camera Mode**:

### 1. Portrait Mode (Reference-Point)

* **Focus:** One specific point (e.g., the Transcription Start Site, **TSS**).
* **Action:** We stand at the TSS and look 4kb upstream and 4kb downstream.
* **Use Case:** Great for seeing promoter activity (**H3K9ac**, Transcription Factors).

### 2. Panorama Mode (Scale-Regions)

* **Focus:** The entire gene body.
* **Action:** Since genes are different lengths (short vs long), we stretch or compress them all to fit the same "frame" (e.g., 5000bp).
* **Use Case:** Great for seeing broad marks that cover the whole gene (**H3K27me3**, H3K36me3).

**The Process:**

1. **Prepare the Map (BED):** Define where the genes are.
2. **Create the Blueprint (Matrix):** `computeMatrix` calculates the numbers.
3. **Take the Photo (Plot):** `plotHeatmap` or `plotProfile` draws the picture.

---

## 2. Requirements (The Files)

We need:

1. **BigWig Files:** Generated in Tutorial 12 (in `bigwigs/`).
2. **BED Files:** From Reference Annotation (Tutorial 01).
    * `tss.bed` (Start sites only)
    * `genes.bed` (Full gene bodies)

> [!NOTE]
> **Output Directory:** All results will go into `bw_plot/` to keep things organized.

---

## 3. Execution (The Blueprint)

We will calculate the matrices using two simple commands.

### Step 3.1: Portrait Mode (TSS Matrix)

This command calculates the signal ±4kb around the TSS.

```bash
mkdir -p bw_plot

# Calculate signals around TSS
computeMatrix reference-point \
    --referencePoint TSS \
    -b 4000 -a 4000 \
    -R tss.bed \
    -S bigwigs/*.bw \
    --skipZeros \
    -o bw_plot/matrix_all_bw_TSS.gz \
    --binSize 1000 \
    --numberOfProcessors 4
```

### Step 3.2: Panorama Mode (Gene Body Matrix)

This command stretches all genes to 5kb to compare them side-by-side.

```bash
# Calculate signals across entire gene bodies
computeMatrix scale-regions \
    -R genes.bed \
    -S bigwigs/*.bw \
    --regionBodyLength 5000 \
    -b 1000 -a 1000 \
    --binSize 1000 \
    --skipZeros \
    -o bw_plot/matrix_all_bw_scalar.gz \
    --numberOfProcessors 4
```

---

## 4. Plotting (The Photo)

Now that we have the "Blueprints" (matrices), we can develop the actual images.

### 4.1 Plot Profiles (Line Graphs)

This shows the *Average* signal across all genes.

```bash
# 1. TSS Profile (Portrait)
plotProfile \
  -m bw_plot/matrix_all_bw_TSS.gz \
  -out bw_plot/profile_TSS.pdf \
  --perGroup \
  --plotTitle "TSS Enrichment Profile" \
  --dpi 600

# 2. Gene Body Profile (Panorama)
plotProfile \
  -m bw_plot/matrix_all_bw_scalar.gz \
  -out bw_plot/profile_Body.pdf \
  --perGroup \
  --plotTitle "Gene Body Enrichment Profile" \
  --dpi 600
```

### 4.2 Plot Heatmaps (Rich Detail)

This shows the signal for *Every Single Gene* stacked on top of each other.

```bash
# Heatmap for TSS
plotHeatmap \
  -m bw_plot/matrix_all_bw_TSS.gz \
  -out bw_plot/heatmap_TSS.pdf \
  --colorMap jet \
  --missingDataColor "#FFF6EB" \
  --refPointLabel "TSS" \
  --dpi 600

# Heatmap for Gene Body
plotHeatmap \
  -m bw_plot/matrix_all_bw_scalar.gz \
  -out bw_plot/heatmap_genebody.pdf \
  --colorMap jet \
  --missingDataColor "#FFF6EB" \
  --refPointLabel "Gene Start" \
  --dpi 600
```

### 4.3 Advanced: K-means Clustering

This groups genes into clusters based on their signal patterns (useful for finding distinct gene classes).

```bash
plotProfile \
  -m bw_plot/matrix_all_bw_scalar.gz \
  --perGroup \
  --kmeans 2 \
  --plotType heatmap \
  -out bw_plot/profile_kmeans_heatmap.pdf
```

---

## 5. Reading the Results

### 5.1 The TSS Profile (Portrait)

This tells us about **Promoter Activity**.

<img width="900" height="271" alt="Screenshot 2025-12-15 at 7 10 31 PM" src="https://github.com/user-attachments/assets/294305c8-a488-45b8-bf90-d1c0a0d6a1be" />

* **Input (Background):** Flat line (~1.2-1.5). No enrichment. This is good quality control.
* **H3K9ac (Active Mark):** **Huge spike** right at the TSS (reaching ~12-25!). This confirms H3K9ac is strongly associated with active promoters.
* **H3K27me3 (Repressive Mark):** Broad, lower hill (~2.5-3.0). It doesn't spike at the TSS but sits broadly around it.
* **CEBPA (Factor):** Small but distinct bump at the TSS (~1.5). Transcription factors bind specific spots, so the average signal is lower than histone marks but still distinct from input.

### 5.2 The Gene Body Profile (Panorama)

This tells us about **Domain Structure**.

<img width="1022" height="323" alt="Screenshot 2025-12-10 at 12 14 14 PM" src="https://github.com/user-attachments/assets/0c10ed9f-6ff9-43a1-82d6-0c029973c56f" />

* **H3K9ac:** Spike at the start (TSS), then rapid drop. It's a "Starting Gun" mark.
* **H3K27me3:** Stays elevated across the **entire gene body**, slowly decaying. It's a "Blanket" mark that covers the whole region to silence it.

| Mark | Portrait (TSS) | Panorama (Gene Body) |
| :--- | :--- | :--- |
| **Input** | Flat. | Flat. |
| **H3K9ac** | **Sharp Spike** | Spike then Drop. |
| **H3K27me3** | Broad Hill | **Broad Blanket** (Entire Gene) |

> [!NOTE]
> **Up Next:** We have visualized the signals and confirmed they look biologically correct. Now, finally, we will use **MACS2** to mathematically identify the exact genomic coordinates of these peaks.
