# Tutorial 09: Visualization (The Camera Angles)

##  1. Basic Concept (Camera Modes)

We want to visualize the "Average" pattern of our protein across all genes. To do this, we need to choose our **Camera Mode**:

1.  **Portrait Mode (Reference-Point):**
    *   **Focus:** One specific point (e.g., the Transcription Start Site, **TSS**).
    *   **Action:** We stand at the TSS and look 3kb upstream and 10kb downstream.
    *   **Use Case:** Great for seeing promoter activity (H3K9ac, Transcription Factors).

2.  **Panorama Mode (Scale-Regions):**
    *   **Focus:** The entire gene body.
    *   **Action:** Since genes are different lengths (short vs long), we stretch or compress them all to fit the same "frame" (e.g., 5000bp).
    *   **Use Case:** Great for seeing broad marks that cover the whole gene (H3K27me3, H3K36me3).

**The Process:**
1.  **Prepare the Map (BED):** Define where the genes are.
2.  **Create the Blueprint (Matrix):** `computeMatrix` calculates the numbers.
3.  **Take the Photo (Plot):** `plotHeatmap` or `plotProfile` draws the picture.

---

## 2. The Blueprint & The Photo - Basic requirement

### Step 1: Prep BED Files (The Map)
We need clean BED files defining the TSS and the Gene Body. We extract these from the annotation (GTF) file.

**A. Extract TSS (Portrait Mode)**
Run this `awk` script to get a 1bp position for every Transcript Start Site.
```bash
awk 'BEGIN{OFS="\t"} $3=="transcript" {                # keep only transcript rows  
  if($7=="+")                                         # plus strand: TSS at start
    print $1, $4-1, $4, $12, ".", $7;                 # [start-1, start) = 1 bp
  else                                                # minus strand: TSS at end
    print $1, $5-1, $5, $12, ".", $7                  # [end-1, end) = 1 bp
}' gencode.v49.annotation.gtf |  
tr -d '";' |                                          # remove quotes and semicolons
sort -k1,1V -k2,2n > tss.bed                          # sort by chrom, then start
```

**B. Extract Gene Bodies (Panorama Mode)**
Run this to get the full gene intervals.
```bash
awk 'BEGIN{OFS="\t"} $3=="gene" {                     # keep only gene rows
  print $1, $4-1, $5, $10, ".", $7                    # full gene body
}' gencode.v47.annotation.gtf | 
tr -d '";' |                                          # clean attributes
sort -k1,1V -k2,2n > genes.bed
```

### Build the Matrix (The Blueprint)
This command calculates the coverage scores for plotting.

**Mode A: Reference-Point (TSS)**
```bash
computeMatrix reference-point \
  --referencePoint TSS \
  -R "$REGIONS" \
  -S \
    "$BW/ceb_ENCFF327JFG.RPGC.bw" \
    "$BW/ceb_ENCFF744SVA.RPGC.bw" \
  -b 3000 -a 3000 \
  --binSize 250 \
  --numberOfProcessors 4 \
  -o bw_plot/ceb_TSS.mat.gz
```

**Mode B: Scale-Regions (Gene Body)**
```bash

computeMatrix scale-regions \
  -R "$REGIONS" \
  -S \
    "$BW/ceb_ENCFF327JFG.RPGC.bw" \
    "$BW/ceb_ENCFF744SVA.RPGC.bw" \
  --regionBodyLength 5000 \
  -b 3000 -a 3000 \
  --binSize 250 \
  --numberOfProcessors 4 \
  -o bw_plot/ceb_combined.mat.gz
```

### Step 3: Plotting (The Photo)

Now we turn the `matrix_TSS.gz` and `matrix_Body.gz` into plots.

```
plotProfile \
  -m ceb_TSS.mat.gz \
  --perGroup \
  -out ceb_TSS_profile.pdf
```

```
plotProfile \
  -m ceb_combined.mat.gz \
  --perGroup \
  -out ceb_combined_profile.pdf

````
likwise , you can geberate plots for other IPs and Inputs 

---

## Reading the Pictures


---

<img width="637" height="393" alt="image" src="https://github.com/user-attachments/assets/5c7b404c-af40-481a-b0db-9e16dbdfcbcf" />



---

Each panel shows average signal around transcription start sites (±3 kb), and the differences between tracks are clear.

The input tracks show low, smooth signal with no sharp peak at the TSS. This is what background looks like and confirms there is no real promoter-specific enrichment in the inputs.

CEB ChIP-seq shows weak enrichment relative to input and does not produce strong, sharp peaks. FRiP values are low (∼5–7%), which is typical for sequence-specific transcription factor ChIP-seq. Thus, individual CEBP binding events are difficult to discern at the whole-genome scale.

The H3K27me3 tracks show broader enrichment across the promoter region. The signal rises gradually toward the TSS and spreads over several kilobases, which fits a repressive chromatin mark that acts over domains rather than sharply at promoters.

The H3K9ac tracks are very different. They show strong, narrow peaks centered exactly at the TSS, with signal far higher than all other tracks. This reflects active promoter-associated acetylation and clear separation from background.

---
<img width="1022" height="323" alt="Screenshot 2025-12-10 at 12 14 14 PM" src="https://github.com/user-attachments/assets/0c10ed9f-6ff9-43a1-82d6-0c029973c56f" />

---

Across all genes, the meta-profiles show clear and reproducible mark-specific patterns that are obvious from both the shape of the curves and their amplitude on the Y axis. 

The **input tracks** sit around **~1.22–1.30** along the entire region, with only a very small rise near the TSS, which confirms that the background structure is low and there are no major mapping or copy-number artifacts. 

In contrast, the **ceb ChIP** reaches **~1.45–1.50** at the TSS before dropping back to ~1.30 across the gene body, a pattern and magnitude that match promoter-restricted occupancy of a transcription factor. 

**H3K27me3** shows a broader profile with TSS values rising to **~2.0–2.2** and gradually decaying toward ~1.5 at the TES, consistent with the expected spread of a repressive chromatin domain rather than a sharply localized promoter mark. 

**H3K9ac** displays the highest dynamic range: the TSS peak reaches ~8–10 in one replicate and ~5 in the other, dropping to ~2–3 across the gene body, which is exactly the kind of focused acetylation seen at active promoters. 

The fact that replicates for each mark reach similar Y-axis heights and have nearly identical curve shapes indicates that the enrichment is real, reproducible, and well above the input background.

---


### 3.2 Comparison: Portrait vs Panorama

| Mark | Uses TSS Mode (Portrait) | Uses Scaled Mode (Panorama) |
| :--- | :--- | :--- |
| **Input** | Flat (~1.25). | Flat (~1.25). |
| **CEBPA** | Small bump at TSS (~1.5). | Bump at start, then flat. |
| **H3K27me3** | Broad peak (~2.0-3.0). | Rise at start (~2.3), slope down gene body. |
| **H3K9ac** | Sharp spike (~25). | Spike at start, rapid drop to ~1. |
 

This step takes the per-replicate matrix files generated by computeMatrix and converts each one into a heatmap and plot using a single command. 

```
 plotHeatmap -m mat.gz \                 #`matrix_TSS.gz` or `matrix_Body.gz`
        -out profile_heat.pdf \
        --colorMap jet \
        --missingDataColor "#FFF6EB" \
        --refPointLabel "TSS" \
        --plotFileFormat pdf --dpi 600
```

---
<img width="852" height="436" alt="Screenshot 2025-12-10 at 12 47 38 PM" src="https://github.com/user-attachments/assets/57933355-b628-4b6f-a861-70bf54a60985" />

---


Input shows a flat, diffuse signal with no structured enrichment, consistent with background and technical bias.
H3K9ac displays a sharp, TSS-centered peak with strong promoter-specific enrichment, indicating a successful ChIP and correct biological localization. 





#### More 
Consider trying different parameters of the compute matrix and plotprofile/plotheatmap, and see how much they influence the plotting and representation. Consider having the bedfile of the lines/sines/pEl




