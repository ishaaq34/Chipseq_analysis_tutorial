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

CEB ChIP-seq and  input peaks are not much distinguisnbkle. Looking back to FRiP values are low (∼5–7%), which is typical for sequence-specific transcription factor ChIP-seq. Thus, individual CEBP binding events are difficult to discern at the whole-genome scale.

The H3K27me3 tracks show broader enrichment across the promoter region. The signal rises gradually toward the TSS and spreads over several kilobases, which fits a repressive chromatin mark that acts over domains rather than sharply at promoters.

The H3K9ac tracks are very different. They show strong, narrow peaks centered exactly at the TSS, with signal far higher than all other tracks. This reflects active promoter-associated acetylation and clear separation from background.

---
<img width="1022" height="323" alt="Screenshot 2025-12-10 at 12 14 14 PM" src="https://github.com/user-attachments/assets/0c10ed9f-6ff9-43a1-82d6-0c029973c56f" />

---

The genome-wide metagene profiles show clear differences between input, CEBP, and histone marks across transcription start sites (TSS) and transcription end sites (TES).

The input tracks show low-amplitude signal with mild structure around both TSS and TES, consistent with background chromatin features and technical biases rather than true enrichment. There is no sharp localization to either boundary.

CEB ChIP-seq and  input peaks are not much distinguisnbkle.

In contrast, H3K27me3 shows broad enrichment centered near the TSS that extends across the gene body and decays toward the TES. The signal is moderate in amplitude and spread over several kilobases, consistent with its role as a repressive chromatin mark acting over domains rather than forming sharp peaks.

H3K9ac displays strong, sharp enrichment precisely at the TSS, with signal levels far exceeding input and transcription factor profiles. The signal decreases steadily across the gene body and approaches baseline near the TES, consistent with promoter-focused acetylation associated with active transcription.

---

### Avergae bws 
## IP/Input
# Compute matrix and plots 



The CEBP profile shows a shallow negative dip centered at the TSS, with values remaining close to zero across the region. This indicates that, after normalization to input, CEBP enrichment at promoters is very weak and in places indistinguishable from background.

In contrast, H3K27me3 shows a modest but reproducible positive enrichment flanking the TSS, with a narrow dip exactly at the TSS itself.

H3K9ac displays a strong, clearly positive log2(IP/Input) peak centered at the TSS, reaching well above background and decaying symmetrically away from promoters.



Now, to find the peak of the Ceb peaks , we want to focus son peak that were identifdied in macs and idr produced cosnsesis peaks 


 








