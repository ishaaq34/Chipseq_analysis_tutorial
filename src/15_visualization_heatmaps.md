# Visualization with deepTools (Creating Signal Heatmaps)

## 1. Basic Concept (Camera Modes)

We want to visualize the "Average" pattern of our protein across all genes. To do this, we need to choose our **Camera Mode**:

1. **Portrait Mode (Reference-Point):**
    * **Focus:** One specific point (e.g., the Transcription Start Site, **TSS**).
    * **Action:** We stand at the TSS and look 3kb upstream and 10kb downstream.
    * **Use Case:** Great for seeing promoter activity (H3K9ac, Transcription Factors).

2. **Panorama Mode (Scale-Regions):**
    * **Focus:** The entire gene body.
    * **Action:** Since genes are different lengths (short vs long), we stretch or compress them all to fit the same "frame" (e.g., 5000bp).
    * **Use Case:** Great for seeing broad marks that cover the whole gene (H3K27me3, H3K36me3).

**The Process:**

1. **Prepare the Map (BED):** Define where the genes are.
2. **Create the Blueprint (Matrix):** `computeMatrix` calculates the numbers.
3. **Take the Photo (Plot):** `plotHeatmap` or `plotProfile` draws the picture.

---

## 2. The Blueprint & The Photo - Basic requirement

### Step 1: Download Reference Annotation

We need the GENCODE human genome annotation to define TSS and gene body regions.

**Download the GENCODE v49 GTF file:**

```bash
# Download GENCODE basic annotation (primary assembly)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.basic.annotation.gtf.gz

# Decompress
gunzip gencode.v49.primary_assembly.basic.annotation.gtf.gz
```

This GTF file contains comprehensive gene annotations including transcript start sites (TSS), gene bodies, and other genomic features.

### Step 2: Extract TSS Regions (Portrait Mode)

Run this `awk` script to get a 1bp position for every Transcript Start Site from the GTF:

```bash
awk 'BEGIN{OFS="\t"} $3=="transcript" {                 
  if($7=="+")                                         
    print $1, $4-1, $4, $12, ".", $7;                
  else                                               
    print $1, $5-1, $5, $12, ".", $7                 
}' gencode.v49.annotation.gtf |  
tr -d '";' |                                        
sort -k1,1V -k2,2n > tss.bed                        
```

This command extracts the **transcription start site (TSS)** for each transcript from the GTF annotation file. The script identifies whether each transcript is on the forward (+) or reverse (-) strand: for forward-strand transcripts, the TSS is the start coordinate (column 4), while for reverse-strand transcripts, it's the end coordinate (column 5).
. The output `tss.bed` is a BED file containing chromosome, TSS position (as a 1bp interval), transcript name, and strand information for every annotated transcript.

### Step 3: Extract Gene Bodies (Panorama Mode)

Run this to get the full gene intervals from the same GTF file:

```bash
awk 'BEGIN{OFS="\t"} $3=="gene" {                     
  print $1, $4-1, $5, $10, ".", $7                    
}' gencode.v49.annotation.gtf | 
tr -d '";' |                                          
sort -k1,1V -k2,2n > genes.bed
```

### Step 4: Create Output Directories

```bash
# Create directories for matrices and plots
mkdir -p deeptools_viz/matrices
mkdir -p deeptools_viz/plots
```

### Step 5: Build the Matrix (The Blueprint)

This command calculates the coverage scores for plotting.

**Mode A: Reference-Point (TSS)**

```bash
computeMatrix reference-point \
  --referencePoint TSS \
  -R tss.bed \
  -S \
    bigwigs/H3K9ac_ENCFF534IPX.bw \
    bigwigs/H3K9ac_ENCFF193NPE.bw \
  -b 3000 -a 3000 \
  --binSize 250 \
  --numberOfProcessors 4 \
  -o deeptools_viz/matrices/H3K9ac_TSS.mat.gz
```

**Mode B: Scale-Regions (Gene Body)**

```bash
computeMatrix scale-regions \
  -R genes.bed \
  -S \
    bigwigs/H3K9ac_ENCFF534IPX.bw \
    bigwigs/H3K9ac_ENCFF193NPE.bw \
  --regionBodyLength 5000 \
  -b 3000 -a 3000 \
  --binSize 250 \
  --numberOfProcessors 4 \
  -o deeptools_viz/matrices/H3K9ac_genes.mat.gz
```

### Step 6: Plotting (The Photo)

Now we turn the matrices into plots.

```bash
plotProfile \
  -m deeptools_viz/matrices/H3K9ac_TSS.mat.gz \
  --perGroup \
  -out deeptools_viz/plots/H3K9ac_TSS_profile.pdf
```

**Generating plots for all samples:**

Repeat the above `computeMatrix` and `plotProfile` commands for all your samples (Input, H3K27me3, CEBPA) to create comprehensive visualization of ChIP-seq signal patterns.

---

## Reading the Pictures

### TSS-Centered Profiles: Individual Replicates

The following plots show average ChIP-seq signal around transcription start sites (TSS ±3 kb) for individual replicates of Input, H3K27me3, and H3K9ac:

---

<img alt="image" src="./images/bw1.png"/>

---

Each panel shows average signal around transcription start sites (±3 kb), and the differences between tracks are clear.

The input tracks show low, smooth signal with no sharp peak at the TSS. This is what background looks like and confirms there is no real promoter-specific enrichment in the inputs.

The H3K27me3 tracks show broader enrichment across the promoter region. The signal rises gradually toward the TSS and spreads over several kilobases, which fits a repressive chromatin mark that acts over domains rather than sharply at promoters.

The H3K9ac tracks are very different. They show strong, narrow peaks centered exactly at the TSS, with signal far higher than all other tracks. This reflects active promoter-associated acetylation and clear separation from background.

---

### TSS/TES Metagene Profiles: Genome-Wide Patterns

The following plots show average ChIP-seq signal across transcription start sites (TSS) and transcription end sites (TES):

<img alt="Screenshot 2025-12-10 at 12 14 14 PM" src="./images/bw2.png" />

---

The genome-wide metagene profiles show clear differences between input, CEBP, and histone marks across transcription start sites (TSS) and transcription end sites (TES).

The input tracks show low-amplitude signal with mild structure around both TSS and TES, consistent with background chromatin features and technical biases rather than true enrichment. There is no sharp localization to either boundary.

In contrast, H3K27me3 shows broad enrichment centered near the TSS that extends across the gene body and decays toward the TES. The signal is moderate in amplitude and spread over several kilobases, consistent with its role as a repressive chromatin mark acting over domains rather than forming sharp peaks.

H3K9ac displays strong, sharp enrichment precisely at the TSS, with signal levels far exceeding input and transcription factor profiles. The signal decreases steadily across the gene body and approaches baseline near the TES, consistent with promoter-focused acetylation associated with active transcription.

---

## Average Signal Analysis

These are individual replicates of the different IPs. Now let's see how they look cumulatively.

### Average BigWigs

```
bigwigAverage -b H3K9ac_ENCFF193NPE.RPGC.bw H3K9ac_ENCFF534IPX.RPGC.bw \
  -o bw_mean/H3K9ac_mean.bw -p 6


bigwigAverage -b H3K27me3_ENCFF164ALR.RPGC.bw H3K27me3_ENCFF532DQH.RPGC.bw \
  -o bw_mean/H3K27me3_mean.bw -p $THREADS




```

## After normalisation to Inpouts , how these plots look like

```
bigwigCompare -b1 bw_mean/H3K9ac_mean.bw -b2 bw_mean/Input_mean.bw \
  --operation log2 --pseudocount 1 -p $THREADS \
  -o bw_log2/H3K9ac_log2IPoverInput.bw

bigwigCompare -b1 bw_mean/H3K27me3_mean.bw -b2 bw_mean/Input_mean.bw \
  --operation log2 --pseudocount 1 -p $THREADS \
  -o bw_log2/H3K27me3_log2IPoverInput.bw

```

# Compute matrix and plots TSS and Scalar plot (for H3k9ac)

```
computeMatrix reference-point \
    --referencePoint TSS \
    -b 3000 -a 3000 \
    --binSize 250 \
    -R $TSS_BED \
    -S bw_log2/${MARK}_log2IPoverInput.bw \
    --skipZeros --missingDataAsZero \
    -p $THREADS \
    -o matrix/${MARK}_TSS_log2.gz


plotProfile \
    -m matrix/${MARK}_TSS_log2.gz \
    --refPointLabel TSS \
    --yAxisLabel "log2(IP/Input)" \
    --plotTitle "${MARK} log2(IP/Input)" \
    -out plots/${MARK}_TSS_log2_profile.pdf

```

```
computeMatrix scale-regions \
    -b 3000 -a 3000 \
    --regionBodyLength 5000 \
    --binSize 250 \
    -R $GENES_BED \
    -S bw_mean/${MARK}_mean.bw \
    --skipZeros --missingDataAsZero \
    -p $THREADS \
    -o matrix/${MARK}_genes_norm.gz

plotProfile \
    -m matrix/${MARK}_genes_norm.gz \
    --yAxisLabel "RPGC-normalized signal" \
    --plotTitle "${MARK} gene-body (normalized)" \
    -out plots/${MARK}_genes_norm_profile.pdf
done
```

## CEBPA Peak-Focused Analysis

Now, to focus on peaks identified by MACS3 and validated by IDR, we'll identify promoters overlapping CEBPA consensus peaks.

### Identifying Promoters Overlapping CEBPA Peaks

```bash
# Count total CEBPA IDR-passed peaks
wc -l ceb_idr_passed.bed
# Output: 9468 ceb_idr_passed.bed

bedtools intersect \
  -a TSS.bed \
  -b ceb_idr_passed.bed \
  -u > cebpa_peak_promoters.bed


(chip) rajaishaqnabikhan@Mac 1.bigwig_smoothlength %  wc -l cebpa_peak_promoters.bed
    5792 cebpa_peak_promoters.bed
############################################
# 4. computeMatrix: TSS-centered
############################################

 Computing TSS-centered signal matrix 

computeMatrix reference-point \
  --referencePoint TSS \
  -b 2000 -a 2000 \
  -R cebpa_peak_promoters.bed \
  -S bw_log2/ceb_log2IPoverInput.bw \
  --binSize 25 \
  --skipZeros \
  -p 8\
  -o ceb_idr_TSS.mat.gz




############################################
# 5A. Heatmap WITHOUT clustering
############################################

 Plotting heatmap 

plotHeatmap \
  -m ceb_idr_TSS.mat.gz \
  --colorMap RdBu_r \
  --refPointLabel TSS \
  --dpi 600 \
  -out cebpa_peakPromoters_heatmap_noKmeans.pdf

echo



 ```

<img alt="Screenshot 2025-12-10 at 12 14 14 PM" src="./images/bw3.png" />

```


(chip) rajaishaqnabikhan@Mac 1.bigwig_smoothlength % 
computeMatrix reference-point \
  --referencePoint TSS \
  -b 2000 -a 2000 \
  -R cebpa_peak_promoters.bed \
  -S bw_log2/ceb_log2IPoverInput.bw \
     bw_log2/H3K9ac_log2IPoverInput.bw \
     bw_log2/H3K27me3_log2IPoverInput.bw \
  --binSize 25 \
  -p 8 \
  -o cebpa_promoters_TSS.mat.gz

(chip) rajaishaqnabikhan@Mac 1.bigwig_smoothlength %                     
plotProfile \                  
  -m cebpa_promoters_TSS.mat.gz \
  --perGroup \                      
  --plotTitle "CEBPA-bound promoters (±2 kb TSS)" \
  --dpi 600 \
-out cebpa_promoters_profile.pdf

(chip) rajaishaqnabikhan@Mac 1.bigwig_smoothlength % plotHeatmap \
  -m cebpa_promoters_TSS.mat.gz \
  --colorMap RdBu_r \
  --dpi 600 \
  -out cebpa_promoters_heatmap.pdf
(chip) rajaishaqnabikhan@Mac 1.bigwig_smoothlength % 
```

<img alt="Screenshot 2025-12-10 at 12 14 14 PM" src="./images/bw4.png" />
(chip) rajaishaqnabikhan@Mac 1.bigwig_smoothlength % plotHeatmap \
  -m cebpa_promoters_TSS.mat.gz \
  --colorMap RdBu_r \
  --dpi 600 \
  -out cebpa_promoters_heatmap.pdf
(chip) rajaishaqnabikhan@Mac 1.bigwig_smoothlength %
<img alt="Screenshot 2025-12-10 at 12 14 14 PM" src="./images/bw5.png" />```

CEBPA-bound promoters are enriched for active chromatin (H3K9ac) and depleted for repressive chromatin (H3K27me3).
