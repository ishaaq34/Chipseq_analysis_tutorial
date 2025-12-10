# Tutorial 07: Advanced QC with deepTools (The "Health Checkup")

## Basic Concept (The Health Check)

Before we call peaks, we must perform a **Health Checkup** on our data.
*   **The Census (Fingerprint):** Are the reads spread out evenly (Socialist/Input) or concentrated in specific spots (Capitalist/ChIP)?
*   **The Coverage Check:** Do we have enough reads? Are they duplicates?
*   **The Family Tree (Correlation & PCA):** Do biological replicates (siblings) look similar? Are they distinct from the control?

We use a suite of tools called **deepTools** to generate these reports.

---

## Execution (Running the Tests)

### 2.1 Fingerprint Plot
Checks if the IP worked (Enrichment).
```bash
plotFingerprint \
  -b sample1.bam sample2.bam sample3.bam \              # BAM files (ChIP + input) to evaluate
  --skipZeros \                       # Ignore bins with zero coverage (speeds up, reduces noise)
  --numberOfSamples 50000 \           # Number of random genomic positions for sampling
  -T "Fingerprints of all BAM samples" \  # Plot title
  --plotFile fingerprints.pdf \  # Output PDF
  --plotFileFormat pdf \
  2>&1 | tee plotFingerprint.log  # Log stderr+stdout to a log file
```

### 2.2 Coverage Plot
Checks sequencing depth and duplication levels.
```bash
plotCoverage \
  -b  sample1.bam sample2.bam sample3.bam \                           # BAM files to summarize
  -o  coverage_histogram.pdf \                             # Output PDF with coverage histogram(s)
  --plotFileFormat pdf \
  --smartLabels \                             # Auto-generate labels from filenames (if labels not passed)
  --numberOfSamples 1000000 \                 # Number of genomic positions to sample
  --ignoreDuplicates \                        # Exclude duplicate reads (avoid PCR bias)
  --minMappingQuality 30 \                    # Only count high-quality alignments (MAPQ ≥ 30)
  --outRawCounts coverage_counts.txt          # Raw coverage counts for downstream inspection
```

### 2.3 Summary Matrix
We need a "Count Matrix" to compare samples. This counts reads in bins across the whole genome.
```bash
multiBamSummary bins \
  -b  sample1.bam sample2.bam sample3.bam \                          
  --numberOfProcessors 4 \        # Parallel processing
  -o matrix.npz \                              # Output compressed matrix (used by other deepTools commands)
  --outRawCounts matrix.tab                   # Tab-delimited matrix of counts per bin per sample
```

### 2.4 Correlation & PCA
Using the matrix from Step 2.3, we compare the samples.

**Correlation Heatmap:**
```bash
plotCorrelation \
  -in matrix.npz \                                   # Input matrix from multiBamSummary
  --corMethod spearman \                           # Spearman correlation (rank-based, robust to outliers)
  --skipZeros \                                    # Ignore bins with zero signal in all samples
  --whatToPlot heatmap \                           # Produce a heatmap
  --plotNumbers \                                  # Print correlation coefficients on the heatmap
  -o  spearman_corr_plot.pdf \                            # Output PDF
  --outFileCorMatrix spearman_corr_plot.tab                 # Save the correlation matrix to a text file
```

**PCA (Principal Component Analysis):**
```bash
plotPCA \
  -in matrix.npz \                                   # Same matrix as for plotCorrelation
  -o pca.pdf \                                # Output PCA plot
  --transpose \                                    # Treat samples as variables (standard deepTools setting here)
  --plotWidth 10 \                                 # Width of figure (in inches)
  --plotHeight 8 \                                 # Height of figure (in inches)
  --plotFileFormat pdf \                           # Output format
  --outFileNameData pca.tab                 # Save underlying PCA coordinates
```

---

## Level 3: Advanced Analysis (Reading the Charts)

### 3.1 Interpreting the Fingerprint (The Census)
This plot shows the cumulative read distribution across the genome. Good ChIP libraries show a clear separation between ChIP and input samples, with ChIP curves rising earlier due to enriched regions. Flat, overlapping curves usually indicate poor enrichment or over-background signal.

---
<img width="755" height="574" alt="Screenshot 2025-12-10 at 12 02 40 PM" src="https://github.com/user-attachments/assets/db1e7f33-6775-4705-b6b3-62a4c0bf7405" />

---


| Sample Type | Interpretation |
| :--- | :--- |
| **Input** | Close to the diagonal. Reads are uniformly distributed, behaving like ideal background. |
| **H3K9ac (Active)** | Strong "Elbow". Reads concentrated in a small fraction of bins (focal peaks), showing a strong bend away from the diagonal. |
| **H3K27me3 (Broad)** | In between. Enrichment is present but more diffuse, so the curve is less sharp than H3K9ac. |
| **ceb ChIP** | Slightly lower than Input. Modest enrichment, close to background. |

### 3.2 Interpreting Coverage
Next, we look at the overall coverage distribution in each BAM using plotCoverage. This reveals whether some samples are globally under-sequenced, dominated by a few high-coverage regions, or heavily affected by duplicated reads. We restrict to high-quality, non-duplicate reads to make the distributions comparable.

**Plot A: The Drop-off**
*   **Inputs :** The Input tracks sit higher at low coverage because they spread their reads across the genome without enrichment. That’s why both Input samples show a large fraction of bases at coverage 0 and 1, then taper off more slowly as coverage increases.
  
*   **ChIPs :** In contrast, every IP sample collapses more sharply; the curves drop faster after coverage 1 because most genomic positions in a ChIP experiment receive almost no reads. Only a small portion of the genome — the actual binding or modification sites — reaches deeper coverage, and that fraction is tiny enough that the tail beyond coverage 2 nearly vanishes.

---
<img width="861" height="578" alt="Screenshot 2025-12-10 at 12 03 51 PM" src="https://github.com/user-attachments/assets/26404eb5-7041-4674-bf94-d44d0d9edc8b" />

---


**Plot B: The Tail**
Zooming in reveals the difference. Input covers more of the genome at 1x depth, while ChIP focuses on peaks. The Input curves decline more slowly because a larger fraction of their genome maintains at least some measurable coverage. The IP curves fall off earlier and more steeply, which reflects the enrichment pattern: most positions have essentially no reads, and only a very small subset of bases in true peak regions sustain higher coverage.

---

<img width="731" height="501" alt="Screenshot 2025-12-10 at 12 04 14 PM" src="https://github.com/user-attachments/assets/b8e24d14-d41b-4c55-a268-9982b49026c5" />

---

### 3.3 Interpreting Correlation (The Family Tree)
Using the binned count matrix, we compute pairwise correlations between samples. A Spearman correlation heatmap shows whether biological replicates cluster together and whether inputs are distinct from ChIP samples. Poor clustering or scattered correlations usually indicate sample swaps, failed IPs, or inconsistent library prep.

**The Heatmap:**
*   **Clustering:** H3K27me3 samples cluster together and show moderately high mutual correlations (around 0.6), which is exactly what you expect for a broad repressive mark. The H3K9ac samples also correlate strongly with each other (0.76–1.0), forming a clean sub-cluster that is distinct from H3K27me3.
*   **Separation:** Active marks (H3K9ac) should look different from Repressive marks (H3K27me3).

---
<img width="667" height="575" alt="Screenshot 2025-12-10 at 12 05 28 PM" src="https://github.com/user-attachments/assets/652741bb-8b63-4fc2-88d8-eca5508e5938" />

---



### 3.4 Interpreting PCA
Finally, we perform PCA on the same binned count matrix. PCA reduces the data to a few dimensions that capture most of the variance. In a good ChIP-seq dataset, biological replicates cluster together in PCA space, and distinct conditions or marks separate along major components. PCA is a convenient visual check for batch effects, sample swaps, and outlier libraries.

**The Map:**
*   **PC1 (X-axis):** The PCA shows clear separation of samples by assay type. PC1 captures most of the variance and cleanly splits the H3K27me3 group from the H3K9ac group, which is expected because these marks have very different genomic distributions.
*   **Clustering:** The two “ceb” samples cluster tightly together, indicating consistent coverage patterns within that group. The H3K27me3 replicates are also tightly paired, which matches their broad and uniform enrichment profile. The H3K9ac replicates sit on the opposite side of PC1, with one replicate shifted slightly on PC2, hinting at a mild difference in coverage distribution but nothing severe. If one replicate is far away, it might be an outlier/bad sample.

---

<img width="983" height="488" alt="Screenshot 2025-12-10 at 12 05 56 PM" src="https://github.com/user-attachments/assets/d2d785d1-4da3-4cc4-a3e0-1ee612ad3d38" />

---

## Summary
1.  **Fingerprint:** Confirms your IP worked (Sharp elbow).
2.  **Coverage:** Confirms sequencing depth (Inputs = broad, ChIP = peaky).
3.  **PCA/Correlation:** Confirms your replicates match (Siblings cluster together).
