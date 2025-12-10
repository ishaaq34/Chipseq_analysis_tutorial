# Tutorial 08: Peak Calling with MACS2 (The Summit Search)

## Basic Concept (The "Heap" Hunt)

### What is a Peak?
Imagine you are looking for hidden treasure on a long beach (the Genome).
*   **Random Reads:** Sand scattered everywhere evenly (Background noise).
*   **A "Peak":** A huge heap of sand in one specific spot.

This heap means thousands of protein molecules were bound to that exact spot of DNA.

### Signal vs. Noise
*   **Signal:** The Heap (Peak).
*   **Noise:** The random thin layer of sand everywhere else (Input/Control).

**MACS2** is the software that scans the beach, measures the height of the sand pile, and calculates if it is "significantly" higher than the background.

---

## Execution (Step-by-Step)

We will not use complex loops. We will run the command for each sample type specifically.

### Step 1: Create Output Directory
Keep your workspace clean.
```bash
mkdir -p macs2_results
```

### Step 2: H3K9ac (Narrow Peak)
H3K9ac (Acetylation) creates sharp, tall spikes. This is the default mode for MACS2.

**The Command:**
```bash
macs2 callpeak \
  -t H3K9ac.bam \
  -c Input.bam \
  -f BAM \
  -g hs \
  -n H3K9ac \
  -q 0.01 \
  --outdir macs2_results
```

**Explanation:**
*   `-t`: Treatment file (Signal).
*   `-c`: Control file (Input/Background).
*   `-f`: Format (BAM).
*   `-g`: Scalable Genome Size (`hs` for human, `mm` for mouse).
*   `-n`: Name of the output prefix.
*   `-q 0.01`: The cutoff. Only keep peaks with a q-value (FDR) better than 0.01.

### Step 3: H3K27me3 (Broad Peak)
H3K27me3 (Methylation) creates wide, gentle hills. MACS2 needs to be told to look for "Broad" regions.

**The Command:**
```bash
macs2 callpeak \
  -t H3K27me3.bam \
  -c Input.bam \
  -f BAM \
  -g hs \
  -n H3K27me3 \
  --broad \
  --broad-cutoff 0.1 \
  --outdir macs2_results
```

**Explanation:**
*   `--broad`: Tells MACS2 to link nearby peaks together into one big region.
*   `--broad-cutoff 0.1`: We accept a slightly weaker signal (0.1) because the signal is spread out over a wider area.

### Step 4: Transcription Factors (e.g., CEBPA)
TFs bind very tightly to tiny spots. We treat them exactly like H3K9ac (Narrow).

**The Command:**
```bash
macs2 callpeak \
  -t CEBPA.bam \
  -c Input.bam \
  -f BAM \
  -g hs \
  -n CEBPA \
  -q 0.01 \
  --outdir macs2_results
```

---

## The Outputs

After running the commands, look in the `macs2_results` folder. You will see these files:

1.  **`NAME_peaks.narrowPeak`**:
    *   **What is it?** The final list of peaks.
    *   **Format:** BED format (Chr, Start, End, Name, Score...).
    *   **Use:** Open this in IGV to see the bars under your peaks.

2.  **`NAME_peaks.xls`**:
    *   **What is it?** An Excel-friendly table.
    *   **Contents:** Coordinates, fold-enrichment, p-values, and q-values.

3.  **`NAME_summits.bed`**:
    *   **What is it?** The single highest point (1bp) of each peak.
    *   **Use:** Important for finding motifs (the exact DNA sequence the protein touched).
