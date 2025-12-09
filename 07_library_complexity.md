# Tutorial 05: Library Complexity (The "Street Photographer")

## Level 1: Basic Concept (The Photographer)

Imagine you are a **Street Photographer** in a busy city. Your goal is to capture the diversity of the population.
*   **High Complexity Library (Good):** You take 100 photos, and every photo shows a different person. You have captured the true variety of the city.
*   **Low Complexity Library (Bad):** You take 100 photos, but it's just the same person 100 times. You wasted your film (sequencing reads) on duplicates.

In ChIP-seq, if we keep sequencing the exact same DNA fragment over and over (PCR duplicates), we aren't learning anything new. We want a "Complex" library with many unique fragments.

---

## Level 2: Execution (The Calculator)

We calculate complexity using metrics called **NRF** and **PBC**.
This calculation is a bit complex, so we use a script that combines `bedtools`, `sort`, and `awk`.

**Run this command block for one sample:**

```bash
# 1. Create a "5-prime" BED file
# (Convert BAM to simple coordinates, keeping only the start position of each read)
bedtools bamtobed -i bowalign_with_RG/Sample1.sorted.bam \
  | awk 'BEGIN{OFS="\t"} ($6=="+"){print $1,$2,$2+1} ($6=="-"){print $1,$3-1,$3}' \
  | sort -k1,1 -k2,2n \
  > QC_results/Sample1.read5.bed

# 2. Compute NRF, PBC1, PBC2
# (Count how many times each position appears)
uniq -c QC_results/Sample1.read5.bed \
  | awk '{c=$1; total+=c; uniq++; if(c==1) single++; if(c==2) double++;} \
    END{ if(total==0){print "NRF=NA\tPBC1=NA\tPBC2=NA"; exit} \
    NRF=uniq/total; \
    PBC1=single/uniq; \
    PBC2=(double? single/double:"Inf"); \
    printf("NRF=%.3f\tPBC1=%.3f\tPBC2=%s\n", NRF, PBC1, PBC2); }' \
  > QC_results/Sample1.pbc.txt
```

**Check the result:**
```bash
cat QC_results/Sample1.pbc.txt
```

---

## Level 3: Advanced Analysis (The Standards)

### Understanding the Metrics
*   **NRF (Non-Redundant Fraction):**  
    `Unique Reads / Total Reads`. (Ideal: > 0.8)
*   **PBC1 (PCR Bottleneck Coefficient 1):**  
    `Genomic positions with 1 read / Genomic positions with ≥1 read`. (Ideal: > 0.8)
*   **PBC2 (PCR Bottleneck Coefficient 2):**  
    `Positions with 1 read / Positions with 2 reads`. (Ideal: > 3.0)

### ENCODE Guidelines
How good is your library? Use this chart from ENCODE to grade your data.

<img width="1126" height="275" alt="Screenshot 2025-11-26 at 11 20 56 AM" src="https://github.com/user-attachments/assets/730d3662-6973-4a4e-a9f9-fe565583ed32" />

[*Source: ENCODE Data Standards*](https://www.encodeproject.org/data-standards/terms/#library)

### Before vs. After Deduplication
**Crucial Concept:**  
Raw data often looks "Low Complexity" just because of PCR duplicates. This is misleading. Once you remove the duplicates, the remaining data reveals the *true* quality of your library.

**Example Comparison:**

| Stage | NRF | PBC1 | PBC2 | Interpretation |
| :--- | :--- | :--- | :--- | :--- |
| **Before** Removal | 0.668 | 0.695 | 3.3 | Appears "Moderate/Low" quality due to duplicates. |
| **After** Removal | **0.952** | **0.949** | **18.6** | Use these values! True library is **High Quality**. |

**Lesson:** Don't panic if your raw NRF is low. Remove duplicates first, then check again.

---

## Summary
1.  **Goal:** Ensure we have many unique DNA fragments (High Complexity).
2.  **Action:** Run the NRF/PBC calculation script.
3.  **Result:** Compare your numbers against the ENCODE chart to validate your experiment.
