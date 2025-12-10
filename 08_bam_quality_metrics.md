# Tutorial 06: Experiment Design & BAM Quality Control

## Level 1: Basic Concept (The Experiment & The File)

### 1.1 The Experiment Story
This tutorial uses real data from **BLaER1 cells** (human immune cells).
*   **Treatment:** Cells were treated for 18 hours with Estradiol, Interleukin-3, and CSF1 to activate specific genes.
*   **The Targets:**
    1.  **CEBPA:** A Transcription Factor (The "Driver" that turns on genes).
    2.  **H3K27me3:** A Histone Mark for **Closed/Repressed** DNA (The "Stop Sign").
    3.  **H3K9ac:** A Histone Mark for **Open/Active** DNA (The "Go Sign").
    4.  **Input:** Random background DNA (The "Noise" control).

### 1.2 The "Zip File" Analogy (BAM vs SAM)
*   **SAM File (Sequence Alignment Map):** This is a huge, readable text file. It's like a 1000-page printed manuscript.
*   **BAM File (Binary Alignment Map):** This is the **Compressed Zip File** version. It contains the exact same info but is smaller and faster for the computer to read.
    *   *Rule:* We always work with BAM files to save space and time.

---

## Level 2: Execution (The Checklist)

### 2.1 Sample Table
Here are the files we are analyzing. In real life, you should make a table like this to track your work.

| Biosample Accession | ChIP Type        | Target     | Custom BAM Filename          |
|---------------------|------------------|------------|-------------------------------|
| ENCFF327JFG       | TF ChIP-seq      | CEBPA      | ceb_ENCFF327JFG.bam         |
| ENCFF744SVA        | TF ChIP-seq      | CEBPA      | ceb_ENCFF744SVA.bam          |
| ENCFF164ALR        | Histone ChIP-seq | H3K27me3   | H3K27me3_ENCFF164ALR.bam     |
| ENCFF532DQH      | Histone ChIP-seq | H3K27me3   | H3K27me3_ENCFF532DQH.bam      |
| ENCFF193NPE        | Histone ChIP-seq | H3K9ac     | H3K9ac_ENCFF193NPE.bam       |
| ENCFF534IPX      | Histone ChIP-seq | H3K9ac     | H3K9ac_ENCFF534IPX.bam       |
| ENCFF110SOB        | Control ChIP-seq | Input      | Input_ENCFF110SOB.bam         |
| ENCFF919XCV      | Control ChIP-seq | Input      | Input_ENCFF919XCV.bam        |

### 2.2 Basic Checks
Before processing, we verify the BAM files are healthy.

**1. Create a smaller test file (Optional)**
Working with full genomes takes time. For testing, we can extract just chromosome 11 and 12:
```bash
samtools view -b -h sample.bam chr11 chr12 | samtools sort -o sample.chr11_12.bam
```

**2. Get Alignment Stats**
```bash
# Quick summary
samtools flagstat sample.bam > sample.flagstat.txt

# Detailed stats
samtools stats sample.bam > sample.stats.txt
```

---

## Quality Control

### 3.1 Interpreting `flagstat`
What do the numbers mean?

```text
2565563 + 0 in total          # Total raw reads.
2565563 + 0 mapped (100%)     # Success! 100% of reads found a home on the genome.
0 + 0 duplicates              # 0 duplicates found (if un-marked).
0 + 0 paired in sequencing    # These are Single-End reads (0 paired).
```
*   **Goal:** High mapping % (>80%).
*   **Warning:** If mapping is low (<50%), you may have the wrong organism or bad sequencing.

### 3.2 Multimapping & The "Lost GPS"
Sometimes a read is repetitive (e.g., "ATATATAT"). It fits in 50 different places on the genome.
The aligner doesn't know which spot is correct, so it gives it a **Low MAPQ Score** (Mapping Quality).

*   **MAPQ = 0:** "I have NO clue where this goes. It fits in many places." (Multimapped)
*   **MAPQ > 30:** "I am highly confident this read belongs EXACTLY here." (Unique)

**How to check MAPQ distribution:**
```bash
samtools view sample.bam | awk '{print $5}' | sort -n | uniq -c
```

**Understanding the Output:**

*   **If you see lots of 0s:** Your file includes "Lost GPS" reads (Multimappers).
*   **If your scores start at 30+:** Your file has **Already Filtered** the bad reads.

**Example of a Clean File (Filtered):**
```text
123964  30  <-- Lowest score is 30 (Good confidence)
1928 31
20741 32
1477 33
4040 34
34434 35
14294 36
53329 37
30665 38
88528 39
77123 40
2960636 42  <-- Highest score is 42 (Perfect confidence)
```
*Verdict: This BAM file contains only uniquely mapped, high-quality reads.*

---

## Summary
1.  **Context:** We are analyzing active/repressed marks in BLaER1 cells.
2.  **Files:** BAMs are compressed alignment maps.
3.  **QC:** We use `samtools flagstat` and check **MAPQ scores** to ensure we aren't analyzing "lost" multimapping reads.
