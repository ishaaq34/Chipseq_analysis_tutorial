# Library Complexity (The "Street Photographer")

`library-complexity` `NRF` `PBC` `PCR-duplicates` `bedtools` `unique-reads` `QC-metrics`

## Level 1: Basic Concept (The Photographer)

Imagine you are a **Street Photographer** in a busy city. Your goal is to capture the diversity of the population.

* **High Complexity Library (Good):** You take 100 photos, and every photo shows a different person. You have captured the true variety of the city.
* **Low Complexity Library (Bad):** You take 100 photos, but it's just the same person 100 times. You wasted your film (sequencing reads) on duplicates.

In ChIP-seq, if we keep sequencing the exact same DNA fragment over and over (PCR duplicates), we aren't learning anything new. We want a "Complex" library with many unique fragments.

---

## Level 2: Execution (The Calculator)

We calculate complexity using metrics called **NRF** and **PBC**.
This calculation is a bit complex, so we use a script that combines `bedtools`, `sort`, and `awk`.

**Input files needed:**

```text
chipseq_tutorial/
├── picard_dedup_bam/
│   ├── H3K27me3_IP_rep1.dedup.bam
│   └── ...
├── QC_results/
│   ├── H3K27me3_IP_rep1.pbc.txt
│   └── ...
└── sample_id.txt
```

**Run this command block for one sample:**

```bash
mkdir -p QC_results

# 1. Create a "5-prime" BED file
# (Convert BAM to simple coordinates, keeping only the start position of each read)
bedtools bamtobed -i picard_dedup_bam/H3K27me3_IP_rep1.dedup.bam \
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

```text
NRF=0.997 PBC1=0.997 PBC2=360.061
```

---

### Understanding the Metrics

* **NRF (Non-Redundant Fraction):**  
    `Unique Reads / Total Reads`. (Ideal: > 0.8)
* **PBC1 (PCR Bottleneck Coefficient 1):**  
    `Genomic positions with 1 read / Genomic positions with ≥1 read`. (Ideal: > 0.8)
* **PBC2 (PCR Bottleneck Coefficient 2):**  
    `Positions with 1 read / Positions with 2 reads`. (Ideal: > 3.0)

### ENCODE Guidelines

How good is your library? Use this chart from ENCODE to grade your data.

<img width="1126" height="275" alt="Screenshot 2025-11-26 at 11 20 56 AM" src="./images/library_complexity_plot.png" />

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

1. **Goal:** Ensure we have many unique DNA fragments (High Complexity).
2. **Action:** Run the NRF/PBC calculation script.
3. **Result:** Compare your numbers against the ENCODE chart to validate your experiment.

> [!NOTE]
> **Up Next:** We'll dive deeper into BAM quality metrics to assess alignment quality and fragment sizes.

---

## Pipeline Summary: Pre-processing Workflow Complete

Congratulations! You have completed the core ChIP-seq pre-processing workflow from raw sequencing reads to analysis-ready BAM files.

**Pre-processing Pipeline:**

* **[FASTQ Acquisition](./02_geo_fastq_download.md):** Downloaded from GEO/SRA (C. elegans H3K27me3 ChIP-seq, single-end reads)
* **[Sample Manifest](./03_sample_list_creation.md):** Created sample ID list for batch processing
* **[Quality Control](./04_fastq_concepts.md):** Read quality assessment and adapter trimming (fastp)
* **[Genome Alignment](./05_alignment_bowtie2.md):** Mapped reads to reference assembly (Bowtie2, single-end mode)
* **[Duplicate Removal](./06_duplicate_removal_qc.md):** PCR duplicate identification and removal (Picard MarkDuplicates)
* **[Library Complexity](./07_library_complexity.md):** QC assessment with NRF, PBC1, PBC2 metrics

**Analysis-ready outputs:**

* Deduplicated, coordinate-sorted BAM files with index (.bai)
* Per-sample QC metrics (alignment rate, duplication rate, complexity)
* Quality-filtered reads ready for peak calling

---

## Transition to ENCODE BAM Files

> [!IMPORTANT]
> **For downstream analysis**, we will use **pre-processed BAM files from ENCODE** instead of the files we generated during pre-processing.
>
> **Why?**
>
> * ENCODE provides high-quality, standardized ChIP-seq data
> * Files are already aligned, deduplicated, and quality-controlled
> * This allows us to focus on downstream analysis (peak calling, visualization, annotation)
>
> **What this means:**
>
> * The pre-processing workflow taught you complete data preparation from raw reads to analysis-ready BAM
> * Downstream analysis will use ENCODE BAM files for demonstration
> * You can apply both approaches: process your own data OR use public ENCODE data

The skills you learned in the pre-processing workflow are valuable for processing your own ChIP-seq data. For the remaining tutorials, we'll demonstrate downstream analysis using ENCODE's curated datasets.

> [!TIP]
> **Challenge yourself!** Before moving to downstream analysis, consider applying the downstream analysis steps (peak calling, visualization, annotation) to your C. elegans H3K27me3 data from the pre-processing workflow. This hands-on practice will solidify your understanding and give you complete end-to-end ChIP-seq analysis experience. You can then compare your results with the ENCODE workflow!
