# Experiment Design & BAM Quality Control

`ENCODE` `BAM-files` `samtools` `flagstat` `MAPQ` `quality-control` `alignment-QC` `deepTools` `single-end`

## 1. Basic Concept (The Experiment & The File)

### 1.1 The Experiment Story

This tutorial uses real data from **BLaER1 cells** (human immune cells).

* **Treatment:** Cells were treated for 18 hours with Estradiol, Interleukin-3, and CSF1 to activate specific genes.
* **The Targets:**
    1. **CEBPA:** A Transcription Factor (The "Driver" that turns on genes).
    2. **H3K27me3:** A Histone Mark for **Closed/Repressed** DNA (The "Stop Sign").
    3. **H3K9ac:** A Histone Mark for **Open/Active** DNA (The "Go Sign").
    4. **Input:** Random background DNA (The "Noise" control).

### 1.2 The "Zip File" Analogy (BAM vs SAM)

* **SAM File (Sequence Alignment Map):** This is a huge, readable text file. It's like a 1000-page printed manuscript.
* **BAM File (Binary Alignment Map):** This is the **Compressed Zip File** version. It contains the exact same info but is smaller and faster for the computer to read.
  * *Rule:* We always work with BAM files to save space and time.

---

**Data Availability:**

* [ENCODE Cart](https://www.encodeproject.org/carts/ca521f95-7835-4369-88a9-b89f98fb39ad/)

## 2. Data used in the tutorial

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

---

### 2.2 Simplest and easiet way to donwload all the bam files in the working folder

```bash
#!/bin/bash
set -euo pipefail

cat <<EOF | xargs -n 2 -P 4 wget -c -O
ceb_ENCFF327JFG.bam        https://www.encodeproject.org/files/ENCFF327JFG/@@download/ENCFF327JFG.bam
ceb_ENCFF744SVA.bam        https://www.encodeproject.org/files/ENCFF744SVA/@@download/ENCFF744SVA.bam
H3K27me3_ENCFF164ALR.bam   https://www.encodeproject.org/files/ENCFF164ALR/@@download/ENCFF164ALR.bam
H3K27me3_ENCFF532DQH.bam   https://www.encodeproject.org/files/ENCFF532DQH/@@download/ENCFF532DQH.bam
H3K9ac_ENCFF193NPE.bam     https://www.encodeproject.org/files/ENCFF193NPE/@@download/ENCFF193NPE.bam
H3K9ac_ENCFF534IPX.bam     https://www.encodeproject.org/files/ENCFF534IPX/@@download/ENCFF534IPX.bam
Input_ENCFF110SOB.bam      https://www.encodeproject.org/files/ENCFF110SOB/@@download/ENCFF110SOB.bam
Input_ENCFF919XCV.bam      https://www.encodeproject.org/files/ENCFF919XCV/@@download/ENCFF919XCV.bam
EOF
```

## 3. Basic Quality Checks

Before processing, we verify the BAM files are healthy. For detailed explanations of these metrics, see [Section 05: Alignment QC](./05_alignment_bowtie2.md).

**3.1. Create a smaller test file (Optional)**
Working with full genomes takes time. For testing, we can extract just chromosome 11 and 12:

```bash
samtools view -b -h ceb_ENCFF327JFG.bam chr11 chr12 | samtools sort -o ceb_ENCFF327JFG.chr11_12.bam
```

**3.2. Get Alignment Stats**

For detailed explanation of `flagstat` output, see [Understanding BAM File Structure](./05_alignment_bowtie2.md#understanding-bam-file-structure).

```bash
# Quick summary
samtools flagstat ceb_ENCFF327JFG.bam > ceb_ENCFF327JFG.flagstat.txt
cat ceb_ENCFF327JFG.flagstat.txt
```

**Output:**

```text
2565563 + 0 in total (QC-passed reads + QC-failed reads)
2565563 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
2565563 + 0 mapped (100.00% : N/A)
2565563 + 0 primary mapped (100.00% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

**Verdict: 100% mapping rate - High-quality alignment, single-end data, no duplicates**

**3.3. Check MAPQ Distribution**

For detailed explanation of MAPQ scores and multimapping, see [Multimapping & The "Lost GPS"](./05_alignment_bowtie2.md#multimapping--the-lost-gps).

```bash
samtools view ceb_ENCFF327JFG.bam | awk '{print $5}' | sort -n | uniq -c
```

**Output:**

```text
123964  30  ← Lowest score is 30 (Good confidence)
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
2960636 42  ← Highest score is 42 (Perfect confidence)
```

**Verdict: MAPQ ≥ 30 - File contains only uniquely mapped, high-quality reads**

---

## Directory Structure After BAM QC

```text
chipseq_tutorial/
├── encode_bam/                  ← ENCODE BAM files
│   ├── ceb_ENCFF327JFG.bam
│   ├── H3K9ac_ENCFF193NPE.bam
│   └── ...                      (8 BAM files total)
├── encode_bam_qc/               ← QC metrics for ENCODE BAMs
│   ├── ceb_ENCFF327JFG.flagstat.txt
│   ├── ceb_ENCFF327JFG.mapq.txt
│   ├── H3K9ac_ENCFF193NPE.flagstat.txt
│   ├── H3K9ac_ENCFF193NPE.mapq.txt
│   └── ...                      (16 QC files total)
└── sample_id.txt
```

---

## Summary

1. **Context:** We are analyzing active/repressed marks in BLaER1 cells.
2. **Files:** BAMs are compressed alignment maps.
3. **QC:** We use `samtools flagstat` and check **MAPQ scores** to ensure we aren't analyzing "lost" multimapping reads.

> [!NOTE]
> **Up Next:** We'll perform cross-correlation analysis to estimate fragment length and assess signal quality.
