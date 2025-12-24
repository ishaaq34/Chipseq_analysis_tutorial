# Quality Metrics - FRiP (Fraction of Reads in Peaks)

`FRiP` `quality-control` `ChIP-seq` `ENCODE-standards` `peak-calling` `QC-metrics`

## Background

Having called peaks with MACS3, we now quantify ChIP enrichment quality using **FRiP** (Fraction of Reads in Peaks). FRiP measures the proportion of sequencing reads that overlap called peak regions versus background.

In successful ChIP-seq experiments, immunoprecipitation selectively enriches DNA fragments bound by the target protein, causing a substantial fraction of reads to concentrate in peak regions. The remainder represents genomic background from non-specific binding or incomplete washing. By contrast, failed experiments show most reads distributed evenly across the genome, yielding low FRiP values.

FRiP thus provides a rapid, quantitative metric for immunoprecipitation efficiency—high FRiP indicates successful enrichment, while low FRiP flags potential technical problems before investing effort in downstream analysis.

---

## 1. What is FRiP?

FRiP measures how much of your ChIP signal is concentrated in called peaks. It's a key ENCODE quality metric.

### Calculating FRiP

**Formula:** FRiP = (Reads in peaks) / (Total mapped reads)

**Input files:**

* BAM file: `encode_bam/ceb_ENCFF327JFG.bam` (deduplicated, MAPQ filtered)
* Peak file: `macs3_results/ceb_ENCFF327JFG_peaks.narrowPeak`

**Complete FRiP calculation script:**

```bash
# Step 1: Count total mapped reads
TOTAL_READS=$(samtools view -c encode_bam/ceb_ENCFF327JFG.bam)
```

This command counts the total number of alignments in the BAM file. Since the BAM is already deduplicated and MAPQ-filtered, this value represents the total number of mapped reads, which serves as the denominator in the FRiP calculation.

```bash
# Step 2: Merge overlapping peak regions (avoid double-counting)
mkdir -p frip_analysis
sort -k1,1 -k2,2n macs3_results/ceb_ENCFF327JFG_peaks.narrowPeak | \
  bedtools merge -i stdin > frip_analysis/ceb_ENCFF327JFG.peaks.merged.bed
```

 MACS3 peak files may contain overlapping or adjacent peak regions. To avoid double-counting reads that overlap multiple nearby peaks, all peak intervals are first sorted and then merged into a non-redundant set of regions. This step ensures that each read contributes at most once to the FRiP numerator.

```bash
# Step 3: Count reads overlapping peaks (unique reads only)
READS_IN_PEAKS=$(
  samtools view -b encode_bam/ceb_ENCFF327JFG.bam | \
  bedtools intersect -u -a stdin -b frip_analysis/ceb_ENCFF327JFG.peaks.merged.bed | \
  samtools view -c
)
```

Aligned reads from the BAM file are intersected with the merged peak regions. The `-u` option ensures that each read is counted only once, even if it overlaps more than one peak interval. The resulting count corresponds to the number of unique reads that fall within CEBPA peak regions.

```bash
# Step 4: Calculate the FRiP score
FRIP=$(awk -v n="${READS_IN_PEAKS}" -v d="${TOTAL_READS}" \
  'BEGIN {printf "%.5f", n/d}')
```

The FRiP score is calculated as the ratio of reads overlapping peak regions to the total number of mapped reads. The result is formatted to five decimal places for clarity and consistency.

### Complete FRiP Calculation Script

For convenience, here's the complete executable script combining all steps:

```bash
#!/bin/bash
# Complete FRiP calculation for ceb_ENCFF327JFG
mkdir -p frip_analysis

# Count total reads
TOTAL_READS=$(samtools view -c encode_bam/ceb_ENCFF327JFG.bam)

# Merge overlapping peaks
sort -k1,1 -k2,2n macs3_results/ceb_ENCFF327JFG_peaks.narrowPeak | \
  bedtools merge -i stdin > frip_analysis/ceb_ENCFF327JFG.peaks.merged.bed

# Count reads in peaks
READS_IN_PEAKS=$(samtools view -b encode_bam/ceb_ENCFF327JFG.bam | \
  bedtools intersect -u -a stdin -b frip_analysis/ceb_ENCFF327JFG.peaks.merged.bed | \
  samtools view -c)

# Calculate FRiP
FRIP=$(awk -v n="${READS_IN_PEAKS}" -v d="${TOTAL_READS}" 'BEGIN {printf "%.5f", n/d}')

# Report results
echo "Sample: ceb_ENCFF327JFG"
echo "Total mapped reads : ${TOTAL_READS}"
echo "Reads in peaks     : ${READS_IN_PEAKS}"
echo "FRiP               : ${FRIP}"
```

---

### Expected Output

**Results we got:**

```text
sample total_reads reads_in_peaks FRiP
ceb_ENCFF327JFG 27863042 1511077 0.05423
ceb_ENCFF744SVA 16814770 1187486 0.07062
H3K27me3_ENCFF164ALR 39952273 11528033 0.28855
H3K27me3_ENCFF532DQH 35184227 10987930 0.31230
H3K9ac_ENCFF193NPE 34567738 14283991 0.41322
H3K9ac_ENCFF534IPX 39792428 6625200 0.16649

```

### Interpreting These Results

**Key Observations:**

* **H3K9ac:** Replicate discordance (41% vs 17%) needs IDR validation
* **H3K27me3:** Excellent concordance (~30%) indicates robust data  
* **CEBPA:** Low (5-7%) but acceptable for sparse transcription factors

---

### FRiP Quality Standards

ENCODE guidelines ([Landt et al. 2012](https://doi.org/10.1101/gr.136184.111)) flag transcription factor experiments with FRiP below ~1%, as such values often reflect weak enrichment or technical failure, while well-defined point-source factors like CTCF or REST frequently achieve FRiP values of 20–50%. These expectations are not universal: factors with few binding sites or low occupancy can yield FRiP <1% yet remain biologically valid.

Histone marks generally show higher FRiP due to their broader genomic coverage, though acceptable values differ markedly between narrow and broad marks.

Because FRiP depends strongly on sequencing depth, peak-calling parameters, and the biological target, it should always be interpreted in conjunction with complementary QC metrics such as fingerprint plots, strand cross-correlation, and coverage profiles, not in isolation.

---

## Directory Structure After FRiP Calculation

```text
chipseq_tutorial/
├── encode_bam/                  ← Input BAM files
│   ├── ceb_ENCFF327JFG.bam
│   └── ...
├── macs3_results/               ← Peak files from MACS3
│   ├── ceb_ENCFF327JFG_peaks.narrowPeak
│   ├── H3K9ac_ENCFF534IPX_peaks.narrowPeak
│   └── ...
└── frip_analysis/               ← FRiP working files (optional)
    ├── ceb_ENCFF327JFG.peaks.merged.bed
    └── frip_scores.txt
```

---

## Summary

You've learned how to:

* Calculate FRiP (Fraction of Reads in Peaks) quality metrics
* Interpret FRiP values according to ENCODE guidelines
* Understand what constitutes good vs. poor ChIP-seq data

> [!NOTE]
> **Up Next:** Assess reproducibility between replicates using IDR (Irreproducible Discovery Rate) and identify consensus peaks for downstream motif analysis.

---
