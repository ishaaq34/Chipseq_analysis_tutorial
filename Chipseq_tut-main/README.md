# Practical ChIP-Seq Analysis Tutorial

A comprehensive, beginner-friendly guide to performing ChIP-seq analysis, from raw data download to peak calling and visualization.

## Project Background
This tutorial series utilizes real-world data from the **ENCODE** consortium to investigate gene regulation mechanisms. Specifically, we analyze the genomic binding patterns of the transcription factor **CEBPA** and correlate them with two distinct histone modifications:
*   **H3K9ac (Histone 3 Lysine 9 Acetylation):** A marker of active gene promoters.
*   **H3K27me3 (Histone 3 Lysine 27 Trimethylation):** A marker of gene silencing (polycomb repression).

By processing these datasets in parallel, we learn to distinguish between sharp peaks (Transcription Factors), narrow peaks (Active Histones), and broad domains (Repressive Histones).

## Prerequisites
*   **Conda**: Ensure you have Anaconda or Miniconda installed.
*   **Environment**: Create the environment using the provided YAML file:
    ```bash
    conda env create -f 00_chip_updated.yml
    conda activate chip_seq
    ```

## Table of Contents

### Phase 1: Setup & Data Prep
*   **[00_setup_environment.md](00_setup_environment.md)**: Setting up your Conda environment.
*   **[01_geo_fastq_download.md](01_geo_fastq_download.md)**: How to download raw data (FASTQ) from GEO/SRA.
*   **[02_bash_automation.ipynb](02_bash_automation.ipynb)**: **(Essential)** Learning to automate tasks with Bash loops and arrays.
*   **[03_sample_list_creation.md](03_sample_list_creation.md)**: Creating a clean list of samples for processing.

### Phase 2: Alignment & QC
*   **[04_fastq_concepts.md](04_fastq_concepts.md)**: Understanding FASTQ file format and quality.
*   **[05_alignment_bowtie2.md](05_alignment_bowtie2.md)**: Aligning reads to the genome using Bowtie2.
*   **[06_duplicate_removal_qc.md](06_duplicate_removal_qc.md)**: Removing PCR duplicates with Picard/Samtools.
*   **[07_library_complexity.md](07_library_complexity.md)**: assessing library complexity.
*   **[08_bam_quality_metrics.md](08_bam_quality_metrics.md)**: Generating mapping statistics tables.

### Phase 3: Advanced QC & Peak Calling
*   **[09_strand_cross_correlation.md](09_strand_cross_correlation.md)**: Assessing signal-to-noise ratio (NSC/RSC).
*   **[10_bam_summary_fingerprint.md](10_bam_summary_fingerprint.md)**: Fingerprint plots and multiBAM summary.
*   **[11_macs2_peak_calling.md](11_macs2_peak_calling.md)**: **(Core)** Calling peaks (Narrow vs Broad) using MACS2.

### Phase 4: Visualization & Annotation
*   **[12_bigwig_generation.md](12_bigwig_generation.md)**: converting BAM to BigWig for IGV visualization.
*   **[13_visualization_heatmaps.md](13_visualization_heatmaps.md)**: Creating Heatmaps and Profile plots with deeptools.
*   **[14_chipseeker_annotation.Rmd](14_chipseeker_annotation.Rmd)**: Annotating peaks with genes using R/ChIPseeker.

---
*Created for the ChIP-seq Analysis Tutorial Series.*

## Dataset Used
The analysis in this tutorial follows a specific set of **ENCODE** data, analyzing **CEBPA** (Transcription Factor), **H3K27me3** (Broad Histone), and **H3K9ac** (Narrow Histone).

| Biosample Accession | ChIP Type | Target | Custom BAM Filename |
| :--- | :--- | :--- | :--- |
| **ENCFF327JFG** | TF ChIP-seq | CEBPA | `ceb_ENCFF327JFG.bam` |
| **ENCFF744SVA** | TF ChIP-seq | CEBPA | `ceb_ENCFF744SVA.bam` |
| **ENCFF164ALR** | Histone ChIP-seq | H3K27me3 | `H3K27me3_ENCFF164ALR.bam` |
| **ENCFF532DQH** | Histone ChIP-seq | H3K27me3 | `H3K27me3_ENCFF532DQH.bam` |
| **ENCFF193NPE** | Histone ChIP-seq | H3K9ac | `H3K9ac_ENCFF193NPE.bam` |
| **ENCFF534IPX** | Histone ChIP-seq | H3K9ac | `H3K9ac_ENCFF534IPX.bam` |
| **ENCFF110SOB** | Control ChIP-seq | Input | `Input_ENCFF110SOB.bam` |
| **ENCFF919XCV** | Control ChIP-seq | Input | `Input_ENCFF919XCV.bam` |

## Diagrams
A set of PDF diagrams illustrating the workflow concepts (e.g., Strand Cross-Correlation, Peak Calling logic) can be found in the `images/` directory.
