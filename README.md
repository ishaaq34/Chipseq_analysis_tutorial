# ChIP-seq Analysis Tutorial

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![mdBook](https://img.shields.io/badge/mdBook-Documentation-blue)](https://ishaaq34.github.io/Chipseq_analysis_tutorial/)

**[ Read the Tutorial Book Online](https://ishaaq34.github.io/Chipseq_analysis_tutorial/)**


=======
A comprehensive, hands-on tutorial for ChIP-seq (Chromatin Immunoprecipitation Sequencing) data analysis, covering the complete workflow from raw FASTQ files to publication-ready visualizations and biological insights.

## Overview

This tutorial provides a practical, step-by-step guide to analyzing ChIP-seq data using established bioinformatics tools and ENCODE best practices. It's designed for researchers, graduate students, and bioinformaticians who want to learn or improve their ChIP-seq analysis skills.

**Key Features:**

- Complete end-to-end workflow (FASTQ → peaks → visualization → annotation)
- Real datasets from GEO and ENCODE
- Reproducibility-focused with IDR analysis
- Quality control at every step
- Extensive visualization techniques
- Bash automation examples

## Tutorial Structure

The tutorial consists of 16 comprehensive chapters plus appendices:

### Setup & Data Acquisition

**00. [Welcome & Introduction](src/00_introduction.md)**

- Overview of ChIP-seq technology and applications
- Tutorial objectives and workflow outline
- Prerequisites and learning outcomes

**01. [Environment Setup](src/01_setup_environment.md)**

- Installing required software (conda, mamba)
- Creating the analysis environment
- Tool verification and troubleshooting

**02. [Bash Automation Fundamentals](src/02_bash_automation.md)**

- Essential bash scripting for bioinformatics
- Loops for batch processing
- Creating sample ID lists from filenames
- Automation best practices

**03. [GEO/FASTQ Download](src/03_geo_fastq_download.md)**

- Navigating NCBI GEO database
- Downloading ChIP-seq datasets
- Understanding SRA and FASTQ formats
- Using sra-toolkit and wget

**04. [FASTQ Concepts & QC](src/04_fastq_concepts.md)**

- FASTQ format deep dive
- Quality control with fastp
- Adapter trimming and filtering
- Batch processing multiple samples

### Alignment & Initial QC

**05. [Alignment with Bowtie2](src/05_alignment_bowtie2.md)**

- Building genome indices
- Aligning reads to reference genome
- Understanding SAM/BAM formats
- Sorting and indexing BAM files

**06. [Duplicate Removal & QC](src/06_duplicate_removal_qc.md)**

- PCR duplicate detection
- Using Picard MarkDuplicates
- Filtering for quality and uniqueness
- Creating final analysis-ready BAM files

**07. [Library Complexity Assessment](src/07_library_complexity.md)**

- Estimating library complexity
- Preseq analysis for saturation
- Interpreting complexity curves
- Quality metrics interpretation

**08. [BAM Quality Metrics](src/08_bam_quality_metrics.md)**

- Computing alignment statistics
- MAPQ score analysis
- Insert size distributions
- Generating MultiQC reports

**09. [Strand Cross-Correlation](src/09_strand_cross_correlation.md)**

- Fragment length estimation
- NSC and RSC quality metrics
- phantompeakqualtools analysis
- Interpreting cross-correlation plots

**10. [BAM Summary & Fingerprint Plots](src/10_bam_summary_fingerprint.md)**

- deepTools multiBamSummary
- Fingerprint plot generation
- Sample correlation analysis
- PCA visualization

### Peak Calling & Reproducibility

**11. [MACS3 Peak Calling](src/11_macs3_peak_calling.md)**

- Narrow vs broad peak calling
- MACS3 parameters and options
- Understanding peak output formats
- Model building and interpretation

**12. [FRiP Quality Metrics](src/12_frip_quality_metrics.md)**

- Fraction of Reads in Peaks (FRiP)
- Calculating signal-to-noise ratios
- Quality thresholds (ENCODE standards)
- Batch FRiP calculation scripts

**13. [IDR & Consensus Peaks](src/13_idr_consensus_motifs_rk_corrected.md)**

- Irreproducible Discovery Rate (IDR) analysis
- Reproducibility assessment between replicates
- Creating consensus peak sets
- HOMER motif discovery
- Biological interpretation

### Visualization & Annotation

**14. [BigWig Generation](src/14_bigwig_generation.md)**

- BAM to BigWig conversion
- RPGC normalization
- Effective genome size concepts
- Creating genome browser tracks

**15. [Visualization with deepTools](src/15_visualization_heatmaps.md)**

- TSS-centered profile plots
- Gene body heatmaps
- BigWig averaging and normalization
- Log2(IP/Input) ratios
- Multi-mark chromatin state analysis
- Understanding pseudocount in normalization

**16. [Peak Annotation with ChIPseeker](src/16_chipseeker_annotation.md)**

- Genomic feature annotation
- Promoter/exon/intron distribution
- Functional enrichment analysis
- Biological interpretation of binding patterns

## Datasets

The tutorial uses two publicly available datasets:

### Dataset 1: C. elegans (Preprocessing Chapters 1-10)

- **GEO Accession:** [GSE115704](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115704)
- **Study:** Histone modifications in C. elegans sperm, oocytes, and early embryos
- **Organism:** Caenorhabditis elegans
- **DOI:** [10.1093/nar/gky1229](https://doi.org/10.1093/nar/gky1229)
- **Used for:** Teaching alignment, QC, and preprocessing workflows

### Dataset 2: Human K562 ENCODE (Peak Calling & Analysis Chapters 11-16)

- **Source:** [ENCODE BLaER1 ChIP-seq](https://www.encodeproject.org/carts/ca521f95-7835-4369-88a9-b89f98fb39ad/)
- **Cell Line:** K562 (human erythroleukemia)
- **Marks:**
  - CEBPA (transcription factor)
  - H3K9ac (active promoter mark)
  - H3K27me3 (repressive mark)
- **Controls:** Input DNA
- **Replicates:** Biological replicates for each mark
- **Used for:** Peak calling, IDR analysis, visualization, and annotation

## Contact

**Raja Ishaq Nabi Khan**

- GitHub: [@ishaaq34](https://github.com/ishaaq34)
- Research Focus: RNA Modifications & Epigenetics

## Citation

If you use this tutorial in your research or teaching, please cite:

```
Khan, R. I. N. (2025). ChIP-seq Analysis Tutorial: A Practical Guide to 
Chromatin Immunoprecipitation Sequencing Data Analysis (Version 1.0). 
GitHub. https://github.com/ishaaq34/Chipseq_analysis_tutorial
```

---

*Last updated: December 2025 | Version 1.0*


