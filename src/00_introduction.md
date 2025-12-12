# Welcome to the Practical ChIP-seq Tutorial

**From Raw Reads to Biological Insight: A "Zero-to-Hero" Guide.**

This book is a comprehensive, hands-on course designed to teach you **Chromatin Immunoprecipitation sequencing (ChIP-seq)** analysis. We replicate a full-scale **ENCODE** analysis, processing real biological data to uncover gene regulation mechanisms.

---

## Why This Tutorial?

We built this course to solve the common frustrations of bioinformatics learning. Here is what makes it different:

### 1. The "Tiered" Learning Method
We believe you shouldn't just run code—you should understand it. Every chapter is broken into three distinct levels:
*   **Level 1: Basic Concept (The "Why"):** We use simple English and real-world analogies (like "Jigsaw Puzzles" for Alignment) to explain the biology.
*   **Level 2: Execution (The "How"):** We provide the exact code, line-by-line, to run on your machine.
*   **Level 3: Interpretation (The "So What?"):** We teach you how to read the output. A plot is useless if you don't know what a "good" result looks like.


### 2. A Complete Pipeline
We don't skip steps. You will go from:
1.  **Downloading** raw FASTQ files from GEO.
2.  **Aligning** reads with Bowtie2.
3.  **Calling Peaks** with MACS2.
4.  **Visualizing** data with deepTools.
5.  **Annotating** genomic features with ChIPseeker.

---

### Experimental Context
These experiments were performed using the human **BLaER1** cell line. The cells were treated for 18 hours with three molecules:
*   **100 nM 17β-estradiol**
*   **10 ng/mL Interleukin-3**
*   **10 ng/mL CSF1**

These treatments activate specific signaling pathways in the cells and create a controlled biological state for studying transcription factors and histone modifications. The biosamples are described as **isogenic**, which means all replicates come from genetically identical cells and differ only in the experimental treatment or assay. The work was carried out in the laboratory of **Roderic Guigó** at CRG.

**Data Availability:**
*   [ENCODE Cart](https://www.encodeproject.org/carts/ca521f95-7835-4369-88a9-b89f98fb39ad/)

### Study Design: Targets & Controls
To study gene regulation under these conditions, several ChIP-seq experiments were carried out. These include:
*   **Transcription Factor ChIP-seq:** Targeted **CEBPA** to indicate specific transcription factor activity.
*   **Histone Modification ChIP-seq:** Targeted **H3K27me3** (repressed chromatin) and **H3K9ac** (active regulatory regions).

Each of these marks highlights a different aspect of gene regulation. Alongside these, **matching input (control) samples** were generated to measure background DNA fragmentation and sequencing noise. These input libraries are essential because they allow a clear comparison between true enrichment (ChIP signal) and nonspecific background.

### Tutorial Learning Goals
While the biological context is interesting, our primary goal is to **master the tools**. Instead of just looking for biological answers, we will focus on the technical questions that every analyst must answer:

*   **Quality Control:** How do we know if our sequencing data is good enough to use? Why are duplicates bad?
*   **Signal vs. Noise:** How precisely does the "Input" control help us filter out background artifacts?
*   **Peak Characteristics:** How do the calling algorithms differ when looking for a **sharp** peak (CEBPA) versus a **broad** mark (H3K27me3)?

By answering these technical questions, you will build a robust pipeline that can be applied to *any* organism or biological question.

By the end of this book, you will have the skills—and the code—to answer these questions for your own research.

**Let's get started**
