# Alignment (Solving the Jigsaw Puzzle)

`Bowtie2` `alignment` `BAM-files` `single-end` `paired-end` `samtools` `genome-index` `MAPQ` `quality-control` `multimapping` `ChIP-seq` `flagstat` `MultiQC` `filtering`

## Basic Concept (The "Puzzle")

Imagine your Reference Genome is the **Picture on the Puzzle Box** (a complete image of the DNA).
Your Reads (FASTAS) are the millions of tiny **Puzzle Pieces** scattered on the floor.

**Alignment** is simply picking up every piece and finding exactly where it fits on the picture.

* **The Input:** Millions of jumbled reads.
* **The Tool:** **Bowtie2** (The Puzzle Solver).
* **The Output:** A **BAM File**. This is the digital record of where every piece belongs.

---

## Execution (Solving It)

### Step 1: The Index (Building the Map)

Before Bowtie2 can work, it needs to process the reference genome into a format it can search quickly. This is called an **Index**.

```bash
# Example syntax: bowtie2-build [genome.fa] [prefix_name]
bowtie2-build genome_index/ce.fa genome_index/ce_index            # ce.fa = C. elegans reference genome
```

*You only do this once!*

You can download the C. elegans reference genome from Ensembl using the following link:  [C. elegans genome](https://ftp.ensembl.org/pub/release-115/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz).  Place it in the **genome_index** directory, decompress it, and rename the file to `ce.fa`

### After Indexing (Bowtie2)

```text
chipseq_tutorial/
├── fastq_raw/
│   └── ...
├── fastq_cleaned/                ← Fastp cleaned reads
│   └── ...
├── fastp_reports/
│   └── ...
├── genome_index/           ← Bowtie2 index files
│   ├── ce_index.1.bt2
│   ├── ce_index.2.bt2
│   ├── ce_index.3.bt2
│   ├── ce_index.4.bt2
│   ├── ce_index.rev.1.bt2
│   └── ce_index.rev.2.bt2
└── sample_id.txt
```

### Step 2: Single Sample Alignment

We use a "Pipe" (`|`) to connect two tools: `bowtie2` aligns the data, and `samtools` sorts it immediately.

**Why Sort?** A puzzle is useless if the pieces are in random order. We sort them by chromosome location (Left to Right) so we can look at them later.

**Run Bowtie2 alignment for a single-end sample**

**Input files needed:**

* Reference genome index: `genome_index/ce_index` (created in previous step)
* Cleaned FASTQ file: `fastq_cleaned/H3K27me3_IP_rep1.clean.fastq.gz`

```bash
mkdir -p bowalign

bowtie2 -x genome_index/ce_index \
  -U fastq_cleaned/H3K27me3_IP_rep1.clean.fastq.gz \
  -p 6 --no-unal \
  2> bowalign/H3K27me3_IP_rep1.log | \
  samtools sort -@ 6 -o bowalign/H3K27me3_IP_rep1.sorted.bam

samtools index bowalign/H3K27me3_IP_rep1.sorted.bam
```

**What this does:**

1. **bowtie2** aligns reads to the reference genome using the `-U` flag for single-end data
2. **-p 6** uses 6 CPU threads for faster processing
3. **--no-unal** suppresses unaligned reads from output (reads that don't match the genome are excluded, saving disk space and processing time downstream)
4. **2>** saves alignment statistics to a log file
5. **samtools sort** sorts alignments by genomic position (required for downstream analysis)
6. **samtools index** creates an index file (.bai) for fast random access

> [!NOTE]
> **If you have paired-end data**, use `-1` and `-2` flags instead of `-U`:

```bash
bowtie2 -x genome_index/ce_index \
  -1 fastq_cleaned/Sample_R1.clean.fastq.gz \
  -2 fastq_cleaned/Sample_R2.clean.fastq.gz \
  -p 6 --no-unal \
  2> bowalign/Sample.log | samtools sort -@ 6 -o bowalign/Sample.sorted.bam

samtools index bowalign/Sample.sorted.bam
```

**Final Touch:**
We always "index" the BAM file. Think of this as creating a **Table of Contents** so the computer can jump to any chromosome instantly.

```bash
samtools index bowalign/H3K27me3_IP_rep1.sorted.bam
```

### **Output Structure**

After running this step, your directory should look like:

```text
bowalign/
├── H3K27me3_IP_rep1.log
├── H3K27me3_IP_rep1.sorted.bam
└── H3K27me3_IP_rep1.sorted.bam.bai
```

Once this single run completes successfully, you can confidently automate for all samples.

### Step 3: Automation Loop

After creating `sample_id.txt` (see [Step 2.4 in file 02](02_geo_fastq_download.md#24-creating-sample-id-list)), here is the script to run this for all your samples:

```bash
#!/bin/bash
set -euo pipefail

mkdir -p bowalign bowalign_log

while read -r sample; do
  echo "Aligning $sample"

  bowtie2 -x genome_index/ce_index \
    -U "fastq_cleaned/${sample}.clean.fastq.gz" \
    -p 6 --no-unal \
    2> "bowalign_log/${sample}.bowtie2.log" \
    | samtools sort -@ 6 -o "bowalign/${sample}.sorted.bam"

  samtools index "bowalign/${sample}.sorted.bam"

done < sample_id.txt
```

> [!NOTE]
> **Paired-end while loop** (if your samples have _R1 and_R2 files):

```bash
while read -r sample; do
  echo "Aligning $sample"

  bowtie2 -x genome_index/ce_index \
    -1 "fastq_cleaned/${sample}_R1.clean.fastq.gz" \
    -2 "fastq_cleaned/${sample}_R2.clean.fastq.gz" \
    -p 6 --no-unal \
    2> "bowalign_log/${sample}.log" | \
  samtools sort -@ 6 -o "bowalign/${sample}.sorted.bam"

  samtools index "bowalign/${sample}.sorted.bam"

done < sample_id.txt
```

---

## Optimization

### Optimization: Threads vs. Jobs

You have a limited number of CPU cores (computers brains). You can use them in two ways:

1. **Multi-Threading (`-p 6`):** One sample uses 6 cores. It finishes very fast, but you only do **one sample at a time**.
    * *Best for:* Large genomes, low memory.
2. **Parallel Jobs:** You run **3 samples** at once, and each sample uses **2 cores**.
    * *Best for:* Many small samples (RNA-seq, small genomes).

**Rule of Thumb:** `bowtie2` stops getting faster after about **8 threads**. Don't give it 50 threads; it's a waste!

[Benchmarking Bowtie2 Threading - Jeff Kaufman (2023)](https://www.jefftk.com/p/benchmarking-bowtie2-threading)

[BOWTIE2 - HPCC Wiki](https://wiki.csi.cuny.edu/HPCCWiki/BOWTIE2)

[Guidance with using multiple threads with samtools - GitHub](https://github.com/samtools/samtools/issues)

**Example with GNU Parallel**:

```bash
#!/bin/bash
set -euo pipefail

mkdir -p bowalign bowalign_log

parallel -j 2 '
  bowtie2 \
    -x genome_index/ce_index \
    -U fastq_cleaned/{}.clean.fastq.gz \
    -p 4 \
    --no-unal \
    2> bowalign_log/{}.log \
  | samtools sort \
      -@ 2 \
      -m 1G \
      -o bowalign/{}.sorted.bam

  samtools index bowalign/{}.sorted.bam
' :::: sample_id.txt
```

**Effective CPU usage (implied by this setup):**

* 2 parallel samples (`-j 2`)
* Per sample: bowtie2 (4 threads, `-p 4`) + samtools sort (2 threads, `-@ 2`)
* Total ≈ 12 threads → safe on a 16-core machine

```text
chipseq_tutorial/
├── fastq_raw/
├── fastq_cleaned/
├── fastp_reports/
├── genome_index/
├── bowalign/                ← Aligned BAM files
│   ├── H3K27me3_IP_rep1.sorted.bam
│   ├── H3K27me3_IP_rep1.sorted.bam.bai
│   └── ...
├── bowalign_log/            ← Alignment statistics
│   ├── H3K27me3_IP_rep1.log
│   └── ...
└── sample_id.txt
```

### Understanding BAM File Structure

Before we check alignment quality, let's look at what's inside a BAM file. BAM files store alignment information in a structured format with **11 mandatory columns**:

```bash
# View the first few alignments in human-readable format
samtools view bowalign/H3K27me3_IP_rep1.sorted.bam | head -3
```

**Example output (single-end data):**

```text
SRR7297994.1    0    I      15200    42    76M    *    0    0    ACGTACGT...    IIIIIIII...
SRR7297994.2    0    II     98432    30    76M    *    0    0    TGCATGCA...    IIIHHHHH...
SRR7297994.3   16    III    45123     0    76M    *    0    0    GATTACA...     IIIHHHGG...
```

**Column Breakdown (First 5 of 11 columns):**

| Column | Name | Example | Description |
|--------|------|---------|-------------|
| 1 | QNAME | SRR7297994.1 | Read/Query name |
| 2 | **FLAG** | 0, 16 | **Bitwise FLAG** (flagstat reads this!) |
| 3 | RNAME | I, II, III | Reference sequence name (chromosome) |
| 4 | POS | 15200 | Leftmost mapping position |
| 5 | **MAPQ** | 42, 30, 0 | **Mapping Quality score** (0-255) |

**What Each QC Tool Evaluates:**

**`samtools flagstat`** reads **Column 2 (FLAG)**:

* Bit 4: Read unmapped?
* Bit 256: Secondary alignment?
* Bit 1024: PCR/optical duplicate?
* Bit 2048: Supplementary alignment?
* And more... (see FLAG section below)

> [!TIP]
> **Decode FLAGS:** Use [Explain SAM Flags](https://broadinstitute.github.io/picard/explain-flags.html) to understand bitwise FLAG values.

### Multimapping & The "Lost GPS"

Sometimes a read is repetitive (e.g., "ATATATAT"). It fits in 50 different places on the genome.
The aligner doesn't know which spot is correct, so it gives it a **Low MAPQ Score** (Mapping Quality).
**MAPQ filtering** (`-q 30`) reads **Column 5 (MAPQ)**:

* **MAPQ = 0:** "I have NO clue where this goes. It fits in many places." (Multimapped)
* **MAPQ > 30:** "I am highly confident this read belongs EXACTLY here." (Unique)

### Quality Check: `samtools flagstat` & `samtools stats`

Did the alignment work? Let's check the score using two tools:

* **flagstat**: Quick alignment summary (mapped reads, pairs, etc.)
* **stats**: Comprehensive metrics including MAPQ distribution (for MultiQC visualization)

```bash
# Create QC output directory for alignment metrics
mkdir -p bowalign_qc

# Quick summary with flagstat
samtools flagstat bowalign/H3K27me3_IP_rep1.sorted.bam > bowalign_qc/H3K27me3_IP_rep1.flagstat.txt

# Comprehensive stats (includes MAPQ distribution for MultiQC)
samtools stats bowalign/H3K27me3_IP_rep1.sorted.bam > bowalign_qc/H3K27me3_IP_rep1.stats.txt
```

**Batch processing for all samples:**

```bash
#!/bin/bash
set -euo pipefail

mkdir -p bowalign_qc

while read -r sample; do
  echo "Running QC for: $sample"
  
  # Generate flagstat summary
  samtools flagstat bowalign/${sample}.sorted.bam > bowalign_qc/${sample}.flagstat.txt
  
  # Generate comprehensive stats (for MultiQC)
  samtools stats bowalign/${sample}.sorted.bam > bowalign_qc/${sample}.stats.txt
  
done < sample_id.txt

# Generate MultiQC report with all QC results
echo "Generating MultiQC report..."
multiqc bowalign_qc/ -o bowalign_qc/ -n alignment_qc_report
echo "MultiQC report saved: bowalign_qc/alignment_qc_report.html"
```

This creates an interactive HTML report (`alignment_qc_report.html`) with visualizations of alignment metrics across all samples, including MAPQ distribution plots.

---

**What to look for:**

* **Mapping Rate:** Ideally **>80-90%**. If it's <50%, your DNA might be contaminated (e.g., bacteria in a sample).
* **Goal:** High mapping percentage indicates good quality alignment.
* **Warning:** If mapping is low (<50%), you may have the wrong organism or bad sequencing.

**The Sensitivity-Specificity Trade-off:**

Multi-mapping reads critically influence the balance between sensitivity and specificity in ChIP-seq peak detection. Excluding multi-mappers, standard practice in most peak callers, improves specificity by preventing artificial signal inflation in repetitive regions but reduces sensitivity for genuine binding events within those regions ([ENCODE Project Consortium, 2012](https://www.encodeproject.org/)). This trade-off is acceptable for transcription factors and chromatin features predominantly in unique genomic loci. However, repetitive and transposable elements (REs/TEs) constitute significant genome portions with important regulatory roles, and discarding multi-mapped reads substantially underrepresents regulatory events in these regions ([He et al., 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/)). For studies involving repetitive elements or broad chromatin domains, alternative strategies include paired-end sequencing (reduces mapping ambiguity), probabilistic allocation methods like *Allo* and *MuMRescueLite* (allocate reads based on uniquely mapped read distributions), or consensus sequence mapping (family-level TE insights) ([Chung et al., 2011](https://academic.oup.com/); [He et al., 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/)). Caution is advised as improper use can lead to double-counting and false positives.

#### Quick Manual MAPQ Check

For a **quick manual inspection** of MAPQ distribution (complementary to the MultiQC visualization):

```bash
samtools view bowalign/H3K27me3_IP_rep1.sorted.bam | awk '{print $5}' | sort -n | uniq -c > bowalign_qc/H3K27me3_IP_rep1.mapq.txt
```

**Understanding the Manual Output:**

* **If you see lots of 0s:** Your file includes "Lost GPS" reads (Multimappers).
* **If your scores start at 30+:** Your file has **Already Filtered** the bad reads.

**Output format: `count` `MAPQ_score`**

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

### Filtering Multi-mapping Reads

For standard ChIP-seq analysis (TFs, histone marks in unique regions), filter BAM files to retain only high-quality uniquely mapped reads:

```bash
# Create directory for filtered BAM files
mkdir -p bowalign_filtered

# Filter out multi-mappers (keep only MAPQ ≥ 30)
samtools view -b -q 30 bowalign/H3K27me3_IP_rep1.sorted.bam > bowalign_filtered/H3K27me3_IP_rep1.filtered.bam

# Index the filtered BAM
samtools index bowalign_filtered/H3K27me3_IP_rep1.filtered.bam
```

**Batch processing for all samples:**

```bash
#!/bin/bash
set -euo pipefail

mkdir -p bowalign_filtered

while read -r sample; do
  echo "Filtering multi-mappers for: $sample"
  samtools view -b -q 30 bowalign/${sample}.sorted.bam > bowalign_filtered/${sample}.filtered.bam
  samtools index bowalign_filtered/${sample}.filtered.bam
done < sample_id.txt
```

> [!NOTE]
> The `-q 30` flag filters out reads with MAPQ < 30 (removes multi-mappers and low-confidence alignments). Most downstream tools (MACS2, deepTools) can also filter by MAPQ, so this step is optional but recommended for cleaner data.

---

### Directory Structure After Alignment QC

```text
chipseq_tutorial/
├── fastq_raw/
├── fastq_cleaned/
├── genome_index/
├── bowalign/                    ← Aligned BAM files
│   ├── H3K27me3_IP_rep1.sorted.bam
│   ├── H3K27me3_IP_rep1.sorted.bam.bai
│   └── ...
├── bowalign_log/                ← Alignment statistics
│   ├── H3K27me3_IP_rep1.log
│   └── ...
├── bowalign_qc/                 ← Quality control metrics
│   ├── H3K27me3_IP_rep1.flagstat.txt
│   ├── H3K27me3_IP_rep1.stats.txt
│   ├── H3K27me3_IP_rep1.mapq.txt
│   ├── alignment_qc_report.html
│   ├── alignment_qc_report_data/
│   └── ...
├── bowalign_filtered/           ← Filtered BAM files (MAPQ ≥ 30)
│   ├── H3K27me3_IP_rep1.filtered.bam
│   ├── H3K27me3_IP_rep1.filtered.bam.bai
│   └── ...
└── sample_id.txt
```

---

## Summary

1. **Analogy:** Alignment is placing puzzle pieces onto the reference picture.
2. **Action:** Use `bowtie2` to align and `samtools sort` to organize.
3. **Result:** A **Sorted BAM** file (the solved puzzle), ready for peak calling.

> [!NOTE]
> **Up Next:** Before peak calling, we need to remove PCR duplicates and assess library quality.
