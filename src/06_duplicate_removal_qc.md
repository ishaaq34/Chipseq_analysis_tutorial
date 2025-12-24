# Handling Duplicates & Quality Control

`PCR-duplicates` `optical-duplicates` `Picard` `samtools` `MarkDuplicates` `deduplication` `read-groups` `MultiQC`

## 1. Basic Concept: The "Photocopier" Analogy

Imagine you are trying to read a rare, handwritten manuscript (your DNA sample). You want to digitize it, so you take photos (sequencing reads) of different pages.

* **Real Signal (Enriched Regions):** If many people are taking photos of the *same important page* because it's interesting, that's good! In ChIP-seq, this happens when a protein binds strongly to a specific DNA spot. We see many reads there because the biological signal is strong.
* **Duplicates (Artifacts):** Now, imagine the photocopier gets stuck and prints 100 copies of a *random, unimportant page* just because of a machine error. These copies don't mean that page is 100 times more important; they are just **junk**. In sequencing, this is called **PCR duplication**—where the chemistry accidentally over-copies a single DNA fragment.

**Goal:** We want to keep the "popular pages" (real biological signal) but throw away the "accidental machine copies" (PCR duplicates) so they don't trick us into thinking a random spot is important.

> [!NOTE]
> **Two Tools, Same Goal:** This chapter covers two methods for duplicate removal. **Picard** (Section 3) offers fine-grained control with read groups and optical duplicate detection—ideal for complex multi-replicate studies. **Samtools** (Section 4) provides a simpler, single-pipeline approach—sufficient for most ChIP-seq analyses. Choose the method that fits your workflow.

---

## 2. Understanding the Details

For those who want to understand the "under the hood" mechanics, here is what are different types of duplicates and  the Read Groups .

### PCR vs. Optical Duplicates

| Type | Origin | Typical Cause | Implication |
|------|---------|----------------|-------------|
| **PCR Duplicates** | Library Prep | Over-cycling (amplifying DNA too much) | Shows low library complexity (not enough unique DNA to start with). |
| **Optical Duplicates** | Sequencer | Camera errors reads the same cluster as two | Technical glitch on the flow cell. |

### The Role of Read Groups (RG)

A **Read Group** is a tag that tells the software "this read came from Sample A, Run 1."

* **Why is it vital?** If you merge two different samples (e.g., Replicate 1 and Replicate 2), you might have two different reads that coincidentally map to the same spot.
* Without Read Groups, Picard acts blindly: "These look identical! Delete one!" -> **Data Loss.**
* With Read Groups, Picard sees: "Oh, one is from Replicate 1 and one is from Replicate 2. They are different samples. Keep both!"

***Full hierarchy in RGing**

```text
RGSM  → biological sample (Control H3k9ac)
  └── RGLB → library prep (usually lib1 however if Different library preps (even from same sample) RGLB=lib1 RGLB=lib2 )
        └── RGPU → flowcell + lane (from header of Fastq file )
              └── RGID → unique ID tying it all together (replicates: H3K27me3_IP_rep1, H3K27me3_IP_rep2)


```

## 3. Marking & Removing Duplicates: Why Picard is Unique

**Picard** stands out among duplicate-marking tools because it:

* Distinguishes **optical vs. PCR duplicates** during QC reporting.  
* Leverages **read group (RG) identifiers** to avoid over-collapsing reads when combining multiple replicates or lanes.  
* Provides metrics for each library or sample, enabling fine-grained quality control.

Picard’s RG-aware logic ensures duplicates are flagged **within**, but not **across**, biological replicates or lanes — preventing false duplicate marking when BAMs are merged.

---

# 3.1 Adding RGs to each Bam file

**Minimal augmented RG setup (merge BAM replicates)**

This read-group configuration is designed for the simplest defensible case: merging biological or technical ChIP-seq replicates and proceeding directly to peak calling.

**Input files needed:**

* Filtered BAM file: `bowalign_filtered/H3K27me3_IP_rep1.filtered.bam` (from section 05)

```bash
mkdir -p picard_rg_bam

picard AddOrReplaceReadGroups \
  I=bowalign_filtered/H3K27me3_IP_rep1.filtered.bam \
  O=picard_rg_bam/H3K27me3_IP_rep1.RG.bam \
  RGID=H3K27me3_IP_rep1 \
  RGSM=H3K9ac
```

**What this does:**

* **I=** specifies the input BAM file from the alignment step
* **O=** specifies the output BAM file with read-group tags added
* **RGID=** sets a unique Read Group ID (one per replicate)
* **RGSM=** sets the biological sample name (used for grouping replicates)

**Minimal augmented RG setup (optical duplicates enabled)**

This configuration extends the previous one by adding the minimum required metadata to make optical duplicate detection meaningful.
The addition of RGPU, encoding the flowcell and lane, defines the physical neighborhood in which optical duplicates can occur.

```bash
picard AddOrReplaceReadGroups \
  I=bowalign_filtered/H3K27me3_IP_rep1.filtered.bam \     # Input: MAPQ-filtered BAM from section 05
  O=picard_rg_bam/H3K27me3_IP_rep1.RG.bam \           # Output with RG tags
  RGID=H3K27me3_IP_rep1 \                      # Read Group ID
  RGSM=H3K9ac \                          # Biological sample
  RGPL=ILLUMINA \                        # REQUIRED for optical duplicate logic
  RGPU=CA0TUACXX.1                       # REQUIRED: flowcell.lane (from FASTQ header)
```

More from [GATK: Read Groups](https://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups)
           [Picard markduplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard)

### Step 3.2: Mark Duplicates (Be Careful!)

First, we will just **mark** the duplicates but **keep them** in the file. This is like highlighting the duplicate pages in yellow but not throwing them in the trash yet. This is important for Quality Control (QC) to see how bad the duplication problem is.

**Input files needed:**

* Filtered BAM file with read groups: `picard_rg_bam/H3K27me3_IP_rep1.RG.bam`

```bash
# ----- 3.3 Mark duplicates (keep all reads, just mark) -----

mkdir -p picard_markdup picard_markdup_metrics

# Run Picard MarkDuplicates
picard MarkDuplicates \
  I=picard_rg_bam/H3K27me3_IP_rep1.RG.bam \
  O=picard_markdup/H3K27me3_IP_rep1.marked.bam \
  M=picard_markdup_metrics/H3K27me3_IP_rep1.metrics.txt \
  REMOVE_DUPLICATES=false

# Index the new file
samtools index picard_markdup/H3K27me3_IP_rep1.marked.bam
```

**What this does:**

* **I=** specifies input BAM file with read groups
* **O=** specifies output BAM file with duplicates marked (not removed)
* **M=** creates metrics file reporting duplicate statistics
* **REMOVE_DUPLICATES=false** marks duplicates but keeps all reads in the file
* **samtools index** creates index for marked BAM file

> [!IMPORTANT]
> **Why keep marked files?** The Picard metrics file (`.metrics.txt`) contains duplicate rates that **MultiQC** can aggregate into summary reports. However, for **downstream analysis (peak calling, BigWig generation)**, always use the **deduplicated BAMs** to prevent PCR artifacts from inflating your signal.
>
> *Note:* While MACS can handle duplicates internally via `--keep-dup` options, using pre-deduplicated BAMs is the recommended practice—it gives you explicit control and avoids relying on MACS's position-based duplicate detection, which differs from Picard/samtools fragment-level detection. [[GitHub #33](https://github.com/macs3-project/MACS/issues/33)] [[Biostars](https://www.biostars.org/p/318974/#319028)]

**Aggregate metrics with MultiQC:**

```bash
# Run MultiQC on Picard metrics folder
multiqc picard_dedup_metrics/ -o multiqc_report/
```

This generates an interactive HTML report showing duplicate rates across all samples.

### Step 3.3: Remove Duplicates (Clean Up)

Now that we've checked the quality, we create a "clean" version of our data for analysis. We remove the marked duplicates so they don't interfere with peak calling.

```bash
# ----- 3.4 Create a duplicate-removed BAM by filtering marked BAM -----

mkdir -p picard_dedup_bam picard_dedup_metrics

# Run Picard to REMOVE duplicates
picard MarkDuplicates \
  I=picard_rg_bam/H3K27me3_IP_rep1.RG.bam \
  O=picard_dedup_bam/H3K27me3_IP_rep1.dedup.bam \
  M=picard_dedup_metrics/H3K27me3_IP_rep1.metrics.txt \
  REMOVE_DUPLICATES=true

# Index the new file
samtools index picard_dedup_bam/H3K27me3_IP_rep1.dedup.bam
```

* `REMOVE_DUPLICATES=true`: This time, we actually delete the highlighted duplicates.

### After Picard Deduplication

```text
chipseq_tutorial/
├── fastq_raw/
├── fastq_cleaned/
├── genome_index/
├── bowalign/                    ← Aligned BAMs
├── bowalign_filtered/           ← Filtered BAMs
├── picard_rg_bam/               ← BAMs with Read Groups
│   ├── H3K27me3_IP_rep1.RG.bam
│   ├── H3K27me3_IP_rep1.RG.bam.bai
│   └── ...
├── picard_dedup_bam/            ← Final Deduplicated BAMs
│   ├── H3K27me3_IP_rep1.dedup.bam
│   ├── H3K27me3_IP_rep1.dedup.bam.bai
│   └── ...
├── picard_metrics/              ← Duplication Reports
│   ├── H3K27me3_IP_rep1.metrics.txt
│   └── ...
└── sample_id.txt
```

---

## 4. Samtools (Simple Alternative)

Samtools provides a lightweight and reliable way to handle duplicates without the complexity of read groups. Choose this option if you want a straightforward workflow.

### Step 4.1: Prepare BAM for Duplicate Detection

```bash
# Paired-end mates are placed next to each other. Groups reads by read name
samtools collate -o temp/H3K27me3_IP_rep1.collate.bam bowalign_filtered/H3K27me3_IP_rep1.filtered.bam

# Synchronize mate flags and add mate-related tags
samtools fixmate -m temp/H3K27me3_IP_rep1.collate.bam temp/H3K27me3_IP_rep1.fixmate.bam

# Coordinate sort for duplicate marking
samtools sort -o temp/H3K27me3_IP_rep1.positionsort.bam temp/H3K27me3_IP_rep1.fixmate.bam
```

### Step 4.2: Mark Duplicates (Keep for QC)

```bash
# Mark duplicates but KEEP them in the output
samtools markdup temp/H3K27me3_IP_rep1.positionsort.bam samtools_markdup/H3K27me3_IP_rep1.marked.bam

# Index the file
samtools index samtools_markdup/H3K27me3_IP_rep1.marked.bam
```

### Step 4.3: Remove Duplicates (Clean for Analysis)

```bash
# Remove duplicates with -r flag
samtools markdup -r temp/H3K27me3_IP_rep1.positionsort.bam samtools_dedup_bam/H3K27me3_IP_rep1.dedup.bam

# Index the file
samtools index samtools_dedup_bam/H3K27me3_IP_rep1.dedup.bam
```

This samtools-based workflow is simpler than Picard and avoids the need for read-group metadata, while still providing accurate duplicate detection for paired-end data.

### Step 4.4: Automation Loop (All Samples)

```bash
#!/bin/bash
set -euo pipefail

mkdir -p samtools_dedup_bam

while read -r sample; do
  echo "Processing $sample..."

  # All steps piped together - no temp files needed!
  samtools collate -u -O "bowalign_filtered/${sample}.filtered.bam" | \
    samtools fixmate -m -u - - | \
    samtools sort -u - | \
    samtools markdup -r - "samtools_dedup_bam/${sample}.dedup.bam"

  # Index the final output
  samtools index "samtools_dedup_bam/${sample}.dedup.bam"

  echo "Finished $sample"
done < sample_id.txt
```

**Why use pipes?**

```
collate → fixmate → sort → markdup
  ↓         ↓        ↓       ↓
 pipe     pipe     pipe    file
```

* **`-u` flag**: Outputs uncompressed BAM for faster streaming between commands
* **`-` symbol**: Represents stdin/stdout, allowing data to flow through the pipeline
* **No temp files**: Everything processes in memory, saving disk space and I/O time
* **Faster**: No waiting for disk writes between each step

More [samtools](https://www.htslib.org/algorithms/duplicate.html)

### After Samtools Deduplication

```text
chipseq_tutorial/
├── bowalign_filtered/           ← Input: MAPQ-filtered BAM files (from section 05)
│   ├── H3K27me3_IP_rep1.filtered.bam
│   └── ...
├── samtools_markdup/            ← Duplicates marked (not removed)
│   ├── H3K27me3_IP_rep1.marked.bam
│   └── ...
├── samtools_dedup_bam/          ← Duplicates removed
│   ├── H3K27me3_IP_rep1.dedup.bam
│   ├── H3K27me3_IP_rep1.dedup.bam.bai
│   └── ...
└── sample_id.txt
```

---

## Summary

1. **Understand:** PCR duplicates inflate signal artificially; we mark/remove them.
2. **Action:** Use Picard (with read groups) or samtools to handle duplicates.
3. **Result:** Clean, deduplicated BAM files ready for peak calling.

> [!NOTE]
> **Up Next:** We'll assess library complexity to understand how deeply our libraries were sequenced.
