# 05. Handling Duplicates & Quality Control

## 1. Basic Concept: The "Photocopier" Analogy

Imagine you are trying to read a rare, handwritten manuscript (your DNA sample). You want to digitize it, so you take photos (sequencing reads) of different pages.

*   **Real Signal (Enriched Regions):** If many people are taking photos of the *same important page* because it's interesting, that's good! In ChIP-seq, this happens when a protein binds strongly to a specific DNA spot. We see many reads there because the biological signal is strong.
*   **Duplicates (Artifacts):** Now, imagine the photocopier gets stuck and prints 100 copies of a *random, unimportant page* just because of a machine error. These copies don't mean that page is 100 times more important; they are just **junk**. In sequencing, this is called **PCR duplication**—where the chemistry accidentally over-copies a single DNA fragment.

**Goal:** We want to keep the "popular pages" (real biological signal) but throw away the "accidental machine copies" (PCR duplicates) so they don't trick us into thinking a random spot is important.

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
*   **Why is it vital?** If you merge two different samples (e.g., Replicate 1 and Replicate 2), you might have two different reads that coincidentally map to the same spot.
*   Without Read Groups, Picard acts blindly: "These look identical! Delete one!" -> **Data Loss.**
*   With Read Groups, Picard sees: "Oh, one is from Replicate 1 and one is from Replicate 2. They are different samples. Keep both!"




***Full hierarchy in RGing**

```
RGSM  → biological sample (Control H3k9ac)
  └── RGLB → library prep (usually lib1 however if Different library preps (even from same sample) RGLB=lib1 RGLB=lib2 )
        └── RGPU → flowcell + lane (from header of Fastq file )
              └── RGID → unique ID tying it all together (replicates of Control H3k9ac: Control H3k9ac_R1, ControlH3k9ac_R2 )


```

##  3. Marking & Removing Duplicates: Why Picard is Unique

**Picard** stands out among duplicate-marking tools because it:

- Distinguishes **optical vs. PCR duplicates** during QC reporting.  
- Leverages **read group (RG) identifiers** to avoid over-collapsing reads when combining multiple replicates or lanes.  
- Provides metrics for each library or sample, enabling fine-grained quality control.

Picard’s RG-aware logic ensures duplicates are flagged **within**, but not **across**, biological replicates or lanes — preventing false duplicate marking when BAMs are merged.

---
# 3.1 Adding RGs to each Bam file 

**Minimal augmented RG setup (merge BAM replicates)**

This read-group configuration is designed for the simplest defensible case: merging biological or technical ChIP-seq replicates and proceeding directly to peak calling.

```
picard AddOrReplaceReadGroups \
  I=H3K9ac_R1.bam \                 # Input BAM file (aligned, typically coordinate-sorted)
  O=H3K9ac_R1.RG.bam \              # Output BAM with read-group (RG) tags added/replaced
  RGID=H3K9ac_R1 \                  # Read Group ID: corresponds to: one replicate , If a replicate spans multiple lanes, you will have multiple RGIDs per replicate.
  RGSM=H39kac \                     # biological sample

```




**Minimal augmented RG setup (optical duplicates enabled)**

This configuration extends the previous one by adding the minimum required metadata to make optical duplicate detection meaningful. 
The addition of RGPU, encoding the flowcell and lane, defines the physical neighborhood in which optical duplicates can occur.


```
picard AddOrReplaceReadGroups \
  I=H3K9ac_R1.bam \
  O=H3K9ac_R1.RG.bam \
  RGID=H3K9ac_R1 \
  RGSM=H3K9ac \
  RGPL=ILLUMINA \                  # REQUIRED for optical duplicate logic
  RGPU=CA0TUACXX.1                 # REQUIRED: flowcell.lane (from FASTQ header)
```

More from [GATK: Read Groups](https://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups)
           [Picard markduplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard)

 


### Step 3.2: Mark Duplicates (Be Careful!)

First, we will just **mark** the duplicates but **keep them** in the file. This is like highlighting the duplicate pages in yellow but not throwing them in the trash yet. This is important for Quality Control (QC) to see how bad the duplication problem is.

```bash
# ----- 3.3 Mark duplicates (keep all reads, just mark) -----
# Create a folder for results
mkdir -p marked_reads

# Run Picard MarkDuplicates
picard MarkDuplicates \
  I=sample01.bam \
  O=marked_reads/sample01.marked.bam \
  M=marked_reads/metrics.txt \
  REMOVE_DUPLICATES=false
```

*   `I="$bam"`: Input BAM file.
*   `O=...`: Output BAM file (with duplicates marked).
*   `M=...`: Metrics file (report card showing how many duplicates were found).
*   `REMOVE_DUPLICATES=false`: Crucial! identifying them, not deleting them.

### Step 3.3: Remove Duplicates (Clean Up)

Now that we've checked the quality, we create a "clean" version of our data for analysis. We remove the marked duplicates so they don't interfere with peak calling.

```bash
# ----- 3.4 Create a duplicate-removed BAM by filtering marked BAM -----
# Create a folder for final results
mkdir -p final_clean_bam

# Run Picard to REMOVE duplicates
picard MarkDuplicates \
  I=sample01.bam \
  O=final_clean_bam/sample01.dedup.bam \
  M=final_clean_bam/dedup_metrics.txt \
  REMOVE_DUPLICATES=true

# Index the new file
samtools index final_clean_bam/sample01.dedup.bam
echo "Duplicates removed. Clean file is ready."
```

*   `REMOVE_DUPLICATES=true`: This time, we actually delete the highlighted duplicates.

---

## 4. Samtools

Samtools provides a lightweight and reliable way to handle duplicates without relying on read groups. The following steps use samtools to correctly prepare paired-end reads and mark or remove duplicate fragments.

```
# Paired-end mates are placed next to each other. It Groups reads by read name
samtools collate -o namecollate.bam example.bam

# it operates on the adjacent read pairs to synchronize mate flags , compute fragment-level information and add mate-related tags

samtools fixmate -m namecollate.bam fixmate.bam

# Coordinate sort for duplicate marking

samtools sort -o positionsort.bam fixmate.bam

# Mark or remove duplicates (-r after the markdup ) 
samtools markdup positionsort.bam markdup.bam

```

This samtools-based workflow is simpler than Picard and avoids the need for read-group metadata, while still providing accurate duplicate detection for paired-end data. It is well suited for ChIP-seq and other enrichment-based analyses where straightforward duplicate handling is sufficient.


https://www.htslib.org/algorithms/duplicate.html
