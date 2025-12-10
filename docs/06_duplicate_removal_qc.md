# 05. Handling Duplicates & Quality Control

## 1. Basic Concept: The "Photocopier" Analogy

Imagine you are trying to read a rare, handwritten manuscript (your DNA sample). You want to digitize it, so you take photos (sequencing reads) of different pages.

*   **Real Signal (Enriched Regions):** If many people are taking photos of the *same important page* because it's interesting, that's good! In ChIP-seq, this happens when a protein binds strongly to a specific DNA spot. We see many reads there because the biological signal is strong.
*   **Duplicates (Artifacts):** Now, imagine the photocopier gets stuck and prints 100 copies of a *random, unimportant page* just because of a machine error. These copies don't mean that page is 100 times more important; they are just **junk**. In sequencing, this is called **PCR duplication**—where the chemistry accidentally over-copies a single DNA fragment.

**Goal:** We want to keep the "popular pages" (real biological signal) but throw away the "accidental machine copies" (PCR duplicates) so they don't trick us into thinking a random spot is important.

---

## 2. Execution: Marking & Removing Duplicates

We use a tool called **Picard** to find these "accidental copies." It looks for reads that start and end at the exact same position. If it finds 10 identical reads, it keeps the best one and marks the other 9 as duplicates.

### Step 3.3: Mark Duplicates (Be Careful!)

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

### Step 3.4: Remove Duplicates (Clean Up)

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

## 3. Understanding the Details

For those who want to understand the "under the hood" mechanics, here is how Picard distinguishes different types of duplicates and why Read Groups matter.

### PCR vs. Optical Duplicates

Picard is smart enough to tell *why* a duplicate happened:

| Type | Origin | Typical Cause | Implication |
|------|---------|----------------|-------------|
| **PCR Duplicates** | Library Prep | Over-cycling (amplifying DNA too much) | Shows low library complexity (not enough unique DNA to start with). |
| **Optical Duplicates** | Sequencer | Camera errors reads the same cluster as two | Technical glitch on the flow cell. |

### PCR Duplicates vs. Biological Signal

How do we know it's a duplicate and not a strong signal?
*   **Duplicates:** Reads start and end at the *exact* same base pair.
*   **Real Signal:** Reads pile up in the same *region* (peak), but they have slightly different start/end points (staggered).

### The Role of Read Groups (RG)

A **Read Group** is a tag that tells the software "this read came from Sample A, Run 1."
*   **Why is it vital?** If you merge two different samples (e.g., Replicate 1 and Replicate 2), you might have two different reads that coincidentally map to the same spot.
*   Without Read Groups, Picard acts blindly: "These look identical! Delete one!" -> **Data Loss.**
*   With Read Groups, Picard sees: "Oh, one is from Replicate 1 and one is from Replicate 2. They are different samples. Keep both!"

**Picard's Superpower:** It checks for duplicates *within* a Read Group, never *across* them.

### References
*   [Biostars: Duplicate marking](https://www.biostars.org/p/318974/)
*   [GATK: Read Groups](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups)
*   [Picard Documentation](https://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups)
