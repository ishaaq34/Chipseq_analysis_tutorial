# Tutorial 04: Alignment (Solving the Jigsaw Puzzle)

## Level 1: Basic Concept (The "Puzzle")

Imagine your Reference Genome is the **Picture on the Puzzle Box** (a complete image of the DNA).
Your Reads (FASTAS) are the millions of tiny **Puzzle Pieces** scattered on the floor.

**Alignment** is simply picking up every piece and finding exactly where it fits on the picture.
*   **The Input:** Millions of jumbled reads.
*   **The Tool:** **Bowtie2** (The Puzzle Solver).
*   **The Output:** A **BAM File**. This is the digital record of where every piece belongs.

---

## Level 2: Execution (Solving It)

### Step 1: The Index (Building the Map)
Before Bowtie2 can work, it needs to process the reference genome into a format it can search quickly. This is called an **Index**.

```bash
# Example syntax: bowtie2-build [genome.fa] [prefix_name]
bowtie2-build hg38.fa hg38_index
```
*You only do this once!*

### Step 2: Single Sample Alignment

We use a "Pipe" (`|`) to connect two tools: `bowtie2` aligns the data, and `samtools` sorts it immediately.

**Why Sort?** A puzzle is useless if the pieces are in random order. We sort them by chromosome location (Left to Right) so we can look at them later.

**Run Bowtie2 alignment for paired end (R1 and R2) sample**
```bash
mkdir -p bowalign

bowtie2 -x hg38_index \                  # 1. The Reference Map
  -1 trim/Sample1_R1.clean.fq.gz \       # 2. Input Read 1
  -2 trim/Sample1_R2.clean.fq.gz \       # 3. Input Read 2
  -p 6 --no-unal \                       # 4. Use 6 threads,  `--no-unal` — suppresses unaligned reads in the output.
  2> bowalign/Sample1.log |              # 5. Save the log file...
  samtools sort -@ 6 -o bowalign/Sample1.sorted.bam  # 6. ...and sort the result!
```

**Run Bowtie2 alignment for a single sample**

```
# Run Bowtie2 alignment for a single sample
bowtie2 -x trim/ssindex \
  -u trim/Sample1.clean.fq.gz \
  -p 6 --no-unal \
  2> bowalign/Sample1.log | samtools sort -@ 6 -o bowalign/Sample1.sorted.bam

```

**Final Touch:**
We always "index" the BAM file. Think of this as creating a **Table of Contents** so the computer can jump to any chromosome instantly.
```bash
samtools index bowalign/Sample1.sorted.bam
```



### **Output Structure**
After running this step, your directory should look like:
```
trim/
 ├── ssindex.1.bt2
 ├── ssindex.2.bt2
 └── ... (Bowtie2 index files)

bowalign/
 ├── Sample1.log
 ├── Sample1.sorted.bam
 └── Sample1.sorted.bam.bai
```

Once this single run completes successfully, you can confidently automate for all samples.


### Step 3: Automation Loop 

Through [Sample list section](https://github.com/ishaaq34/Chipseq_analysis_tutorial/blob/main/src/03_sample_list_creation.md#ready-for-use)  Here is the script to run this for all your samples:

```bash
#!/bin/bash
mkdir -p bowalign

# Loop through every ID in our "Attendance Sheet"
for sample in $(cat sample_id.txt); do
  echo "Aligning $sample..."
  
  bowtie2 -x hg38_index \
    -1 trim/${sample}_R1.clean.fq.gz \
    -2 trim/${sample}_R2.clean.fq.gz \
    -p 6 --no-unal \
    2> bowalign/${sample}.log | samtools sort -@ 6 -o bowalign/${sample}.sorted.bam

  samtools index bowalign/${sample}.sorted.bam
done
```

---

## Level 3: Optimization & QC

### Optimization: Threads vs. Jobs
You have a limited number of CPU cores (computers brains). You can use them in two ways:

1.  **Multi-Threading (`-p 6`):** One sample uses 6 cores. It finishes very fast, but you only do **one sample at a time**.
    *   *Best for:* Large genomes, low memory.
2.  **Parallel Jobs:** You run **3 samples** at once, and each sample uses **2 cores**.
    *   *Best for:* Many small samples (RNA-seq, small genomes).

**Rule of Thumb:** `bowtie2` stops getting faster after about **8 threads**. Don't give it 50 threads; it's a waste!

### Quality Check: `samtools flagstat`
Did the alignment work? Let's check the score.

```bash
samtools flagstat bowalign/Sample1.sorted.bam
```

**What to look for:**
*   **Mapping Rate:** Ideally **>80-90%**. If it's <50%, your DNA might be contaminated (e.g., bacteria in a human sample).
*   **Properly Paired:** Ideally **high**. This means R1 and R2 pointed towards each other at the correct distance.

---

## Summary
1.  **Analogy:** Alignment is placing puzzle pieces onto the reference picture.
2.  **Action:** Use `bowtie2` to align and `samtools sort` to organize.
3.  **Result:** A **Sorted BAM** file (the solved puzzle), ready for peak calling.
