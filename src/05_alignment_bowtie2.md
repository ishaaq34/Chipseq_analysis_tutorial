# 04: Alignment (Solving the Jigsaw Puzzle)

## Basic Concept (The "Puzzle")

Imagine your Reference Genome is the **Picture on the Puzzle Box** (a complete image of the DNA).
Your Reads (FASTAS) are the millions of tiny **Puzzle Pieces** scattered on the floor.

**Alignment** is simply picking up every piece and finding exactly where it fits on the picture.
*   **The Input:** Millions of jumbled reads.
*   **The Tool:** **Bowtie2** (The Puzzle Solver).
*   **The Output:** A **BAM File**. This is the digital record of where every piece belongs.

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
 

## After Indexing (Bowtie2)

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

**Why Sort?** A puzzle is useless if the pieces are in random order. We sort them by chromosome location (Left to Right) so we can look at them later. And we are working on clean pieces **(fastq files in the fastq_cleaned folder)** for alignment 

**Run Bowtie2 alignment for paired end (R1 and R2) sample**
```bash
mkdir -p bowalign

bowtie2 -x  genome_index/ce_index  \                  # 1. The Reference Map
  -1 fastq_cleaned/Sample1_R1.clean.fastq.gz \       # 2. Input Read 1 assume - SRR7297994_R1.clean.fastq.gz
  -2 fastq_cleaned/Sample1_R2.clean.fastq.gz \       # 3. Input Read 2 - SRR7297994_R2.clean.fastq.gz
  -p 6 --no-unal \                       # 4. Use 6 threads,  `--no-unal` — suppresses unaligned reads in the output.
  2> bowalign/Sample1.log |              # 5. Save the log file...
  samtools sort -@ 6 -o bowalign/Sample1.sorted.bam  # 6. ...and sort the result!
```

**Run Bowtie2 alignment for a single-end sample**

```
# Run Bowtie2 alignment for a single sample
bowtie2 -x genome_index/ce_index \
  -u fastq_cleaned/SRR7297994.clean.fastq.gz \
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

bowalign/
 ├── SRR7297994.log
 ├── SRR7297994.sorted.bam
 └── SRR7297994.sorted.bam.bai
```

Once this single run completes successfully, you can confidently automate for all samples.


### Step 3: Automation Loop 

Through [Sample list section](https://github.com/ishaaq34/Chipseq_analysis_tutorial/blob/main/src/03_sample_list_creation.md#ready-for-use) ,  here is the script to run this for all your samples:

```
#!/bin/bash
set -euo pipefail

mkdir -p bowalign bowalign_log


while read -r sample; do
  echo "Aligning $sample"

  bowtie2 -x genome_index/ce_index \
    -1 fastq_cleaned/${sample}_R1.clean.fq.gz" \
    -2 fastq_cleaned/${sample}_R2.clean.fq.gz" \
    -p 6 --no-unal \
    2> bowalign_log/${sample}.bowtie2.log"
  | samtools sort -@ 6 -o bowalign/${sample}.sorted.bam"

  samtools index "${OUTDIR}/${sample}.sorted.bam"

done < sample_id.txt

```

---

## Optimization 

### Optimization: Threads vs. Jobs
You have a limited number of CPU cores (computers brains). You can use them in two ways:

1.  **Multi-Threading (`-p 6`):** One sample uses 6 cores. It finishes very fast, but you only do **one sample at a time**.
    *   *Best for:* Large genomes, low memory.
2.  **Parallel Jobs:** You run **3 samples** at once, and each sample uses **2 cores**.
    *   *Best for:* Many small samples (RNA-seq, small genomes).

**Rule of Thumb:** `bowtie2` stops getting faster after about **8 threads**. Don't give it 50 threads; it's a waste!

 [Benchmarking Bowtie2 Threading - Jeff Kaufman (2023)](https://www.jefftk.com/p/benchmarking-bowtie2-threading)

 [BOWTIE2 - HPCC Wiki](https://wiki.csi.cuny.edu/HPCCWiki/BOWTIE2)

 [Guidance with using multiple threads with samtools - GitHub](https://github.com/samtools/samtools/issues)


**Example with GNU Parallel**:

```
#!/bin/bash
set -euo pipefail
mkdir -p bowalign bowalign_log

cat sample_id.txt | parallel -j 3 '                     # -j 3 → runs 3 samples (jobs) in parallel

  bowtie2 -x genome_index/ce_index \
    -1 fastq_cleaned/${sample}_R1.clean.fq.gz  \
    -2 fastq_cleaned/${sample}_R2.clean.fq.gz  \
    -p 4 --no-unal \                                   # -p 4 → bowtie2 uses 4 CPU threads per sample
    2>  bowalign_log/{}.log | samtools sort -@ 2 \           # -@ 2 → samtools sort uses 2 CPU threads per sample
    -o bowalign/{}.sorted.bam

  samtools index bowalign/{}.sorted.bam                 # samtools index is single-threaded
'
```


Effective CPU usage (implied by this setup)

- Per sample:  4 threads (bowtie2) + 2 threads (samtools sort) = 6 threads

- Total at once: 3 parallel jobs × 6 threads = 18 CPU threads

- Adjust based on available CPU cores and memory.


## After Alignment (Bowtie2)

```text
chipseq_tutorial/
├── fastq_fastq_raw/
│   └── ...
├── fastq_fastq_cleaned/
│   └── ...
├── fastp_reports/
│   └── ...
├── genome_index/
│   └── ...
├── bowalign/                ← Aligned BAM files
│   ├── SRR7297994.sorted.bam
│   ├── SRR7297994.sorted.bam.bai
│   ├── SRR7297995.sorted.bam
│   ├── SRR7297995.sorted.bam.bai
│   ├── SRR7297998.sorted.bam
│   ├── SRR7297998.sorted.bam.bai
│   ├── SRR7298003.sorted.bam
│   ├── SRR7298003.sorted.bam.bai
│   └── ...
├── bowalign_log/            ← Alignment statistics
│   ├── SRR7297994.log
│   ├── SRR7297995.log
│   ├── SRR7297998.log
│   ├── SRR7298003.log
│   └── ...
└── sample_id.txt
```







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




