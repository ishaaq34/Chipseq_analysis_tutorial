# Understanding and Cleaning Your Data (FASTQ & fastp)

`FASTQ-format` `fastp` `quality-control` `adapter-trimming` `Phred-scores` `read-filtering` `QC-reports` `quality-assessment`

## 1: Basic Concept (The Anatomy of a Read)

A FASTQ file is just a text file full of DNA sequences. But unlike a simple list of letters, every single read carries extra baggage (its quality score).

Think of every read like a **Luggage Tag** with 4 lines of information:

1. **Line 1 (The Header):** Starts with `@`. This is the **ID Card**. It tells you the machine name, flowcell lane, and coordinates.
2. **Line 2 (The Sequence):** The DNA letters (`ACTG...`). This is the **Content** inside the bag.
3. **Line 3 (The Spacer):** Starts with `+`. Just a divider.
4. **Line 4 (The Quality):** A string of weird characters (`F:F#,,...`). This is the **Trust Score**. Each character represents the probability that the corresponding base in Line 2 is wrong.

**Example Read:**

```text
@SN227:495:CA0TUACXX:1:1106:1159:2114 1:N:0:ATCACG  <-- ID: Read #1
GTAAAAAGATTACATATATATTTAAAGTACACTGTAATTCTTANCA    <-- DNA: "N" means the machine failed to call that base
+                                                  <-- Spacer
FDFFFHHHHHJIJGIJIJJHJJJJJJIIJIJGIHFIIJJIIIIJG        <-- Quality: Each character encodes a Phred quality score (Illumina 1.8+, ASCII offset 33).

```

---

## Level 2: Execution (The Car Wash)

Before we start analyzing, we need to clean our data.

* **The Problem:** Sequencers sometimes make mistakes, especially at the ends of reads. They also leave "adapters" (artificial tags) attached to the DNA.
* **The Solution:** We use a tool called **fastp**. It acts like an automatic car wash: dirty reads go in, clean reads come out.

### 2.1 Basic Cleaning (Single-End)

```bash
# -i: Input (dirty)
# -o: Output (clean)
fastp -i fastq_raw/H3K27me3_IP_rep1.fastq.gz -o fastq_cleaned/H3K27me3_IP_rep1.clean.fastq.gz
```

> [!NOTE]
> **This dataset uses single-end sequencing.** If you have paired-end data, you would use:
>
> ```bash
> fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz
> ```

**What does fastp do automatically?**

* **Quality Filtering:** Drops reads if too many bases have low scores (default: Phred -q < 15).
* **Length Filtering:** Drops reads that become too short after trimming (default: -l < 15bp).
* **Adapter Removal:** Finds and cuts off adapter sequences automatically.

---

## Level 3: Advanced Analysis (The Math)

### 3.1 Quick Stats with AWK

Sometimes you don't want to run a full report; you just want to know "How many reads do I have?"
You can use `awk` (a math tool for text) to count directly from the compressed file.

**Count Total Reads:**

```bash
# A FASTQ record is 4 lines. We count total lines and divide by 4.
gzcat fastq_raw/H3K27me3_IP_rep1.fastq.gz | wc -l | awk '{print $1/4 " reads"}'
```

**Count Total Bases (Coverage):**

```bash
# Sums the length of line 2 (sequence) for every record
gzcat fastq_raw/H3K27me3_IP_rep1.fastq.gz | awk 'NR%4==2 {b+=length($0)} END{print b " bases"}'
```

* *Approximation:* If you have 100 Million bases and your genome is 3 Billion bases (Human), your coverage is roughly 0.03x.

### 3.2 Batch Processing

If you have 50 files, you can use a script to run `fastp` on all of them in parallel.
The `fastp` developers provide a handy script called  [parallel.py](https://github.com/OpenGene/fastp/blob/master/parallel.py):

```bash
# Process 3 files at a time (-f 3), using 2 threads per file (-t 2)
python parallel.py -i /fastq_raw -o /fastq_cleaned -r /fastp_reports -f 3 -t 2
```

**Parameter explanation:**

* `-f 3`: Sets the batch size—the script processes 3 files at a time
* `-t 2`: Sets the number of threads per file—each file is processed with 2 threads

This automatically finds pairs and generates HTML reports for every sample.

---
> [!IMPORTANT]
> [parallel.py](https://github.com/OpenGene/fastp/blob/master/parallel.py) avoids the need to explicitly loop over `sample_id.txt` in a Bash script.

```text
chipseq_tutorial/
├── fastq_raw/                  # Original FASTQ files
│   ├── H3K27me3_IP_rep1.fastq.gz
│   ├── H3K27me3_IP_rep2.fastq.gz
│   └── Input_rep1.fastq.gz
├── fastq_cleaned/              # Cleaned FASTQ files
│   ├── H3K27me3_IP_rep1.clean.fastq.gz
│   ├── H3K27me3_IP_rep2.clean.fastq.gz
│   └── Input_rep1.clean.fastq.gz
├── fastp_reports/              # HTML QC reports
│   ├── H3K27me3_IP_rep1.html
│   ├── H3K27me3_IP_rep2.html
│   └── Input_rep1.html
└── sample_id.txt
```

## Summary

1. **Understand:** FASTQ files have 4 lines per read; line 4 is the quality score.
2. **Action:** Always run `fastp` to trim adapters, low-quality bases, and too-short reads.
3. **Check:** Use `wc -l` or `awk` for instant feedback on your data size.

> [!NOTE]
> **Up Next:** With clean reads in hand, we're ready to align them to a reference genome using Bowtie2.
