# Tutorial 03: Understanding and Cleaning Your Data (FASTQ & fastp)

## Level 1: Basic Concept (The Anatomy of a Read)

A FASTQ file is just a text file full of DNA sequences. But unlike a simple list of letters, every single read carries extra baggage (its quality score).

Think of every read like a **Luggage Tag** with 4 lines of information:
1.  **Line 1 (The Header):** Starts with `@`. This is the **ID Card**. It tells you the machine name, flowcell lane, and coordinates.
2.  **Line 2 (The Sequence):** The DNA letters (`ACTG...`). This is the **Content** inside the bag.
3.  **Line 3 (The Spacer):** Starts with `+`. Just a divider.
4.  **Line 4 (The Quality):** A string of weird characters (`F:F#,,...`). This is the **Trust Score**. Each character represents the probability that the corresponding base in Line 2 is wrong.

**Example Read:**
```text
@SRR7298010.1 1/1                 <-- ID: Read #1
NGATTTCTATTCTTGGAACCATTAAAA...    <-- DNA: "N" means the machine failed to call that base
+SRR7298010.1 1/1                 <-- Spacer
FFFFF:F#,,...                     <-- Quality: Matches the length of the DNA
```

---

## Level 2: Execution (The Car Wash)

Before we start analyzing, we need to clean our data.
*   **The Problem:** Sequencers sometimes make mistakes, especially at the ends of reads. They also leave "adapters" (artificial tags) attached to the DNA.
*   **The Solution:** We use a tool called **fastp**. It acts like an automatic car wash: dirty reads go in, clean reads come out.

### 2.1 Basic Cleaning (Single-End)
```bash
# -i: Input (dirty)
# -o: Output (clean)
fastp -i raw_sample.fastq.gz -o clean_sample.fastq.gz
```

### 2.2 Basic Cleaning (Paired-End)
For paired-end data, we must process R1 and R2 together so they stay synchronized.

```bash
fastp -i in.R1.fq.gz -I in.R2.fq.gz \
      -o out.R1.fq.gz -O out.R2.fq.gz \
      -f 3 -t 2
```
*   `-f 3`: Trim 3 bases from the front (optional, if you know the start is bad).
*   `-t 2`: Trim 2 bases from the tail.

**What does fastp do automatically?**
*   **Quality Filtering:** Drops reads if too many bases have low scores (default: Phred < 15).
*   **Length Filtering:** Drops reads that become too short after trimming (default: < 15bp).
*   **Adapter Removal:** Finds and cuts off adapter sequences automatically.

---

## Level 3: Advanced Analysis (The Math)

### 3.1 Quick Stats with AWK
Sometimes you don't want to run a full report; you just want to know "How many reads do I have?"
You can use `awk` (a math tool for text) to count directly from the compressed file.

**Count Total Reads:**
```bash
# A FASTQ record is 4 lines. We count total lines and divide by 4.
gzcat sample.fastq.gz | wc -l | awk '{print $1/4 " reads"}'
```

**Count Total Bases (Coverage):**
```bash
# Sums the length of line 2 (sequence) for every record
gzcat sample.fastq.gz | awk 'NR%4==2 {b+=length($0)} END{print b " bases"}'
```
*   *Approximation:* If you have 100 Million bases and your genome is 3 Billion bases (Human), your coverage is roughly 0.03x.

### 3.2 Batch Processing
If you have 50 files, you can use a script to run `fastp` on all of them in parallel.
The `fastp` developers provide a handy script called `parallel.py`:

```bash
# Process 3 files at a time (-f 3), using 2 threads per file (-t 2)
python parallel.py -i /raw_data_folder -o /clean_data_folder -r /report_folder -f 3 -t 2
```
This automatically finds pairs and generates HTML reports for every sample.

---

## Summary
1.  **Understand:** FASTQ files have 4 lines per read; line 4 is the quality score.
2.  **Action:** Always run `fastp` to trim adapters, low-quality bases, and too-short reads.
3.  **Check:** Use `wc -l` or `awk` for instant feedback on your data size.
