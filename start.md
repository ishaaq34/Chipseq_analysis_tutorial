# Bash Scripting for ChIP-seq Analysis: Sample List Creation

---

## Introduction: Why Learn Bash for Bioinformatics?

**What is Bash?**

Bash is a way to control a computer using text commands instead of clicking through menus. It lets you run programs, manage files, and automate tasks by writing simple commands or scripts.

**Why is this essential for bioinformatics?**

In bioinformatics, you work with large datasets and many specialized tools that must run in a specific sequence. Bash makes these workflows:

- **Repeatable:** Run the same analysis on new data with one command
- **Transparent:** See exactly what was done in plain text
- **Reliable:** Less prone to manual errors from repetitive clicking

Without Bash automation, analyses become harder to track, reproduce, and trust.

---

## The Foundation: Setting Up Safe Scripts

Every reliable Bash script should start with these two critical lines:

```bash
#!/bin/bash
set -euo pipefail
```

**What do these lines do?**

1. **`#!/bin/bash`** – This "shebang" line tells your computer to use the Bash shell. It ensures your script behaves the same way on any system.

2. **`set -euo pipefail`** – This changes Bash's default behavior from "keep going even if something breaks" to "stop immediately when an error occurs."

This single choice separates casual scripting from **defensible scientific analysis**. You want your script to fail loudly if something goes wrong, not silently produce incorrect results.

**Adding practical elements:**

```bash
#!/bin/bash
set -euo pipefail

mkdir -p output
echo "Script started"
```

- **`mkdir -p output`** – Creates a directory called `output` if it doesn't exist. The `-p` flag prevents errors if the directory is already there. This establishes a predictable place for your results.

- **`echo "Script started"`** – Prints a message to confirm the script is running. This is helpful for logging and debugging, especially when scripts run unattended.

---

## The Big Picture: Manual vs. Automated Analysis

Imagine you have a classroom with 100 students and need to process each student's assignment.

**Two approaches:**

1. **Manual Approach:**
   - Type the command for Student 1
   - Wait for it to finish
   - Type for Student 2
   - Wait again...
   - Repeat 98 more times (takes days, prone to typos)

2. **Automated Approach:**
   - Write a small instruction sheet (a script)
   - Give it to a "robot" (your computer)
   - The robot processes all 100 students while you drink coffee

This tutorial teaches you the automated approach for ChIP-seq data.

---

## Part 1: Understanding Sample Lists

### The "Roll Call" Analogy

Imagine grading 100 students. You don't want to type `"Student_John_Doe_Homework_Final_v2.docx"` every single time. You just want a simple list:

- John
- Sarah
- Mike


In bioinformatics, we feel the same way about our sequencing files.

**The Problem:**

Our sequencing files have long, complex names like:

```text
Control_A_H3K9ac_R1.fastq.gz
Control_B_H3K9ac_R1.fastq.gz
```

**The Goal:**

Create a clean list of sample IDs (like a "roll call") that we can feed into automated scripts:

```text
Control_A_H3K9ac
Control_B_H3K9ac

```

**Why do we need this?**

This simple text file will allow our computer to automatically loop through every sample and process them one by one. Instead of typing commands 100 times, we type it once and let the loop do the work.

---

## Part 2: Creating Your Sample List

### Step 1: Verify Your Files

Before creating any list, always check what files you actually have.

```bash
ls *.fastq.gz 
```

**What this does:**
- `ls` = "list" command
- `*.fastq.gz` = show only files ending in `.fastq.gz`

This explicitly targets FASTQ files and is useful at the very start of any sequencing workflow (RNA-seq, ChIP-seq, ATAC-seq, CUT&RUN) to quickly verify that your input files are named correctly and consistently.

---

### Step 2: Extract Clean Sample Names

Now we need to remove the messy file extensions to get clean sample IDs. We'll use a tool called `sed` (Stream Editor) for this.

---

#### Scenario A: Single-End Reads

If you have **single-end sequencing data** , use this approach:

#### Method 1: Step-by-step

This method breaks down the process into individual steps so you can see what's happening at each stage.

> **Note:** In real workflows, your FASTQ files are typically stored in a subdirectory like `raw/` or `fastq_raw/`. It is better to craete the sample list in the `raw/` folder and copy that list to working directory

**Step 1: List all FASTQ files from the raw folder and save to a text file**

```bash
cd raw/
ls *.fastq.gz > samples.txt
```

This command finds all files ending in `.fastq.gz` inside the `raw/` folder and saves their names to `samples.txt`.

**Step 2: Check what got saved**

```bash
cat samples.txt
```

This displays the contents of `samples.txt` so you can verify the filenames. You should see:

```text
Control_A_H3K9ac.fastq.gz
Control_B_H3K9ac.fastq.gz
```

**Step 3: Remove the path and `.fastq.gz` extension from each line**

```bash
sed 's/.fastq.gz//' > sample_id.txt
```

**What this does:**


- Second `sed 's/.fastq.gz//'` = Remove the `.fastq.gz` extension
- `>` = Redirect the final output to a new file
- `sample_id.txt` = Save the cleaned names here

**Step 4: Verify the final result**

```bash
cat sample_id.txt
```

Now you should see clean sample IDs without paths or extensions:

```text
Control_A_H3K9ac
Control_B_H3K9ac
```

---

#### Method 2: One-line command (recommended)

Once you understand the steps above, you can combine them into one efficient command:

```bash
ls *.fastq.gz |  sed 's/.fastq.gz//' > sample_id.txt
```

**What the pipe (`|`) does:**

Instead of saving to intermediate files, the pipe passes output directly from one command to the next:
1. `ls *.fastq.gz` lists the files
3.  `sed` removes the `.fastq.gz` extension  
4. Final result is saved to `sample_id.txt`

This is faster and cleaner than the step-by-step method.  Then moving the sample_id.txt to working directory. 

```
cd ..
mv raw/sample_id.txt .

```



**Check your result:**

```bash
cat sample_id.txt
```

You should see clean sample names:

```text
Control_A_H3K9ac
Control_B_H3K9ac
```

---

#### Scenario B: Paired-End Reads (Most Common)

Most ChIP-seq experiments use **paired-end sequencing**, which produces *two* files per sample:

```text
Control_A_H3K9ac_R1.fastq.gz
Control_A_H3K9ac_R2.fastq.gz
Control_B_H3K9ac_R1.fastq.gz
Control_B_H3K9ac_R2.fastq.gz
```

**The challenge:** We don't want `Control_A_H3K9ac` to appear twice in our list. We only want each sample name *once*.

**The solution:** Target only the `_R1` files.

```bash
ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//' > sample_id.txt
```

**Why look for R1?**

By listing only the `_R1` files, we get exactly one entry per sample. We then strip off the `_R1.fastq.gz` suffix to get the clean sample name.

The script will automatically know to look for both `_R1` and `_R2` files later when we use this list.

---

### Why Not Create the List Manually?

You *could* just type sample names into a text file yourself. However:

- **Prone to typos:** One wrong character breaks your analysis
- **Not verifiable:** Using `ls` ensures you only list files that actually exist

Automation prevents these errors.

---

## Part 3: Using Your Sample List in Automation

Now that we have `sample_id.txt`, we can use it to automate processing of all samples.

### Understanding the Concept

Since we know the sample ID (e.g., `Control_A_H3K9ac`), we can tell our script:

> "For each sample ID, look for two files: the sample name plus `_R1.fastq.gz` and the sample name plus `_R2.fastq.gz`"

The computer can construct these filenames automatically.

---

### Example 1: Basic Loop to Verify File Pairs

This script reads your sample list and prints out the paired filenames:

```bash
#!/bin/bash
set -euo pipefail

while read -r sample; do
  echo "sample_id: $sample"

  fq1="${sample}_R1.fastq.gz"
  fq2="${sample}_R2.fastq.gz"

  echo "paired end: $fq1 : $fq2"
done < sample_id.txt
```

**Output:**

```text
sample_id: Control_A_H3K9ac
paired end:  Control_A_H3K9ac_R1.fastq.gz  :  Control_A_H3K9ac_R2.fastq.gz

sample_id: Control_B_H3K9ac
paired end:    Control_B_H3K9ac_R1.fastq.gz  :  Control_B_H3K9ac_R2.fastq.gz
```

**What this script does:**

- `while read -r sample` – Reads one line (sample ID) at a time from `sample_id.txt`
- `fq1="${sample}_R1.fastq.gz"` – Constructs the forward read filename
- `fq2="${sample}_R2.fastq.gz"` – Constructs the reverse read filename
- `echo` – Prints the results so you can verify

---

### Example 2: Adding Directory Paths

In real workflows, your raw FASTQ files are usually in a specific folder. Let's account for that:

```bash
#!/bin/bash
set -euo pipefail

RAW_DIR="fastq_raw"

while read -r sample; do
  echo "sample_id: $sample"

  fq1="${RAW_DIR}/${sample}_R1.fastq.gz"
  fq2="${RAW_DIR}/${sample}_R2.fastq.gz"

  echo "inputs:"
  echo "  $fq1"
  echo "  $fq2"

done < sample_id.txt
```

**Output:**

```text
sample_id: Control_A_H3K9ac
inputs:
  fastq_raw/Control_A_H3K9ac_R1.fastq.gz
  fastq_raw/Control_A_H3K9ac_R2.fastq.gz

sample_id: Control_B_H3K9ac
inputs:
  fastq_raw/Control_B_H3K9ac_R1.fastq.gz
  fastq_raw/Control_B_H3K9ac_R2.fastq.gz
```

**Key change:**

- `RAW_DIR="fastq_raw"` – Defines where your input files are located
- `${RAW_DIR}/${sample}_R1.fastq.gz` – Constructs the full path to each file

---


### Example 3: Complete Workflow with Input and Output Directories

This example shows a real bioinformatics workflow pattern: reading files from an **input directory** and writing results to an **output directory**.

```bash
#!/bin/bash
set -euo pipefail

RAW_DIR="fastq_raw"
OUT_DIR="bowtie_align"

mkdir -p "$OUT_DIR"

while read -r sample; do
  echo "sample_id: $sample"

  fq1="${RAW_DIR}/${sample}_R1.fastq.gz"
  fq2="${RAW_DIR}/${sample}_R2.fastq.gz"
  bam="${OUT_DIR}/${sample}.sorted.bam"

  echo "inputs:"
  echo "  $fq1"
  echo "  $fq2"

  echo "output:"
  echo "  $bam"

  echo "bowtie2 command: bowtie2 -x hg38_index -1 $fq1 -2 $fq2 | samtools sort -o $bam"
  echo
done < sample_id.txt
```

**Understanding the directory structure:**

This script follows a crucial bioinformatics principle: **separate input data from output results**.

- **Input directory** (`RAW_DIR="fastq_raw"`): Where your raw sequencing files are stored. This should be read-only to preserve original data.
- **Output directory** (`OUT_DIR="bowtie_align"`): Where processed results will be saved. Created automatically if it doesn't exist.

**Breaking down the script:**

1. **`mkdir -p "$OUT_DIR"`** – Creates the output directory before processing starts. The `-p` flag means "don't error if it already exists."

2. **Input file construction:**
   - `fq1="${RAW_DIR}/${sample}_R1.fastq.gz"` – Forward reads from the raw data folder
   - `fq2="${RAW_DIR}/${sample}_R2.fastq.gz"` – Reverse reads from the raw data folder

3. **Output file construction:**
   - `bam="${OUT_DIR}/${sample}.sorted.bam"` – Aligned, sorted BAM file will go in the output folder

4. **The command preview:**
   - Shows what the actual bowtie2 alignment command would be
   - `bowtie2 -x hg38_index` – Uses the human genome reference
   - `-1 $fq1 -2 $fq2` – Paired-end input files
   - `| samtools sort -o $bam` – Pipes output to samtools to create sorted BAM

**Output:**

```text
sample_id: Control_A_H3K9ac
inputs:
  fastq_raw/Control_A_H3K9ac_R1.fastq.gz
  fastq_raw/Control_A_H3K9ac_R2.fastq.gz
output:
  bowtie_align/Control_A_H3K9ac.sorted.bam
bowtie2 command: bowtie2 -x hg38_index -1 fastq_raw/Control_A_H3K9ac_R1.fastq.gz -2 fastq_raw/Control_A_H3K9ac_R2.fastq.gz | samtools sort -o bowtie_align/Control_A_H3K9ac.sorted.bam

sample_id: Control_B_H3K9ac
inputs:
  fastq_raw/Control_B_H3K9ac_R1.fastq.gz
  fastq_raw/Control_B_H3K9ac_R2.fastq.gz
output:
  bowtie_align/Control_B_H3K9ac.sorted.bam
bowtie2 command: bowtie2 -x hg38_index -1 fastq_raw/Control_B_H3K9ac_R1.fastq.gz -2 fastq_raw/Control_B_H3K9ac_R2.fastq.gz | samtools sort -o bowtie_align/Control_B_H3K9ac.sorted.bam
```

**Why organize this way?**

This directory structure keeps your project clean and organized:

```text
your_project/
├── fastq_raw/              ← Original data (never modified)
│   ├── Control_A_H3K9ac_R1.fastq.gz
│   ├── Control_A_H3K9ac_R2.fastq.gz
│   └── ...
├── bowtie_align/           ← Alignment results (can regenerate)
│   ├── Control_A_H3K9ac.sorted.bam
│   └── ...
├── sample_id.txt           ← Sample list
└── analysis_script.sh      ← This script
```

If something goes wrong with alignment, you can safely delete the `bowtie_align/` folder and rerun the script without touching your original raw data.

---




**Best Practice: Separate Folders for Each Process**

In real ChIP-seq analysis, you should create a **separate output folder for each processing step and file type**. This keeps your project organized and makes troubleshooting easier. A typical ChIP-seq project might have folders like:

- `fastq_raw/` - Original FASTQ files from sequencing
- `fastqc_reports/` - Quality control reports
- `trimmed_fastq/` - Adapter-trimmed reads (if needed)
- `bowtie_align/` - Aligned BAM files
- `dedup_bam/` - Deduplicated BAM files
- `peaks/` - Peak calling results from MACS2
- `bigwig/` - Coverage tracks for genome browsers
- `qc_metrics/` - deepTools quality metrics

This folder structure makes it easy to: (1) track which processing stage created which files, (2) quickly find the data you need, (3) delete and regenerate intermediate files without affecting earlier steps, and (4) share your work with collaborators who can immediately understand your project organization.

---





## Summary

**What you learned:**

1. **Safe scripting:** Start every script with `#!/bin/bash` and `set -euo pipefail`
2. **Sample list creation:** Use `ls` and `sed` to extract clean sample IDs from messy filenames
3. **Single-end:** `ls *.fastq.gz | sed 's/.fastq.gz//' > sample_id.txt`
4. **Paired-end:** `ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//' > sample_id.txt`
5. **Automation:** Use `while read` loops to process all samples automatically

**The key insight:**
This simple text file (`sample_id.txt`) acts as a "roll call" or "attendance sheet" for your entire analysis pipeline. It's the foundation that makes large-scale automation possible.

**Next steps:**

In the following tutorials, you'll use this sample list to automate quality control, alignment, peak calling, and visualization for your ChIP-seq experiment.
