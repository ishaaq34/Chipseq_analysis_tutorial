

Bash is a way to control a computer using text commands instead of clicking through menus. It lets you run programs, manage files, and automate tasks by writing simple commands or scripts.

In bioinformatics, this is essential because analyses involve large datasets and many tools that must run in a specific order. Bash makes these workflows repeatable, transparent, and less prone to manual errors. Without it, analyses are harder to track, reproduce, and trust.

```
#!/bin/bash
set -euo pipefail

```

This is the defensible starting point. The first line fixes the shell so the script behaves the same way everywhere. The second line changes Bash’s default behavior from “keep going no matter what” to “stop when something is wrong.” That single choice separates casual scripting from analysis you can stand behind.

```
#!/bin/bash
set -euo pipefail

mkdir -p output
echo "Script started"
```

mkdir -p output establishes a predictable place for results. That is practical, not fundamental. 

echo "Script started" is simply a trace. It helps when the script runs unattended or produces logs, but it has no effect on correctness.



Imagine you have a classroom with 100 students. You have two options:
1.  **Manual:** Type the command for Sample 1. Wait. Type for Sample 2. Wait... (This takes days and is prone to typos).
2.  **Automated:** Write a small instruction sheet (Script) and give it to a robot. The robot does all 100 samples while you drink coffee.



## Level 1: Basic Concept (The "Roll Call")

Imagine you have a classroom with 100 students. If you want to give everyone a grade, you don't want to type out `"Student_John_Doe_Homework_Final_v2.docx"` every single time. You just want a simple list of names like `"John"`, `"Sarah" `, `"Mike"`.

In bioinformatics, we feel the same way about our files.
*   **The Problem:** Our files have long, messy names like `Sample1_Rep1_R1.fastq.gz`.
*   **The Goal:** We want a clean list of "Sample IDs" (e.g., `Sample1_Rep1`) saved in a text file.
*   **The Why:** We will feed this list to our computer later so it can automatically loop through every sample and process them one by one.

---

## Level 2: Execution (Cleaning the Names)

We will use a command-line tool called `sed` (Stream Editor) to "search and replace" the messy file extensions with nothing, leaving only the clean name.



Start by checking the files in your directory:

```bash
ls *.fastq.gz 
```

This explicitly targets FASTQ files.Useful at the very start of any workflow (RNA-seq, ChIP-seq, ATAC-seq, CUT&RUN) to quickly verify that your input files are named correctly and consistently.

---

# 2. Create a list of sample IDs (single-end)

If you have single-end reads, remove the `.fastq.gz` suffix using **sed**:

```bash
ls *.fastq.gz > samples.txt          # lists all .fastq.gz files and saves the list into samples.txt
cat samples.txt                      # prints the contents of samples.txt so you can check the filenames
sed 's/.fastq.gz//' sample.txt> sample_id.txt   # tries to remove '.fastq.gz' from each line of sample.txt 
                                                 # s/ = substitute & // means “replace it with nothing”
                                                

```
or more simple way
```
ls *.fastq.gz | sed 's/.fastq.gz//' > sample_id.txt  ## lists all .fastq.gz files, pipes them into sed 

```

The file `sample_id.txt` will contain one sample ID per line.



```
cat sample_id.txt
```
will give 

```
Control_A_H3K9ac
Control_B_H3K9ac
```


**What did this do?**
1.  `ls *.fastq.gz`: Listed all the files.
2.  `|`: Passed that list to the next tool.
3.  `sed 's/.fastq.gz//'`: Replaced the text ".fastq.gz" with nothing (effectively deleting it).
4.  `> sample_id.txt`: Saved the result to a new file.

likwise , 
Usually, you need to download many samples. Instead of typing the command 10 times, we put the IDs in a list.

1.  Create a file named `srr_list.txt` with one ID per line:
    ```text
    SRR7297994
    SRR7297995
    SRR7297998
    SRR7298003
    ```

2.  Run this "Loop" to download them one by one:
   

### Scenario B: Paired-End Reads (Most Common)
Paired-end data usually comes in pairs of files for *one* sample:


```
Control_A_H3K9ac_R1.fastq.gz
Control_A_H3K9ac_R2.fastq.gz
Control_B_H3K9ac_R1.fastq.gz
Control_B_H3K9ac_R2.fastq.gz

```

We don't want `"Experiment1" ` to appear twice in our list. We only want to grab the name *once*.

Run this command:
```bash
ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//' > sample_id.txt
```

**Why look for R1?**
By listing only the `_R1` files, we get exactly one entry per sample. We then strip off the `_R1.fastq.gz` part to get the clean sample name.

---


### Creating the List Manually?
You *could* just type the names into a text file yourself. However, that is prone to typos. Using `ls` ensures you only list files that actually exist.


## Ready for Use
thus  `sample_id.txt` file can now be used in downstream pipeline loops for:


* Quality Control (FastQC, MultiQC)

* Adapter Trimming (Trim Galore)

* Alignment (Bowtie2, HISAT2)




Since we know the `Sample ID (e.g., Control_A_H3K9ac)`, we can just tell the robot: `"Look for Control_A_H3K9ac plus _R1 ` and `Control_A_H3K9ac plus _R2 `


```
#!/bin/bash
set -euo pipefail

while read -r sample; do
  echo "sample_id: $sample"

  fq1="${sample}_R1.fastq.gz"
  fq2="${sample}_R2.fastq.gz"

  echo "paired end: $fq1 : $fq2"
done < sample_id.txt
```

  
    sample_id: Control_A_H3K9ac
    paired end:  Control_A_H3K9ac_R1.fastq.gz  :  Control_A_H3K9ac_R2.fastq.gz
    
    sample_id: Control_B_H3K9ac
    paired end:    Control_B_H3K9ac_R1.fastq.gz  :  Control_B_H3K9ac_R2.fastq.gz
   


Now if we have input fastq in the fastq_raw folder. 

```

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

sample_id: Control_A_H3K9ac
inputs:
  fastq_raw/Control_A_H3K9ac_R1.fastq.gz
  fastq_raw/Control_A_H3K9ac_R2.fastq.gz

sample_id: Control_B_H3K9ac
inputs:
  fastq_raw/Control_B_H3K9ac_R1.fastq.gz
  fastq_raw/Control_B_H3K9ac_R2.fastq.gz

    
 

## Summary
1.  **Goal:** Create a clean list of sample names (`sample_id.txt`).
2.  **Tool:** Use `ls` to find files and `sed` to remove extensions.
3.  **Result:** A simple text file that acts as an "Attendance Sheet" for your future automation loops.
