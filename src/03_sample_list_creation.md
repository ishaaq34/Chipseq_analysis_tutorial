# Tutorial 02: Creating a Sample List (The "Attendance Sheet")

> [!NOTE]
> **Project Context:** In our specific analysis (see `README.md`), our files are named with ENCODE IDs (e.g., `ENCFF327JFG`). When you run the commands below on our dataset, your `sample_id.txt` will contain these IDs instead of "Sample1" or "Experiment1".

## Level 1: Basic Concept (The "Roll Call")

Imagine you have a classroom with 100 students. If you want to give everyone a grade, you don't want to type out "Student_John_Doe_Homework_Final_v2.docx" every single time. You just want a simple list of names like "John", "Sarah", "Mike".

In bioinformatics, we feel the same way about our files.
*   **The Problem:** Our files have long, messy names like `Sample1_Rep1_R1.fastq.gz`.
*   **The Goal:** We want a clean list of "Sample IDs" (e.g., `Sample1_Rep1`) saved in a text file.
*   **The Why:** We will feed this list to our computer later so it can automatically loop through every sample and process them one by one.

---

## Level 2: Execution (Cleaning the Names)

We will use a command-line tool called `sed` (Stream Editor) to "search and replace" the messy file extensions with nothing, leaving only the clean name.

### Scenario A: Single-End Reads
If your files look like `SampleA.fastq.gz`, we just want to remove `.fastq.gz`.

```
Control_A_H3K9ac.fastq.gz
Control_B_H3K9ac.fastq.gz

```

Run this command:
```bash
ls *.fastq.gz | sed 's/.fastq.gz//' > sample_id.txt

```

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

### Scenario B: Paired-End Reads (Most Common)
Paired-end data usually comes in pairs of files for *one* sample:


```
Control_A_H3K9ac_R1.fastq.gz
Control_A_H3K9ac_R2.fastq.gz
Control_B_H3K9ac_R1.fastq.gz
Control_B_H3K9ac_R2.fastq.gz

```

We don't want "Experiment1" to appear twice in our list. We only want to grab the name *once*.

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
Your sample_id.txt file can now be used in downstream pipeline loops for:

Quality Control (FastQC, MultiQC)

Adapter Trimming (Trim Galore)

Alignment (Bowtie2, HISAT2)

Peak Calling (MACS2)


Since we know the Sample ID (e.g., Input1), we can just tell the robot: "Look for Input1 plus _R1 and Input1 plus _R2."

```bash

for sample in $(cat sample_id.txt); do

  echo "sample_id: $sample"

  fq1="${sample}_R1.fastq.gz"
  fq2="${sample}_R2.fastq.gz"

  echo "paired end:  $fq1  :  $fq2"

done 
```

    sample_id: Control_A_H3K9ac
    paired end:  Control_A_H3K9ac_R1.fastq.gz  :  Control_A_H3K9ac_R2.fastq.gz
    
    sample_id: Control_B_H3K9ac
    paired end:    Control_B_H3K9ac_R1.fastq.gz  :  Control_B_H3K9ac_R2.fastq.gz
   


  
    
 
    



## Summary
1.  **Goal:** Create a clean list of sample names (`sample_id.txt`).
2.  **Tool:** Use `ls` to find files and `sed` to remove extensions.
3.  **Result:** A simple text file that acts as an "Attendance Sheet" for your future automation loops.
