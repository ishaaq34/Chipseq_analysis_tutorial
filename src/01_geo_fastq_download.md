#  02: Getting the Raw Data (FASTQ)

> [!NOTE]
> **Project Context:** While this tutorial teaches you how to download data from SRA, our specific project utilizes **ENCODE** data (CEBPA, H3K27me3, H3K9ac). The `fastq-dl` tool is excellent for public SRA data, but ENCODE data is often downloaded directly from the ENCODE portal. I have worked on [SRA data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115704#:~:text=GEO%20Accession%20viewer&text=GEO%20help:%20Mouse%20over%20screen%20elements%20for%20information.&text=Here%20we%20report%20that%20the,a%20sperm%2Dspecific%20chromatin%20signature ) to demonstrate steps from data download to the alignemnt step.



Complete information for all sequencing runs associated with this repository is available through the [NCBI SRA Run Selector (PRJNA475794)](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA475794&o=acc_s%3Aa). The Run Selector provides an interactive interface to inspect sequencing metadata, including library layout, platform, read length, and experimental design.

From this interface, you can download the full metadata table as well as a plain accession list containing the SRR identifiers. The accession list can be saved as `srr_list.txt` and used directly for automated data retrieval. Unwanted runs can be removed from this file before download, allowing precise control over which datasets are processed.



## Level 1: Basic Concept 

### The "Library" Analogy
Before we download anything, it helps to understand where the data lives. Think of the public databases like a University Library.
*   **GEO (Gene Expression Omnibus):** This is the **Card Catalog**. It has the descriptions (Metadata) of the experiments—like "H3K27ac in Breast Cancer"—but it usually doesn't hold the actual heavy books (data files).
*   **SRA (Sequence Read Archive) & ENA (European Nucleotide Archive):** These are the **Stacks**. This is where the actual raw data files are stored.

**Your Goal:** You find an interesting study in the Catalog (GEO), get its ID, and then send a runner (our software tool) to the Stacks (SRA/ENA) to fetch the files.

### The Files You Will Encounter
As we process data, we change formats. Think of it like this:
1.  **FASTQ (Raw Reads):** The jumbled words. This is what comes off the sequencer. It’s just millions of short strings of letters (A, C, G, T) with no context.
2.  **BAM (Aligned Reads):** The assembled sentences. We map the words to a reference genome so we know where they belong.
3.  **BED/BigWig (Signals):** The highlighted passages. These are simplified files that show us where the interesting "peaks" (protein binding sites) are.

---

## Level 2: Fetching the data

To download data, we use a tool called [fastq-dl](https://github.com/rpetit3/fastq-dl). It acts like a smart librarian—you just give it the ID number, and it deals with the complicated databases for you.

### 2.1 Download a Single Sample
If you have a Run ID (starts with **SRR** or **ERR**), use this command:

```bash
# Download one sample from SRA
fastq-dl --accession SRR7297994 --provider SRA --cpus 4

# Or from ENA (often faster/more reliable)
fastq-dl --accession SRR7297994 --provider ena
```

### 2.2 Download Multiple Samples (The Loop)
Usually, you need to download many samples. Instead of typing the command 10 times, we put the IDs in a list.

1.  Create a file named `srr_list.txt` with one ID per line:
    ```text
    SRR7297994
    SRR7297995
    SRR7297998
    SRR7298003
    ```


2.  Run this "Loop" to download them one by one:

```
#!/bin/bash
set -euo pipefail

RAW_DIR="fastq_raw"
mkdir -p "$RAW_DIR"

while read -r acc; do
  echo "Downloading accession: $acc"

  fastq-dl \
    --accession "$acc" \
    --provider SRA \
    --cpus 1 \
    --outdir "$RAW_DIR"

  echo "Finished downloading: $acc"
done < srr_list.txt

```



 

### 2.3 Parallel Download (The Fast Way)
If you have a powerful computer, you can download multiple files at the same time using `parallel`.


```
#!/bin/bash
set -euo pipefail

mkdir -p fastq_raw

parallel -j 4 \
  'fastq-dl --accession {} --provider SRA --cpus 1 --outdir fastq_raw' \
  :::: srr_list.txt

```
 Adding echos to remain updated about what is going on!!
```
#!/bin/bash
set -euo pipefail

mkdir -p fastq_raw

parallel -j 4 \
  'echo "Starting download: {}" &&
   fastq-dl --accession {} --provider SRA --cpus 1 --outdir fastq_raw &&
   echo "Finished download: {}"' \
  :::: srr_list.txt
```

### 2.4 Download an Entire Study
You can also download everything associated with a study ID (starts with **SRP** or **PRJNA**):

```bash
fastq-dl --accession SRP115709
```
*Note: Be careful! A whole study might have hundreds of files.*

---


### Connecting GEO to SRA
How do we find the **SRR** numbers?
In GEO, you will see a hierarchy. It's important not to mix these up:

| Level | Prefix (SRA / ENA) | What It Is |
| :--- | :--- | :--- |
| **Project** | PRJNA / PRJEB | The umbrella project (e.g., "Breast Cancer Epigenomics 2024"). |
| **Study** | SRP / ERP | A specific paper or dataset. |
| **Sample** | SRS / ERS | The biological sample (e.g., "Patient 5 Tumor"). |
| **Experiment** | SRX / ERX | The library prep info. |
| **Run** | **SRR / ERR** | **The actual data.** This is what you download. |

### Technical Replicates (Multi-lane)
Sometimes, one biological sample is sequenced across multiple "lanes" of a machine to get more reads.
*   **Result:** You might see `SRR900100` and `SRR900101` for the *same* sample.
*   **Action:** `fastq-dl` will download them. Later, we will merge these FASTQ files together so we have one big file for that sample.

---


## After Downloading FASTQ Files

```text
chipseq_tutorial/
├── raw/                    ← Raw FASTQ files from sequencing
│   ├── Control_A_H3K9ac_R1.fastq.gz
│   ├── Control_A_H3K9ac_R2.fastq.gz
│   ├── Control_B_H3K9ac_R1.fastq.gz
│   ├── Control_B_H3K9ac_R2.fastq.gz
│   └── ...
└── sample_id.txt           ← Sample list for automation
```

**Key Points:**
- Original sequencing data is never modified
- Each sample has paired-end reads (_R1 and _R2)
- `sample_id.txt` contains clean sample names for automation

---







## Summary
1.  **Understand:** GEO is for metadata; SRA/ENA is for data.
2.  **Identify:** Find the **SRR** (Run) IDs for your samples.
3.  **Download:** Use `fastq-dl` with a loop or parallel command to fetch the FASTQ files.
