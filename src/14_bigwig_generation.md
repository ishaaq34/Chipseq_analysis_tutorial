# Bam to BigWig (Visualizing the Signal)

`bamCoverage` `bigWig` `normalization` `RPGC` `effective-genome-size` `visualization`

## 1. Basic Concept (The Traffic Map)

### Why BigWig?

---

A **BAM** file gives you the location of every single "car" (read) on the road. It's massive and slow to load.

A **BigWig** file is like the **Google Maps Traffic View**. It doesn't show you the individual cars; it just shows you a green, yellow, or red line indicating "Volume".

* **Small & Fast:** Compact file size.
* **Visual:** Perfect for viewing on a Genome Browser (IGV/UCSC).

---

## 2. Requirements (Effective Genome Size)

The "Effective" size is the part of the genome that is mappable.

* **Unique Regions:** Easy to map.
* **Repetitive Regions:** Hard to map, but reads can technically align there (Ambiguous).

**Formula 1: Total Raw Length**
$$
\text{Total Genome Length} = \underbrace{\text{Effective Genome}}_{\text{All Mappable}} + \underbrace{\text{N-bases}}_{\text{Gaps/Unsequenced}}
$$

**Formula 2: Breaking Down "Effective Genome"**
$$
\text{Effective Genome} = \underbrace{\text{Unique Genome}}_{\text{Easy}} + \underbrace{\text{Repetitive DNA}}_{\text{Ambiguous}}
$$

### How we calculated it

There are **two main ways** to estimate this:

#### Method 1: faSize (Calculates General Effective Genome)

The **faSize** tool counts all non-N bases. This assumes that if a read aligns to a repeat, it still counts as "mapped".

**Step 1.1: Extract chr11 and chr12**

```bash
samtools faidx genome.fasta chr11 chr12 > chr11_chr12.fasta
```

**Step 1.2: Get statistics**

```bash
faSize chr11_chr12.fasta
```

**Output:**

```
268361931 bases (690373 N's 267671558 real 267671558 upper 0 lower) in 2 sequences inside 1 file
```

**Step 1.3: Calculate total effective size (non-N)**

```bash
awk '{nonN = $2 - $5; sum += nonN} END {print sum}' chr11_chr11_chr12.faSize.txt
# Result: 268361931
```

**Result:** `268,361,931 bp` (Unique + Repetitive)

#### Method 2: khmer (Calculates Unique Effective Genome)

The **khmer** tool counts only unique k-mers. It strictly ignores repetitive zones.

**Command:**

```bash
unique-kmers.py -k 21 chr11_chr12.fasta
```

**Output:**

```
(chip) rajaishaqnabikhan@Mac human % unique-kmers.py -k 21 chr11_chr12.fasta

|| This is the script unique-kmers.py in khmer.
|| You are running khmer version 3.0.0a3
|| You are also using screed version 1.1.3
||
|| If you use this script in a publication, please cite EACH of the following:
||
||   * MR Crusoe et al., 2015. https://doi.org/10.12688/f1000research.6924.1
||   * A. Döring et al. https://doi.org:80/10.1186/1471-2105-9-11
||   * Irber and Brown. https://doi.org/10.1101/056846
||
|| Please see http://khmer.readthedocs.io/en/latest/citations.html for details.

Estimated number of unique 21-mers in chr11_chr12.fasta: 220798375
Total estimated number of unique 21-mers: 220798375
```

### The Result

* **Raw length (Method 1: faSize):** `268,361,931 bp` (Unique + Repetitive)
* **Unique length (Method 2: khmer):** `220,798,375 bp` (Unique Only)
* **The Difference:** `47,563,556 bp` of repetitive/unmappable DNA (Ambiguous).

### The Decision

Which number do we use?

* Therefore, we must use the **khmer** result (`Unique Only`), because our data does not contain reads mapped to repetitive regions.

> [!TIP]
> **Use this number:** `220798375`

---

## 3. Execution (The Converter)

We use `bamCoverage` to convert BAM files into BigWig tracks. We will normalize using **RPGC** (Reads Per Genome Coverage) so all tracks are comparable (1x coverage scale).

### Step 3.1: Create Output Directory

Keep your workspace clean.

```bash
mkdir -p bigwigs
```

### Step 3.2: Run bamCoverage Loop

We process all BAM files using our `sample_id.txt` list.

```bash
# Loop through each sample ID in the text file
cat sample_id.txt | while read id; do
  
  echo "Generating BigWig for: $id"
  
  bamCoverage \
    -b encode_bam/${id}.bam \                    
    -o bigwigs/${id}.bw \                        
    --binSize 10 \                               
    --normalizeUsing RPGC \                      
    --effectiveGenomeSize 2701495761 \            
    --smoothLength 30 \                          
    --numberOfProcessors 4 

done
```

**What this does:**

1. **Reads** the BAM file.
2. **Chops** the genome into 10bp bins (buckets).
3. **Counts** reads in each bin.
4. **Normalizes** the count so 1.0 = 1x genome coverage.
5. **Smooths** the signal (averaging neighbors) to make the peaks look cleaner.
6. **Saves** the result as a `.bw` file.

---

## 4. Fine Tuning (Under the Hood)

### 4.1 Bin Size (Resolution)

* **Concept:** Think of this as the **Resolution** of your image.
* **Small Bin (10bp):** "HD" resolution. You see every tiny peak, but the file is larger and slower to compute.
* **Large Bin (100bp):** "SD" resolution. Good for zooming out and looking at broad trends. Faster to process.
* **Decision:** For transcription factors (sharp peaks), 10bp is great. For broad histone marks (H3K27me3), 50-100bp is often sufficient. We use **10bp** here for high detail.

### 4.2 Smoothing (Blurring the Photo)

* **Concept:** Smoothing averages the signal of nearby bins to reduce "jumpy" noise.
* **Why?** Raw sequencing data can be pixelated. Smoothing applies a slight blur to make the biological signal (the hill) stand out against the background noise (the grass).
* **Our Setting:** `--smoothLength 30`. This averages the signal over a 30bp window.

### 4.3 Normalization (The Currency Exchange)

Samples have different sequencing depths.

* **Sample A:** 40 million reads (Rich)
* **Sample B:** 20 million reads (Poor)

If we don't fix this, Sample A will look huge just because it has more money (reads). **Normalization** puts everyone on the same currency.

* **RPGC (Reads Per Genome Coverage):** *Preferred for ChIP-seq.*
  * **Logic:** "What would this signal look like if we had exactly **1x coverage** of the genome?"
  * **Why:** It creates a standardized "Currency" (1x coverage) that makes biological sense. A value of "5.0" in the track means "5 times more reads than random background".

---

## Summary

1. **Input:** BAM files (`encode_bam/`).
2. **Action:** `bamCoverage` with **RPGC** normalization.
3. **Output:** BigWig files (`bigwigs/`) ready for IGV visualization.

> [!NOTE]
> **Up Next:** Now that we have our signals (BigWigs) and our QC (Fingerprints) done, we are ready to call peaks! (Wait, technically we usually call peaks *before* or *parallel* to visualization, but viewing tracks helps confirm peak calls).

---

## Directory Structure After BigWig Generation

```text
chipseq_tutorial/
├── bam_files_final/            ← Source BAM files
├── macs3_results/              ← MACS3 peak files
├── idr/                        ← IDR outputs
└── bigwigs/                    ← **NEW: Normalized BigWig tracks**
    ├── H3K9ac_ENCFF534IPX.bw
    ├── H3K9ac_ENCFF193NPE.bw
    ├── H3K27me3_ENCFF164ALR.bw
    ├── H3K27me3_ENCFF532DQH.bw
    ├── Input_ENCFF110SOB.bw
    └── Input_ENCFF919XCV.bw
```

---
