# Tutorial 08: Bam to BigWig (Visualizing the Signal)

## Basic Concept (The Traffic Map)

### Why BigWig?
A **BAM** file gives you the location of every single "car" (read) on the road. It's massive and slow to load.
A **BigWig** file is like the **Google Maps Traffic View**. It doesn't show you the individual cars; it just shows you a green, yellow, or red line indicating "Volume".
*   **Small & Fast:** Compact file size.
*   **Visual:** Perfect for viewing on a Genome Browser (IGV/UCSC).

---

## Execution (The Converter)

The `bamCoverage` tool converts a BAM file into a bigWig file. A bigWig file stores the read coverage across the genome in a compact format that is easy to load in genome browsers like IGV or UCSC. Each option in the command controls how the coverage is calculated.

But first, we need to know the **"Effective Genome Size"**.

### Step 1: Calculate Effective Genome Size
The "Effective" size is the part of the genome that is actually mappable (not repetitive). Different species have different effective genome sizes. Total genome size is simply the full length of all chromosomes, including repeats, low-complexity regions, and stretches of Ns that cannot be mapped. Effective genome size refers only to the portion of the genome where reads can be uniquely aligned. Since the analysis is limited to chr11 and chr12, the effective genome size should be calculated using only these two chromosomes.

There are **two main ways** to estimate effective genome size:

*   **Method 1: faSize (Fast & Simple)**
    The **faSize** tool reports:
    - total bases
    - number of N bases
    - A/T/G/C bases
    - GC content

    The **non-N bases** (total − Ns) are a good approximation of the mappable genome.

    **Step 1.1: Extract chr11 and chr12**
    ```bash
    samtools faidx genome.fasta chr11 chr12 > chr11_chr12.fasta
    ```

    **Step 1.2: Get basic FASTA statistics**
    ```bash
    faSize chr11_chr12.fasta
    ```

    **Step 1.3: Get per-chromosome detailed stats**
    ```bash
    faSize -detailed -tab chr11_chr12.fasta > chr11_chr12.faSize.txt
    ```

    **Step 1.4: Calculate the total effective genome size (non-N bases)**
    ```bash
    awk '{nonN = $2 - $5; sum += nonN} END {print sum}' chr11_chr12.faSize.txt
    ```

    **Example**

    ```
    (chip) rajaishaqnabikhan@Mac % faSize chr11_chr12.fasta

    268361931 bases (690373 N's 267671558 real 267671558 upper 0 lower) in 2 sequences in 1 files
    Total size: mean 134180965.5 sd 1280791.7 min 133275309 (chr12) max 135086622 (chr11) median 135086622
    N count: mean 345186.5 sd 293723.0
    U count: mean 133835779.0 sd 987068.7
    L count: mean 0.0 sd 0.0
    %0.00 masked total, %0.00 masked real
    ```

    ```
    (chip) rajaishaqnabikhan@Mac  % awk '{nonN = $2 - $5; sum += nonN} END {print sum}' chr11_chr12.faSize.txt

    268361931

    ```

*   **Method 2: khmer (Accurate for unique reads)**
    If multimapping reads are removed, the **faSize method becomes less accurate**.
    The **khmer** tool uses the reads themselves to estimate the effective genome size.

    **How it works**

    - Breaks reads into **k-mers**
    - Counts **unique** k-mers
    - Estimates the number of real, mappable bases

    **Typical usage**

    ```bash
    unique-kmers.py -k 21 chr11_chr12.fasta
    ```

    **Example**

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
    **Thus**,

    Raw sequence length (non-N)          = 268,361,931 bp

    Unique 21-mer mappable portion       = 220,798,375 bp

    Difference (unmappable/repetitive)   = ~47.6 million bp

    So about 47.6 Mb of chr11+chr12 is repetitive, low complexity (centromeres, segmental duplications, satellite repeats), or otherwise not uniquely mappable with 21-mers.

    When we calculate effective genome size using faSize, the tool removes N regions but still includes repetitive and low-complexity sequence, because these regions contain real A/T/G/C bases and are mappable as long as multimapping reads are allowed. This makes faSize appropriate when the BAM files still contain multimappers. However, if multimapping reads are removed or MAPQ filtering is applied, the mappable genome becomes smaller because repeats and low-complexity regions no longer contribute usable signal. In that case, the khmer tool—especially unique-kmers.py—gives a better estimate because it counts only unique k-mers from the reads and naturally excludes repeats, low-complexity regions, and N stretches. For uniquely mapped data, the khmer-based estimate is the correct effective genome size to use.

    In our BAM file, multimapping reads have already been excluded and low-quality alignments removed by MAPQ filtering. Because these steps bias us toward uniquely mappable regions, we estimate the effective genome size using a k-mer–based approach (e.g. khmer).

### Step 2: Run bamCoverage
Now we make the BigWig.
```bash
bamCoverage \
  -b sample.bam \
  -o sample.bw \
  --normalizeUsing RPGC \
  --effectiveGenomeSize 220798375 \
  --binSize 10 \
  --smoothLength 30 \
  --numberOfProcessors 4 \
  --ignoreDuplicates
```

---

## Fine Tuning 

### 3.1 Bin Size (Resolution)
*   **Concept:** Think of this as the "Resolution" of your image.
A bin is simply a small section of the genome. If the bin size is 100 base pairs, the chromosome is cut into pieces that are each 100 bases long. If the bin size is 1000 base pairs, each piece is much larger. When we make bigWig files, we don’t measure the signal at every single base. Instead, we take each bin and calculate the average signal inside that bin. This makes the data easier to compare and much faster to compute.
*   **10 bp bin:** High definition (HD). You see every detail, but the file is bigger. Small bins give very detailed tracks but create millions of bins, which slows down the analysis.
*   **100 bp bin:** Standard definition. Good for overview, faster to process. Larger bins, like 200 bp or 1000 bp, give less detail but are much faster and still good for things like correlation, PCA, and general quality checks. In most cases, larger bins (100–200 bp for a few chromosomes, around 1000 bp for whole-genome work) give results that are almost the same but cleaner and quicker to process.

### 3.2 Smoothing (Blurring the Photo)
*   **Concept:** Smoothing averages the signal of nearby bins.
Smoothing reduces noise in the coverage track by averaging nearby bins.
*   **Why use it?** Raw data can be "pixelated" (jumpy). Smoothing applies a slight blur to make the biological peaks stand out against the background noise.
Example: If bin size is 20 and smooth length is 60, then each point is the average of a bin plus its neighbors left and right.
*   **Rule:** If smoothLength > binSize, it averages the neighbors. If the smoothing window is smaller than the bin size, it is ignored.

### 3.3 Normalization (The Currency Exchange)
Samples have different sequencing depths (Total Reads). One sample might have 10 million reads, another 20 million.
If we don’t adjust for that, a sample with more reads will always look “stronger,” even if the biology is the same. Normalization fixes this. It makes all samples easier to compare by putting them on the same scale.
If we compare them directly, the 20M sample will look twice as strong. This is unfair.

**Normalization Options:**
You can choose different ways to normalize your data, depending on what you want to compare.
1.  **CPM (Counts Per Million):** The classic "Per Capita" correction. Corrects only for the total number of reads in your sample, which is enough when all bins are the same size.
2.  **RPKM (Reads Per Kilobase per Million):** Corrects for how many reads you have and how big each bin is, so big bins don’t look stronger just because they are longer.
3.  **RPGC (Reads Per Genome Coverage):** *Preferred for ChIP-seq.*
    *   **Logic:** "What would this signal look like if we had exactly 1x coverage of the genome?"
    *   **Why:** It creates a standardized "Currency" (1x coverage) that makes biological sense. RPGC changes the data to show what it would look like if you had exactly 1× genome coverage, which is why it needs the effective genome size.

If you pick “None,” you see the raw, uncorrected coverage.

---

## Summary
1.  **Format:** Convert BAM (Raw Lists) to BigWig (Signal Tracks).
2.  **Calculation:** Use **Effective Genome Size** to handle repeats correctly.
3.  **Normalization:** Use **RPGC** to make fair comparisons between samples.
