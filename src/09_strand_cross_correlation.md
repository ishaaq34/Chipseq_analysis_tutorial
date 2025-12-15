# Tutorial 07: Strand Cross-Correlation (The Echo)

##  1: Basic Concept (The Echo)

### Why do we see two peaks?
When you do ChIP-seq, the DNA fragments are 3D objects, about 200bp long, with the protein in the middle.
However, the sequencer only reads the **Ends** of these fragments (5' ends).
*   **Forward Strand Reads:** Read from the left end (Start of fragment).
*   **Reverse Strand Reads:** Read from the right end (End of fragment).

This creates two piles of reads separated by the fragment length, like two mountains with a valley in between.

### The Echo Analogy
Imagine you shout **"HELLO"** (Forward Reads). A split second later, you hear the echo **"HELLO"** (Reverse Reads).
*   **Cross-Correlation** is measuring exactly how long that delay is.
*   We slide the Forward reads towards the Reverse reads. When they overlap perfectly, the "volume" is loudest (Max Correlation).
*   The distance we slid them telling us the **True Fragment Length**.

---
<img width="639" height="403" alt="Screenshot 2025-12-10 at 11 39 06 AM" src="https://github.com/user-attachments/assets/cfc18aa5-2c0d-4772-ad2a-f4ae8c0f578d" />


*(Source:[Genome-wide analysis of transcription factor binding sites based on ChIP-Seq data](https://www.nature.com/articles/nmeth.1246)*

---

## 2: PhantomPeakQualTools



We use a tool called `run_spp.R` (part of [PhantomPeakQualTools]( https://github.com/kundajelab/phantompeakqualtools) to calculate this. It finds the "Best Match" distance.

### The Command
```bash
# Run PhantomPeakQualTools
Rscript /opt/anaconda3/envs/chip/bin/run_spp.R \
      -c=sample.bam \
      -savp=Sample1_spp.qc.pdf \
      -s=0:5:400 \
      -out=Sample1_spp.qc.txt
```

**Explanation:**
*   `-c`: Input BAM file.
*   `-savp`: Saves the diagnostic PDF plot (The "Echo" graph).
*   `-out`: Output file containing the score numbers (NSC and RSC).
*   `-s=0:5:400`: Scan shift distances from 0 to 400 bp in steps of 5.

---

## Level 3: Analysis (Signal vs Noise)

The output plot (`Sample1_spp.qc.pdf`) usually shows TWO peaks. This is where quality control happens.

### 3.1 The Peaks
1.  **The Real Peak (Red Line):**
    *   **What is it?** The “Echo”. The point where Forward and Reverse reads overlap because they bind the same protein.
    *   **Location:** Usually around 150-250 bp (your fragment size).
    *   **Meaning:** Represents **Biological Signal**.

2.  **The Phantom Peak (Blue Line):**
    *   **What is it?** “Microphone Feedback”. It occurs at the **Read Length** (e.g., 50bp or 100bp).
    *   **Why?** It's caused by mapping artifacts and "sticky" sequences. It happens in *every* experiment, even bad ones.
    *   **Meaning:** Represents **Background Noise**.


  `Sample1_spp.qc.txt` provides numeric QC metrics in colunmns from 1-11

```
H3K9ac_ENCFF193NPE.chr11_12.bam   #1  BAM file name (input to SPP QC)
3411159                          #2  Number of reads used
220                              #3  Estimated fragment length (bp)
0.545418                         #4  Cross-correlation at fragment length
55                               #5  Phantom peak shift (bp)
0.501850                         #6  Cross-correlation at phantom peak
1500                             #7  Maximum strand shift tested (bp)
0.373681                         #8  Minimum cross-correlation (background)
1.459581                         #9  NSC – Normalized Strand Cross-correlation
1.339926                         #10 RSC – Relative Strand Cross-correlation
1                               #11 SPP quality tag (1 = good)
```


  
---

<img width="758" height="504" alt="Screenshot 2025-12-15 at 9 23 05 PM" src="https://github.com/user-attachments/assets/696f703d-a182-466f-aeb3-bf19a07b0272" />


---

### 3.2 The Metrics (NSC & RSC)
We compare the Height of the Real Peak (Signal) to the Phantom Peak (Noise).

**Key Metrics Table:**

| Metric | Full Name | Meaning | Good Threshold |
| :--- | :--- | :--- | :--- |
| **NSC** | Normalized Strand Cross-correlation | **Signal-to-Noise Ratio.** How much higher is the Real Peak (COL4) than the flat background (COL8)? (COL4 / COL8) | **> 1.05** |
| **RSC** | Relative Strand Cross-correlation | **Signal vs Phantom.** Is the Real Peak (COL4) clearly stronger than the Phantom Peak (COL6)? ((COL4 – COL8) / (COL6 – COL8)) | **> 0.8** |

### 3.3 Interpreting the Data (Example Analysis)

---

<img width="943" height="470" alt="image3" src="https://github.com/user-attachments/assets/5b03189d-8fb7-4f4c-8145-5b5930921e05" />



 ---
*   **Inputs (Control):**
    *   NSC values sit almost exactly at background **(1.003 and 1.005)**. RSC values are only **0.62 and 0.64** (Low).
    *   **Interpretation:** Just noise. This is normal for inputs! The numerically high correlation values (0.5509 and 0.5474) are driven by the dominant phantom peak, not real enrichment.
*   **H3K9ac (Active Mark):**
    *   Correlation is **High (0.5454 and 0.4112)**.
    *   **Interpretation:** Strong signal. Acetylation marks usually give huge peaks, reflecting the broad and high-coverage nature of these marks.
*   **H3K27me3 (Repressive Mark):**
    *   Correlation is **Medium (0.3708 and 0.3572)**.
    *   **Interpretation:** Expected. Repressive marks are broad and diffuse, so the "Echo" is quieter, fitting a repressive mark that produces wide but moderate enrichment.
*   **CEBPA (Transcription Factor):**
    *   Correlation is **Lower (0.2876 and 0.1979)**.
    *   **Interpretation:** TF peaks are sharp but rare (small % of genome). Total signal is lower, but the peaks are distinct. This is expected because TF peaks are sharp and occupy a small fraction of the genome, so their genome-wide cross-correlation values are naturally smaller.

---

## Summary
1.  **Cross-Correlation** shifts reads to find the fragment length (The Echo).
2.  **Phantom Peak** is a background artifact at read length (Microphone Feedback).
3.  **RSC > 0.8** means your Signal (Real Peak) is louder than your Noise (Phantom Peak).
