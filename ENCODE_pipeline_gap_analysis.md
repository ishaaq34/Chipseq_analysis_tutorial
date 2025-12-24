# ENCODE ChIP-seq Pipeline - Gap Analysis

**Generated:** 2025-12-22  
**Purpose:** Identify parameters, features, and analytical approaches in the ENCODE ChIP-seq pipeline (github.com/ENCODE-DCC/chip-seq-pipeline2) that are NOT covered in the current tutorial

---

## Executive Summary

This document compares the [ENCODE ChIP-seq Pipeline](https://github.com/ENCODE-DCC/chip-seq-pipeline2) with the current ChIP-seq tutorial to identify missing parameters, quality control metrics, and methodological approaches. The analysis is organized by pipeline stage.

### Key Gaps Identified

1. **Pipeline Infrastructure** (WDL workflow, cloud platforms, containerization)
2. **Advanced Alignment Parameters** (BWA support, custom aligners, trimming options)
3. **Subsampling & Control Depth Management**
4. **BAM Redaction** (Privacy-focused masking)
5. **Pseudoreplication Analysis** (True replicates + pseudo-replicates for optimal IDR)
6. **Advanced Cross-Correlation Parameters** (Exclusion ranges, PE-specific settings)
7. **Peak Calling Flexibility** (SPP support, cap_num_peak, customizable thresholds)
8. **Comprehensive QC Reporting** (Automated HTML/JSON reports)
9. **Blacklist Filtering** (Genomic region exclusion)
10. **Chromosomal Filtering & Regex Patterns**

---

## 1. Pipeline Architecture & Infrastructure

### What ENCODE Has (Missing in Tutorial)

#### 1.1 Workflow Definition Language (WDL)

- **ENCODE Implementation:**
  - Entire pipeline defined in WDL (`chip.wdl`)
  - Task-based parallelization
  - Automatic resource allocation
  - Cloud-native execution (Terra, DNAnexus, AWS, Google Cloud)
  
- **Tutorial Status:** ❌ **Not covered**
  - Uses manual bash scripts and loops
  - Sequential execution (no automatic parallelization)
  - Local execution only

**Missing Tutorial Content:**

```
- WDL workflow syntax
- `caper` (Cromwell wrapper for WDL execution)
- Cloud platform integration
- HPC schedulers (SLURM, SGE, PBS) integration
- Automatic retry mechanisms (maxRetries)
```

#### 1.2 Containerization

- **ENCODE Implementation:**
  - Docker images for all tools
  - Singularity support for HPC
  - Version-locked dependencies
  
- **Tutorial Status:** ❌ **Not covered**
  - Relies on conda environments
  - No containerization discussion

#### 1.3 Genome Database System

- **ENCODE Implementation:**
  - Pre-built genome TSV files for hg38, hg19, mm10, mm9
  - Automated genome database downloader/builder
  - Cloud-hosted genome databases (Google Storage, DNAnexus)
  - Second blacklist support (`chip.blacklist2`)
  
- **Tutorial Status:** ⚠️  **Partially covered**
  - Shows manual genome download
  - No discussion of genome TSV database system
  - No second blacklist (`blacklist2`) parameter

---

## 2. Input Data & Preprocessing

### What ENCODE Has (Missing in Tutorial)

#### 2.1 Flexible Pipeline Entry Points

**ENCODE allows starting from:**

- FASTQ (`chip.fastqs_repX_R1`, `chip.fastqs_repX_R2`)
- BAM (`chip.bams`)
- Filtered/nodup BAM (`chip.nodup_bams`)
- TAG-ALIGN (`chip.tas`)
- **Can mix data types per replicate** (e.g., FASTQ for rep1, BAM for rep2)

**Tutorial Status:** ⚠️  **Partially covered**

- Shows FASTQ → BAM workflow
- Does not discuss TAG-ALIGN format
- No mixing of data types per replicate

#### 2.2 FASTQ Preprocessing Parameters

| Parameter | ENCODE Default | Tutorial Coverage |
|-----------|----------------|-------------------|
| `chip.crop_length` | 0 (disabled) | ❌ Not mentioned |
| `chip.crop_length_tol` | 2 | ❌ Not mentioned |
| `chip.trimmomatic_phred_score_format` | auto | ❌ Not mentioned |

**Purpose of Missing Parameters:**

- `crop_length`: Trim all reads to fixed length (useful for variable-length runs)
- `crop_length_tol`: Minimum read length after cropping (removes very short reads)
- `trimmomatic_phred_score_format`: Override auto-detection of Phred33/64 encoding

**Tutorial Uses:** `fastp` for QC/trimming (different tool, similar purpose)

---

## 3. Alignment Stage

### What ENCODE Has (Missing in Tutorial)

#### 3.1 Multiple Aligner Support

| Aligner | ENCODE Support | Tutorial Coverage |
|---------|----------------|-------------------|
| **Bowtie2** | ✅ Default | ✅ **Covered** |
| **BWA** | ✅ Supported | ❌ Not mentioned |
| **Custom aligner** | ✅ Via Python script | ❌ Not mentioned |

**ENCODE BWA-specific parameters:**

```json
{
  "chip.aligner": "bwa",
  "chip.use_bwa_mem_for_pe": false,
  "chip.bwa_mem_read_len_limit": 70
}
```

**Logic:** For PE data with read length < 70bp → use `bwa aln`, otherwise use `bwa mem`

**Tutorial Status:** ❌ **Not covered** - Only demonstrates Bowtie2

#### 3.2 Bowtie2 Advanced Parameters

| Parameter | ENCODE Support | Tutorial Coverage |
|-----------|----------------|-------------------|
| `chip.use_bowtie2_local_mode` | `--local` flag option | ❌ Not mentioned |

**Purpose:** Local mode allows soft-clipping of read ends (useful for adapters, poor quality ends)

**Tutorial Uses:** Default end-to-end mode only

#### 3.3 Custom Aligner Framework

**ENCODE allows:**

- Custom Python script for alignment (`chip.custom_align_py`)
- Custom aligner index (`chip.custom_aligner_idx_tar`)

**Example use case:** Novel aligner (e.g., STAR for RNA-ChIP, minimap2 for long reads)

**Tutorial Status:** ❌ **Not covered**

---

## 4. BAM Filtering & Quality Control

### What ENCODE Has (Missing in Tutorial)

#### 4.1 Duplicate Marker Choice

| Tool | ENCODE Support | Tutorial Coverage |
|------|----------------|-------------------|
| Picard | ✅ Default | ✅ **Covered** |
| Sambamba | ✅ Alternative | ❌ Not mentioned |

**ENCODE Parameter:**

```json
{
  "chip.dup_marker": "picard"  // or "sambamba"
}
```

**Why Sambamba?** Faster than Picard, useful for very large datasets. Use when Picard fails.

**Tutorial Uses:** Picard + samtools (two methods), no Sambamba

#### 4.2 Skip Duplicate Removal Option

**ENCODE Parameter:**

```json
{
  "chip.no_dup_removal": false
}
```

**Purpose:** For ultra-low input ChIP (single-cell, CUT&RUN), removing duplicates may discard real signal.

**Tutorial Status:** ❌ **Not covered** - Assumes duplicates are always removed

#### 4.3 BAM Redaction

**ENCODE Parameter:**

```json
{
  "chip.redact_nodup_bam": false
}
```

**What it does:**

- Masks positional information in BAM files using `ptools`
- Converts CIGAR strings to simple `M` only (e.g., `76M5D20M` → `101M`)
- Preserves alignment for analysis but removes identifying information

**Purpose:**

- Data sharing compliance (remove patient-identifiable genomic positions)
- Privacy for human data

**Side effects:**

- `bedtools bamtobed` read length calculation becomes approximate
- Slight noise in TAG-ALIGN conversion

**Tutorial Status:** ❌ **Not covered** - No discussion of privacy/data sharing

---

## 5. Subsampling & Read Depth Control

### What ENCODE Has (Missing in Tutorial)

#### 5.1 Manual Subsampling

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `chip.subsample_reads` | 0 (disabled) | Downsample IP reads to fixed depth |
| `chip.ctl_subsample_reads` | 0 (disabled) | Downsample control reads to fixed depth |
| `chip.xcor_subsample_reads` | 15000000 | Subsample reads for cross-correlation ONLY |

**Important distinction:**

- `subsample_reads`: Affects ALL downstream analyses (peak calling, BigWig, etc.)
- `xcor_subsample_reads`: Only for cross-correlation QC (saves compute time)

**ENCODE cross-correlation optimization:**

```json
{
  "chip.xcor_subsample_reads": 15000000  // 15M reads for xcor
}
```

**Why?** Cross-correlation is computationally expensive. 15M reads provide stable estimates without full dataset processing.

**Tutorial Status:** ❌ **Not covered** - No subsampling discussion

#### 5.2 Automatic Control Depth Management

**ENCODE Parameters:**

```json
{
  "chip.ctl_depth_limit": 200000000,           // 200M read hard limit
  "chip.exp_ctl_depth_ratio_limit": 5.0        // Max IP:Control ratio
}
```

**What it does:**

- If control has > 200M reads → subsample to 200M
- If control has > 5× IP depth → subsample control to 5× IP depth

**Purpose:**

- Prevents memory overflow from huge control files
- Avoids over-correction (too-deep controls can mask real signal)

**Tutorial Status:** ❌ **Not covered** - No automatic depth balancing

---

## 6. Cross-Correlation Analysis (PhantomPeakQualTools)

### What ENCODE Has (Missing in Tutorial)

#### 6.1 Advanced Cross-Correlation Parameters

| Parameter | ENCODE Default | Tutorial Coverage |
|-----------|----------------|-------------------|
| `chip.xcor_trim_bp` | 50 | ❌ Not mentioned |
| `chip.use_filt_pe_ta_for_xcor` | false | ❌ Not mentioned |
| `chip.xcor_exclusion_range_min` | -500 | ❌ Not mentioned |
| `chip.xcor_exclusion_range_max` | Auto | ❌ Not mentioned |

**Parameter Details:**

**1. `xcor_trim_bp` (50 bp)**

- Trim R1 FASTQ to 50bp for cross-correlation analysis ONLY
- Why? Faster computation, removes low-quality 3' ends
- Does NOT affect alignment or peak calling

**2. `use_filt_pe_ta_for_xcor`**

- If `true`: Use filtered PE BAM/TAG-ALIGN for cross-correlation (ignores trimming)
- If `false`: Use trimmed R1 FASTQ

**3. `xcor_exclusion_range_min` & `_max`**

- Corresponds to phantompeakqualtools `-x=<min>:<max>` parameter
- Excludes specific shift ranges from correlation calculation
- Default: `-500` to `max(read_len + 10, 50)` for TF, `max(read_len + 10, 100)` for histone

**Purpose:** Avoid spurious peaks near read length

**Tutorial Status:** ⚠️  **Partially covered**

- Shows basic `run_spp.R` usage
- Does NOT explain ENCODE-specific trimming or exclusion ranges

---

## 7. Control Sample Management

### What ENCODE Has (Missing in Tutorial)

#### 7.1 Intelligent Control Selection

**ENCODE Parameters:**

```json
{
  "chip.always_use_pooled_ctl": true,
  "chip.ctl_depth_ratio": 1.2
}
```

**Control selection logic:**

**Case 1: Single control provided**
→ Use it for all replicates

**Case 2: Multiple controls provided**

- If `always_use_pooled_ctl = true` → Pool all controls, use for all replicates
- If `always_use_pooled_ctl = false`:
  - If control depth ratio > 1.2 (one control 20% deeper than others) → Pool controls
  - Otherwise → Pair rep1_IP with rep1_Ctl, rep2_IP with rep2_Ctl

**Purpose:** Maximize statistical power while avoiding batch effect confounding

**Tutorial Status:** ⚠️  **Partially covered**

- Shows 1:1 IP-Control pairing
- No discussion of pooling strategies or automatic selection

---

## 8. Peak Calling

### What ENCODE Has (Missing in Tutorial)

#### 8.1 Peak Caller Choice

| Caller | ENCODE Support | Tutorial Coverage |
|--------|----------------|-------------------|
| **MACS2/MACS3** | ✅ Default for histone | ✅ **Covered (MACS3)** |
| **SPP** | ✅ Default for TF | ❌ Not mentioned |

**ENCODE Parameter:**

```json
{
  "chip.peak_caller": "spp"  // or "macs2"
}
```

**SPP (run_spp.R) for TF ChIP:**

- More conservative than MACS2 for sharp peaks
- Requires control
- Better for low-signal TFs

**Tutorial Uses:** MACS3 for both TF and histone (which is acceptable)

#### 8.2 Peak Capping

**ENCODE Parameters:**

```json
{
  "chip.cap_num_peak_macs2": 500000,
  "chip.cap_num_peak_spp": 300000
}
```

**What it does:** Filters peaks, keeping only top `N` by signal value

**Purpose:**

- Prevents memory overflow in IDR with millions of low-confidence peaks
- ENCODE standard for reproducibility

**Tutorial Status:** ❌ **Not covered** - No peak capping

#### 8.3 FDR Threshold for SPP

**ENCODE Parameter:**

```json
{
  "chip.fdr_thresh": 0.01  // For SPP only
}
```

**Parameter passed to:** `run_spp.R -fdr=0.01`

**Tutorial Uses:** No SPP, so not applicable

#### 8.4 IDR Threshold

**ENCODE Parameter:**

```json
{
  "chip.idr_thresh": 0.05
}
```

**Tutorial Coverage:** ✅ **Covered** (`IDR < 0.05` filtering shown)

---

## 9. Pseudoreplication & Optimal IDR

### What ENCODE Has (Missing in Tutorial)

#### 9.1 True Replicate Only Mode

**ENCODE Parameter:**

```json
{
  "chip.true_rep_only": false
}
```

**If `false` (default):** Pipeline runs FULL ENCODE IDR workflow:

1. **True replicates:** Rep1 vs Rep2
2. **Self-consistency (pseudoreplicates):**
   - Pool Rep1+Rep2 → Shuffle → Split into Pseudo-rep1 and Pseudo-rep2
   - Run IDR on Pseudo-rep1 vs Pseudo-rep2
3. **Compare:** True IDR vs Pseudo IDR
   - If True > Pseudo → Good reproducibility
   - If Pseudo > True → Poor biological replicates

**If `true`:** Skip pseudoreplication (faster, but less rigorous QC)

**Tutorial Status:** ❌ **Not covered**

- Only shows True replicate IDR
- No pseudoreplication or self-consistency checks

#### 9.2 Pseudoreplication Random Seed

**ENCODE Parameter:**

```json
{
  "chip.pseudoreplication_random_seed": 0
}
```

**What it does:**

- If `0`: Use TAG-ALIGN file size (bytes) as random seed (deterministic)
- If `> 0`: Use specified integer as seed

**Purpose:** Reproducible pseudoreplication (same input → same pseudo-reps)

**Tutorial Status:** ❌ **Not covered**

---

## 10. Advanced QC Metrics

### What ENCODE Has (Missing in Tutorial)

#### 10.1 Deeptools JSD (Jensen-Shannon Distance)

**ENCODE Parameter:**

```json
{
  "chip.enable_jsd": true
}
```

**What it measures:**

- `plotFingerprint` from deepTools
- Quantifies enrichment vs input using Jensen-Shannon divergence
- Complements FRiP

**Tutorial Status:** ⚠️  **Partially covered**

- Shows `plotFingerprint` in section 10
- Does NOT explain JSD metric or ENCODE standard interpretation

#### 10.2 GC Bias Calculation

**ENCODE Parameter:**

```json
{
  "chip.enable_gc_bias": true
}
```

**What it does:**

- `deepTools computeGCBias`
- Detects PCR amplification bias toward GC-rich or AT-rich regions

**When to check:**

- High duplicate rates
- Uneven coverage in promoters vs gene bodies

**Tutorial Status:** ❌ **Not covered** - No GC bias discussion

#### 10.3 Count Signal Track

**ENCODE Parameter:**

```json
{
  "chip.enable_count_signal_track": false
}
```

**What it generates:**

- BigWig with raw read counts (not normalized)
- Useful for very low-depth samples where normalization introduces noise

**Tutorial Status:** ❌ **Not covered** - Only shows normalized BigWig

---

## 11. Blacklist & Chromosomal Filtering

### What ENCODE Has (Missing in Tutorial)

#### 11.1 Dual Blacklist Support

**ENCODE Parameters:**

```json
{
  "chip.blacklist": "path/to/ENCODE_blacklist.bed",
  "chip.blacklist2": "path/to/custom_blacklist.bed"
}
```

**What it does:**

- Merges two blacklist files
- Filters out peaks overlapping blacklisted regions
- Common use: ENCODE blacklist + custom problematic regions

**Tutorial Status:** ❌ **Not covered** - No blacklist filtering

#### 11.2 Chromosomal Filtering

**ENCODE Parameters:**

```json
{
  "chip.filter_chrs": ["chrM", "chrY"],
  "chip.regex_bfilt_peak_chr_name": "chr[\\dXY]+"
}
```

**What it does:**

- `filter_chrs`: Remove reads from specified chromosomes during BAM filtering
- `regex_bfilt_peak_chr_name`: Keep peaks only on chromosomes matching regex

**Default regex:** `chr[\\dXY]+` keeps chr1-22, chrX, chrY (removes chrM, alt contigs)

**Common use cases:**

- Remove mitochondrial contamination (`chrM`)
- Female samples: remove `chrY`
- Remove unplaced/alt contigs

**Tutorial Status:** ❌ **Not covered** - No discussion of chromosomal filtering

---

## 12. Fragment Length Handling

### What ENCODE Has (Missing in Tutorial)

#### 12.1 Manual Fragment Length Override

**ENCODE Parameter:**

```json
{
  "chip.fraglen": [150, 180, 200]  // Array: one per replicate
}
```

**When to use:**

- Cross-correlation fails (negative fragment length)
- Very low signal data
- TF with very few peaks

**How it works:** Overrides `run_spp.R` estimated fragment lengths

**Tutorial Status:** ⚠️  **Partially covered**

- Shows cross-correlation analysis
- Does NOT show manual override for failed cases

---

## 13. HTML QC Report Generation

### What ENCODE Has (Missing in Tutorial)

#### 13.1 Automated QC Reports

**ENCODE generates:**

- `qc.html`: Interactive HTML with all metrics, plots, and pass/fail indicators
- `qc.json`: Machine-readable QC metrics for database submission

**Includes:**

- Alignment stats (per replicate)
- Cross-correlation plots (NSC, RSC, fragment length)
- FRiP scores
- IDR plots and peak counts
- Library complexity
- GC bias plots
- Fingerprint plots (JSD)

**Tutorial Uses:** MultiQC for aggregation (different tool, similar result)

**Gap:** ENCODE's QC report is ChIP-seq specific with ENCODE thresholds pre-annotated

---

## 14. Pipeline Modes

### What ENCODE Has (Missing in Tutorial)

#### 14.1 Pipeline Type

**ENCODE Parameter:**

```json
{
  "chip.pipeline_type": "tf"  // or "histone" or "control"
}
```

**Modes:**

| Mode | Purpose | Behavior |
|------|---------|----------|
| `tf` | Transcription factor | Default narrow peak parameters, SPP caller, requires control |
| `histone` | Histone modification | Broad peak support, MACS2 default, can work without control |
| `control` | Map controls only | No peak calling, only alignment/filtering |

**Control mode details:**

- Sets `align_only = true` automatically
- Skips cross-correlation, JSD, GC-bias
- Used for pre-processing controls for later pairing

**Tutorial Status:** ⚠️  **Partially covered**

- Shows TF and histone analysis separately
- No `control`-only mode

#### 14.2 Align-Only Mode

**ENCODE Parameter:**

```json
{
  "chip.align_only": false
}
```

**If `true`:**

- Stops after BAM filtering (no peak calling)
- Generates BAMs and TAG-ALIGNs only

**Use case:** Pre-process data, then call peaks with different parameters later

**Tutorial Status:** ❌ **Not covered**

---

## 15. Output Organization

### What ENCODE Has (Missing in Tutorial)

#### 15.1 Automatic Output Organization

**ENCODE provides:**

- `croo` tool (Cromwell Output Organizer)
- Organizes scattered WDL outputs into structured directories
- Generates file manifest

**Tutorial Uses:** Manual directory structure (user creates folders)

#### 15.2 QC Metrics Spreadsheet

**ENCODE provides:**

- `qc2tsv` tool
- Converts `qc.json` → TSV spreadsheet for Excel/comparison

**Tutorial Status:** ❌ **Not covered**

---

## 16. Summary of Critical Gaps

### High-Priority Missing Features

| Feature | ENCODE Implementation | Tutorial Status | Impact |
|---------|----------------------|-----------------|---------|
| **Pseudoreplication** | Standard workflow | ❌ Missing | High - ENCODE requirement |
| **Blacklist filtering** | Dual blacklist support | ❌ Missing | High - Reduces false positives |
| **Peak capping** | `cap_num_peak` | ❌ Missing | Medium - Prevents IDR overflow |
| **Chromosomal filtering** | `filter_chrs` + regex | ❌ Missing | Medium - Removes contaminants |
| **Control depth management** | Auto-subsampling | ❌ Missing | Medium - Prevents over-correction |
| **SPP peak caller** | TF default | ❌ Missing | Low - MACS3 acceptable alternative |
| **GC bias QC** | `enable_gc_bias` | ❌ Missing | Low - Good for troubleshooting |
| **BWA aligner** | Alternative to Bowtie2 | ❌ Missing | Low - Bowtie2 sufficient |

### Medium-Priority Missing Parameters

- `crop_length` / `crop_length_tol` (FASTQ trimming)
- `ctl_depth_ratio` (control selection logic)
- `xcor_exclusion_range_min/max` (cross-correlation tuning)
- `no_dup_removal` (low-input ChIP option)
- `fraglen` manual override (failed cross-correlation backup)

### Low-Priority (Infrastructure)

- WDL workflow syntax
- Cloud platform execution (Terra, DNAnexus)
- Docker/Singularity containers
- `caper` tool usage

---

## 17. Recommendations for Tutorial Enhancement

### Quick Wins (Can Add Immediately)

1. **Add Blacklist Filtering Section:**

   ```bash
   # Download ENCODE blacklist
   wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
   
   # Filter peaks
   bedtools intersect -v -a peaks.narrowPeak -b hg38-blacklist.v2.bed > peaks.blacklist_filtered.bed
   ```

2. **Add Chromosomal Filtering:**

   ```bash
   # Remove chrM reads during alignment
   samtools view -h sample.bam | \
     awk '$1 ~ /^@/ || $3 != "chrM"' | \
     samtools view -b > sample.no_chrM.bam
   ```

3. **Add Peak Capping:**

   ```bash
   # Keep top 300,000 peaks by signal
   sort -k7,7nr peaks.narrowPeak | head -300000 > peaks.capped.narrowPeak
   ```

### Medium Effort (Requires New Sections)

4. **Pseudoreplication Analysis:**
   - Create section explaining ENCODE self-consistency framework
   - Show how to shuffle and split TAG-ALIGN
   - Run IDR on pseudo-replicates

5. **Control Depth Management:**
   - Add section on checking control vs IP depth
   - Show subsampling commands when control >> IP

6. **GC Bias QC:**
   - Add `deepTools computeGCBias` to QC section
   - Explain interpretation

### Optional (Advanced)

7. **SPP Peak Caller:**
   - Add alternative peak calling with `run_spp.R`
   - Compare SPP vs MACS3 for TFs

8. **Fragment Length Override:**
   - Add troubleshooting section for failed cross-correlation
   - Show manual `--extsize` parameter for MACS3

---

## 18. Parameter-by-Parameter Comparison Table

### Complete ENCODE Input JSON Parameters

#### Pipeline Type & Mode

| Parameter | Default | Tutorial Covered? |
|-----------|---------|-------------------|
| `chip.pipeline_type` | `tf` | ⚠️  Partial (no `control` mode) |
| `chip.align_only` | `false` | ❌ No |
| `chip.true_rep_only` | `false` | ❌ No (pseudorep missing) |
| `chip.redact_nodup_bam` | `false` | ❌ No |

#### Alignment

| Parameter | Default | Tutorial Covered? |
|-----------|---------|-------------------|
| `chip.aligner` | `bowtie2` | ✅ Yes (Bowtie2 only) |
| `chip.crop_length` | `0` | ❌ No |
| `chip.crop_length_tol` | `2` | ❌ No |
| `chip.trimmomatic_phred_score_format` | `auto` | ❌ No (uses fastp) |
| `chip.use_bwa_mem_for_pe` | `false` | ❌ No (BWA not covered) |
| `chip.bwa_mem_read_len_limit` | `70` | ❌ No |
| `chip.use_bowtie2_local_mode` | `false` | ❌ No |

#### Filtering

| Parameter | Default | Tutorial Covered? |
|-----------|---------|-------------------|
| `chip.mapq_thresh` | `30` | ✅ Yes |
| `chip.dup_marker` | `picard` | ✅ Yes (+ samtools) |
| `chip.no_dup_removal` | `false` | ❌ No |
| `chip.filter_chrs` | `[]` | ❌ No |

#### Subsampling

| Parameter | Default | Tutorial Covered? |
|-----------|---------|-------------------|
| `chip.subsample_reads` | `0` | ❌ No |
| `chip.ctl_subsample_reads` | `0` | ❌ No |
| `chip.xcor_subsample_reads` | `15000000` | ❌ No |
| `chip.ctl_depth_limit` | `200000000` | ❌ No |
| `chip.exp_ctl_depth_ratio_limit` | `5.0` | ❌ No |

#### Cross-Correlation

| Parameter | Default | Tutorial Covered? |
|-----------|---------|-------------------|
| `chip.xcor_trim_bp` | `50` | ❌ No |
| `chip.use_filt_pe_ta_for_xcor` | `false` | ❌ No |
| `chip.xcor_exclusion_range_min` | `-500` | ❌ No |
| `chip.xcor_exclusion_range_max` | Auto | ❌ No |

#### Controls

| Parameter | Default | Tutorial Covered? |
|-----------|---------|-------------------|
| `chip.always_use_pooled_ctl` | `true` | ⚠️  Partial (1:1 pairing shown) |
| `chip.ctl_depth_ratio` | `1.2` | ❌ No |

#### Peak Calling

| Parameter | Default | Tutorial Covered? |
|-----------|---------|-------------------|
| `chip.peak_caller` | `spp` (TF) / `macs2` (histone) | ⚠️  Partial (MACS3 only) |
| `chip.cap_num_peak_macs2` | `500000` | ❌ No |
| `chip.cap_num_peak_spp` | `300000` | ❌ No |
| `chip.pval_thresh` | `0.01` | ⚠️  Partial (uses `-q`) |
| `chip.idr_thresh` | `0.05` | ✅ Yes |
| `chip.fdr_thresh` | `0.01` | ❌ No (SPP-specific) |

#### Blacklist & Regex

| Parameter | Default | Tutorial Covered? |
|-----------|---------|-------------------|
| `chip.blacklist` | - | ❌ No |
| `chip.blacklist2` | - | ❌ No |
| `chip.regex_bfilt_peak_chr_name` | `chr[\\dXY]+` | ❌ No |

#### Advanced QC

| Parameter | Default | Tutorial Covered? |
|-----------|---------|-------------------|
| `chip.enable_jsd` | `true` | ⚠️  Partial (plotFingerprint shown, but JSD not explained) |
| `chip.enable_gc_bias` | `true` | ❌ No |
| `chip.enable_count_signal_track` | `false` | ❌ No |

#### Fragment Length

| Parameter | Default | Tutorial Covered? |
|-----------|---------|-------------------|
| `chip.fraglen` | Auto (from xcor) | ⚠️  Partial (no manual override) |

#### Pseudoreplication

| Parameter | Default | Tutorial Covered? |
|-----------|---------|-------------------|
| `chip.pseudoreplication_random_seed` | `0` | ❌ No |

---

## 19. References

### ENCODE Pipeline Documentation

- Main repository: <https://github.com/ENCODE-DCC/chip-seq-pipeline2>
- Input parameters: <https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/docs/input.md>
- ENCODE ChIP-seq standards: <https://www.encodeproject.org/chip-seq/transcription_factor/>
- Google Doc specification: <https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c>

### Tools Referenced

- Caper (Cromwell wrapper): <https://github.com/ENCODE-DCC/caper>
- PhantomPeakQualTools: <https://github.com/kundajelab/phantompeakqualtools>
- IDR: <https://github.com/nboley/idr>
- ptools (redaction): <https://github.com/ENCODE-DCC/ptools>
- MACS3: <https://github.com/macs3-project/MACS>
- deepTools: <https://deeptools.readthedocs.io/>

---

## 20. Conclusion

The current tutorial provides excellent coverage of the **core ChIP-seq workflow** (alignment → deduplication → QC → peak calling → IDR → motif analysis). However, it is missing several **ENCODE-specific rigor components**:

### Critical Gaps for ENCODE Compliance

1. **Pseudoreplication** (self-consistency checks)
2. **Blacklist filtering** (standard in ENCODE)
3. **Peak capping** (prevents computational issues)
4. **Chromosomal filtering** (removes chrM contamination)

### Recommended Next Steps

1. Add blacklist filtering section (easiest, high impact)
2. Add pseudoreplication workflow (most important for ENCODE standards)
3. Add chromosomal filtering examples
4. Add troubleshooting section for failed cross-correlation (manual fraglen override)

The tutorial is already comprehensive for **academic research** and **publication-ready analysis**. The missing ENCODE features are mostly for **consortium-level standardization** and **edge-case handling**.

---

**End of Gap Analysis**
