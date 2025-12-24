# nf-core/chipseq - MACS2/MACS3 Parameters

**nf-core ChIP-seq pipeline version:** 2.1.0  
**Peak Caller:** MACS3 (upgraded from MACS2 in recent versions)  
**Generated:** 2025-12-22

---

## Overview

The nf-core/chipseq pipeline uses **MACS3** (formerly MACS2) for peak calling with configurable parameters. Here are the exact parameters and commands used.

---

## 1. Core MACS Parameters in nf-core

### **Peak Calling Mode**

```bash
# Narrow peak mode (default for TFs)
--narrow_peak

# Broad peak mode (for histone marks)
# Omit --narrow_peak and set:
--broad_cutoff 0.1
```

---

### **Significance Thresholds**

| Parameter | Default | Description | MACS Flag |
|-----------|---------|-------------|-----------|
| `--macs_fdr` | `0.05` | Minimum FDR (q-value) cutoff | `-q` |
| `--macs_pvalue` | Not set | P-value cutoff (mutually exclusive with FDR) | `-p` |

**Important:** `--macs_fdr` and `--macs_pvalue` are mutually exclusive. If you use `--macs_pvalue`, q-values will NOT be calculated.

---

### **Genome Size**

```bash
# Using shortcut
--macs_gsize 'hs'   # Human
--macs_gsize 'mm'   # Mouse
--macs_gsize 'ce'   # C. elegans
--macs_gsize 'dm'   # Drosophila

# Using numerical value
--macs_gsize 2.7e9  # Human genome size
```

**MACS flag:** `-g`

---

### **Fragment Size Handling**

By default, MACS3 **auto-estimates** fragment size. However, nf-core allows manual override:

```bash
# Let MACS auto-estimate (default)
# No parameters needed

# Manual override (if auto-detection fails)
# In nf-core config or modules, you can specify:
# --nomodel --extsize <fragment_size>
```

---

## 2. Exact MACS3 Commands Used by nf-core

### **For Narrow Peaks (TFs)**

```bash
macs3 callpeak \
    -t <IP.bam> \
    -c <Input.bam> \
    -f BAM \
    -g hs \
    -n <sample_name> \
    -q 0.05 \
    --keep-dup all \
    --call-summits \
    --outdir peaks/
```

**Key Parameters:**

- `-t`: Treatment/IP BAM file
- `-c`: Control/Input BAM file
- `-f BAM`: Input format
- `-g hs`: Genome size (human)
- `-n`: Name prefix
- `-q 0.05`: FDR cutoff (default)
- `--keep-dup all`: Keep all duplicate reads (nf-core filters duplicates earlier in the pipeline)
- `--call-summits`: Find precise summit positions

---

### **For Broad Peaks (Histone Marks)**

```bash
macs3 callpeak \
    -t <IP.bam> \
    -c <Input.bam> \
    -f BAM \
    -g hs \
    -n <sample_name> \
    -q 0.05 \
    --keep-dup all \
    --broad \
    --broad-cutoff 0.1 \
    --outdir peaks/
```

**Additional Broad Peak Parameters:**

- `--broad`: Enable broad peak calling
- `--broad-cutoff 0.1`: Cutoff for linking nearby enriched regions

---

### **Optional: Create bedGraph Files**

```bash
# If --save_macs_pileup is enabled in nf-core
macs3 callpeak \
    ... (above parameters) ... \
    --bdg \
    --SPMR
```

**Added flags:**

- `--bdg`: Generate bedGraph output
- `--SPMR`: Signal per million reads normalization

---

## 3. Complete nf-core MACS Parameters Table

| nf-core Parameter | Default | MACS3 Flag | Description |
|-------------------|---------|------------|-------------|
| `--narrow_peak` | `false` | (no `--broad`) | Call narrow peaks for TFs |
| `--broad_cutoff` | `0.1` | `--broad-cutoff` | Broad peak cutoff value |
| `--macs_fdr` | `0.05` | `-q` | Minimum FDR (q-value) cutoff |
| `--macs_pvalue` | `null` | `-p` | P-value cutoff (overrides FDR) |
| `--macs_gsize` | Auto-detected | `-g` | Effective genome size |
| `--save_macs_pileup` | `false` | `--bdg --SPMR` | Save normalized bedGraph tracks |
| `--min_reps_consensus` | `1` | N/A (post-processing) | Min replicates for consensus peaks |
| `--skip_peak_qc` | `false` | N/A | Skip MACS peak QC plots |
| `--skip_peak_annotation` | `false` | N/A | Skip HOMER annotation |
| `--skip_consensus_peaks` | `false` | N/A | Skip consensus peak generation |

---

## 4. nf-core Pre-Processing Before MACS

**Important:** nf-core performs extensive filtering BEFORE peak calling:

### BAM Filtering Steps

1. ✅ Remove blacklisted regions
2. ✅ Remove duplicates (via Picard)
3. ✅ Remove non-primary alignments
4. ✅ Remove unmapped reads
5. ✅ Remove multi-mappers
6. ✅ Remove reads with >4 mismatches (BAMTools)
7. ✅ Remove insert size >2kb (paired-end)
8. ✅ Remove reads mapping to different chromosomes (paired-end)
9. ✅ Remove reads not in FR orientation (paired-end)

**Result:** MACS receives **highly filtered, clean BAM files**

This is why nf-core uses `--keep-dup all` in MACS - duplicates already removed upstream!

---

## 5. Comparison: nf-core vs Standard MACS Usage

### **nf-core Approach:**

```bash
# Step 1: Filter BAM extensively (via nf-core modules)
picard MarkDuplicates ...
samtools view -F 1804 ...
bamtools filter -tag "nM:<4" ...

# Step 2: Simple MACS command on clean data
macs3 callpeak -t clean.bam -c clean_input.bam -f BAM -g hs -q 0.05 --keep-dup all
```

### **Standard Tutorial Approach:**

```bash
# Step 1: Basic duplicate removal
picard MarkDuplicates REMOVE_DUPLICATES=true ...

# Step 2: MACS handles filtering
macs3 callpeak -t dedup.bam -c dedup_input.bam -f BAM -g hs -q 0.01 --keep-dup 1
```

**Key Difference:**

- **nf-core:** Extensive pre-filtering → Simple MACS command
- **Tutorial:** Basic pre-filtering → MACS does more filtering

---

## 6. Fragment Length Estimation

### **nf-core Behavior:**

1. **Auto-estimation (default):**
   - MACS3 builds shift model from data
   - Works for most high-quality datasets

2. **Manual override (if needed):**

   ```bash
   # In nf-core config or custom module:
   --nomodel --extsize 200
   ```

3. **Use phantompeakqualtools estimate:**
   - nf-core runs `run_spp.R` for QC
   - Fragment length from SPP can inform manual override if MACS fails

---

## 7. Output Files Generated

### **Narrow Peak Mode:**

```
sample_peaks.narrowPeak    # BED6+4 format with summit
sample_peaks.xls           # Detailed peak info
sample_summits.bed         # Summit positions only
sample_model.r             # R script for shift model plot
sample_treat_pileup.bdg    # Treatment pileup (if --save_macs_pileup)
sample_control_lambda.bdg  # Control lambda (if --save_macs_pileup)
```

### **Broad Peak Mode:**

```
sample_peaks.broadPeak     # BED6+3 format (no summits)
sample_peaks.gappedPeak    # BED12+3 format with gaps
sample_peaks.xls
sample_model.r
```

---

## 8. Example nf-core Command

### **Running the Pipeline:**

```bash
nextflow run nf-core/chipseq \
  --input samplesheet.csv \
  --outdir results \
  --genome GRCh38 \
  --narrow_peak \
  --macs_fdr 0.01 \
  --macs_gsize 2.7e9 \
  --save_macs_pileup \
  -profile docker
```

**This translates to MACS3 command:**

```bash
macs3 callpeak \
  -t filtered_IP.bam \
  -c filtered_Input.bam \
  -f BAM \
  -g 2.7e9 \
  -n sample \
  -q 0.01 \
  --keep-dup all \
  --call-summits \
  --bdg \
  --SPMR \
  --outdir peaks/
```

---

## 9. Advanced: Consensus Peak Generation

After individual MACS peak calling, nf-core creates **consensus peaks**:

```bash
# For each peak across replicates:
# 1. Intersect peaks using bedtools
bedtools intersect -a rep1_peaks.bed -b rep2_peaks.bed -f 0.5 -r

# 2. Merge overlapping peaks
bedtools merge -i all_peaks.bed

# 3. Filter by min_reps_consensus
# Keep peaks present in ≥ min_reps_consensus replicates
```

**Configurable:**

```bash
--min_reps_consensus 2  # Peak must be in at least 2 replicates
```

---

## 10. Summary: Key Takeaways

### **Default nf-core MACS3 Settings:**

| Setting | Value |
|---------|-------|
| **Peak mode** | Narrow (for TF), Broad (for histone) |
| **FDR cutoff** | 0.05 (q-value) |
| **Genome size** | Auto-detected from `--genome` |
| **Duplicate handling** | `--keep-dup all` (pre-filtered by nf-core) |
| **Summit calling** | `--call-summits` (narrow peaks only) |
| **bedGraph output** | Optional (`--save_macs_pileup`) |

### **When to Override:**

1. **Stricter peaks:** `--macs_fdr 0.01` (more stringent)
2. **Non-standard organism:** `--macs_gsize 1.2e9` (custom size)
3. **Failed auto-detection:** Add `--nomodel --extsize 200` in custom config
4. **Very broad marks:** `--broad_cutoff 0.05` (more lenient)

---

## 11. Where to Find More Details

- **nf-core documentation:** <https://nf-co.re/chipseq/2.1.0/parameters>
- **MACS3 documentation:** <https://github.com/macs3-project/MACS>
- **nf-core source code:** <https://github.com/nf-core/chipseq/blob/master/modules/nf-core/macs3/callpeak/main.nf>

---

**End of nf-core MACS Parameters Guide**
