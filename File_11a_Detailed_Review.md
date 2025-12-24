# Review: File 11a - MACS3 Peak Calling Tutorial

**Review Date:** 2025-12-22  
**File:** `src/11a_macs3_peak_calling.md`  
**Reviewer Focus:** Numbering, Keywords, Scripts, Flow, Synchronization with File 10

---

## ✅ 1. NUMBERING ANALYSIS

### **Current State:**

- **Filename:** `11a_macs3_peak_calling.md`
- **Internal Title:** "Tutorial 11: Peak Calling with MACS3"

### **Issues:**

| Issue | Severity | Description |
|-------|----------|-------------|
| **Missing 11b reference** | ⚠️ Medium | File is `11a` but no `11b` exists or is referenced |
| **Inconsistent with sequence** | ⚠️ Medium | Jumps from 10 → 11a → 12, skipping plain "11" |
| **File 10 reference mismatch** | ✅ Good | Line 55 correctly references "Tutorial 10" |

### **Recommendations:**

**Option 1: Rename to 11 (simple)**

```bash
# If no 11b is planned
mv 11a_macs3_peak_calling.md 11_macs3_peak_calling.md
```

**Option 2: Keep 11a, add 11b for BigWig**

```bash
# If BigWig generation is separate
11a_macs3_peak_calling.md     (peak calling)
11b_bigwig_generation.md       (signal tracks)
12_frip_quality_metrics.md     (FRiP)
```

**Verdict:** ✅ **Rename to `11_macs3_peak_calling.md`** unless you plan a dedicated BigWig section as 11b.

---

## ✅ 2. KEYWORDS ANALYSIS

### **Current Keywords (Line 3):**

```
`MACS3` `peak-calling` `ChIP-seq` `narrow-peaks` `broad-peaks` 
`transcription-factors` `histone-modifications` `FRiP` `IDR` 
`reproducibility` `ENCODE`
```

### **Assessment:**

| Aspect | Rating | Notes |
|--------|--------|-------|
| **Relevance** | ✅ Excellent | All keywords directly relate to content |
| **Comprehensiveness** | ⚠️ Good | Missing some technical terms |
| **Searchability** | ✅ Good | Covers main concepts |

### **Missing Keywords to Add:**

Based on actual file content:

- `summit` or `call-summits` (mentioned multiple times)
- `q-value` or `FDR` (key parameter)
- `--keep-dup` (important parameter explained)
- `genome-size` (discussed in context)
- `control` or `input` (critical for methodology)

### **Recommended Keywords:**

```markdown
`MACS3` `peak-calling` `ChIP-seq` `narrow-peaks` `broad-peaks` 
`transcription-factors` `histone-modifications` `summits` `q-value` 
`FDR` `control` `input` `IDR` `reproducibility` `ENCODE` `--keep-dup` `genome-size`
```

**Verdict:** ⚠️ **Add 5-6 more technical keywords** for better discoverability.

---

## ✅ 3. SCRIPT ACCURACY REVIEW

### **Script 1: H3K9ac Narrow Peaks (Lines 79-101)**

```bash
macs3 callpeak \
  -t encode_bam/H3K9ac_ENCFF534IPX.bam \
  -c encode_bam/Input_ENCFF110SOB.bam \
  -f BAM \
  -g hs \
  -n H3K9ac_ENCFF534IPX \
  -p 0.01 \
  --keep-dup all \
  --outdir macs3_results
```

**Issues Found:**

| Line | Issue | Severity | Fix |
|------|-------|----------|-----|
| 85 | Uses `-p 0.01` (p-value) | ⚠️ Medium | Should use `-q 0.01` (FDR) for consistency |
| 86 | `--keep-dup all` | ✅ Correct | Good for deduplicated BAMs |
| Missing | No `--call-summits` | ⚠️ Medium | Should add for narrow peaks |
| Missing | No `-B` or `--bdg` | ℹ️ Info | Optional but useful for bedGraph |

**Recommended Fix:**

```bash
macs3 callpeak \
  -t encode_bam/H3K9ac_ENCFF534IPX.bam \
  -c encode_bam/Input_ENCFF110SOB.bam \
  -f BAM \
  -g hs \
  -n H3K9ac_ENCFF534IPX \
  -q 0.01 \                          # Changed from -p to -q
  --keep-dup all \
  --call-summits \                    # Added for narrow peaks
  --outdir macs3_results
```

---

### **Script 2: H3K27me3 Broad Peaks (Lines 138-161)**

```bash
macs3 callpeak \
  -t encode_bam/H3K27me3_ENCFF532DQH.bam \
  -c encode_bam/Input_ENCFF110SOB.bam \
  -f BAM \
  -g hs \
  -n H3K27me3_ENCFF532DQH \
  --broad -p 0.01\
  --keep-dup all \
  --outdir macs3_results
```

**Issues Found:**

| Line | Issue | Severity | Fix |
|------|-------|----------|-----|
| 144 | `-p 0.01` for broad peaks | ⚠️ Medium | Should be `-q 0.01` for consistency |
| 144 | Missing `--broad-cutoff` | ⚠️ Medium | Default is 0.1, but should be explicit |
| 144 | Formatting issue | ✅ Minor | Space after `--broad` |

**Recommended Fix:**

```bash
macs3 callpeak \
  -t encode_bam/H3K27me3_ENCFF532DQH.bam \
  -c encode_bam/Input_ENCFF110SOB.bam \
  -f BAM \
  -g hs \
  -n H3K27me3_ENCFF532DQH \
  --broad \
  -q 0.01 \                          # Changed from -p
  --broad-cutoff 0.1 \                # Added explicitly
  --keep-dup all \
  --outdir macs3_results
```

---

### **Script 3: CEBPA TF (Lines 178-200)**

```bash
macs3 callpeak \
  -t encode_bam/ceb_ENCFF327JFG.bam \
  -c encode_bam/Input_ENCFF110SOB.bam \
  -f BAM \
  -g hs \
  -n ceb_ENCFF327JFG \
  -p 0.01 \
  --keep-dup all \
  --outdir macs3_results
```

**Same Issues:**

- Uses `-p` instead of `-q`
- Missing `--call-summits` (critical for TFs!)

**Recommended Fix:**

```bash
macs3 callpeak \
  -t encode_bam/ceb_ENCFF327JFG.bam \
  -c encode_bam/Input_ENCFF110SOB.bam \
  -f BAM \
  -g hs \
  -n ceb_ENCFF327JFG \
  -q 0.01 \                          # Changed from -p
  --keep-dup all \
  --call-summits \                    # CRITICAL for TFs!
  --outdir macs3_results
```

---

## ✅ 4. PARAMETER EXPLANATION ACCURACY

### **Line 117-127: Parameter Explanations**

**Current Explanation (Line 124):**

```markdown
* `-q 0.01`: FDR (False Discovery Rate) cutoff...
```

**Issue:** The scripts use `-p 0.01` but explanation says `-q`! **Mismatch!**

**Current Explanation (Line 125):**

```markdown
* `--keep-dup all`: for deduplicated bam files
```

**Issue:** Explanation is correct but could be clearer.

**Recommended Fix:**

```markdown
Parameters:

* `-t`: Treatment file (IP/ChIP BAM - the enriched DNA)
* `-c`: Control file (Input BAM - background DNA)
* `-f BAM`: Input file format specification
* `-g hs`: Effective genome size (`hs` = human ~2.7-2.9Gb)
* `-n`: Output file prefix (recommend using descriptive sample names)
* `-q 0.01`: Minimum FDR cutoff - keep only peaks with q-value < 0.01 (1% false discovery rate)
* `--keep-dup all`: Keep all reads (use this for deduplicated BAMs; duplicates already removed upstream)
* `--call-summits`: Find precise peak summits (narrow peaks only - TF and sharp histone marks)
* `--outdir`: Directory for all output files
```

---

## ✅ 5. FLOW & READABILITY ANALYSIS

### **Strengths:**

| Aspect | Rating | Comments |
|--------|--------|----------|
| **Section 1: Concept** | ✅ Excellent | "Heap Hunt" analogy is brilliant |
| **Section 2: Requirements** | ✅ Good | Clear file listing |
| **Section 3: Execution** | ✅ Good | Step-by-step structure |
| **Section 4: Outputs** | ✅ Good | File format explanations |
| **Visual aids** | ✅ Excellent | Screenshot at line 28 helps |

### **Flow Issues:**

| Issue | Location | Impact | Fix |
|-------|----------|--------|-----|
| **Parameter mismatch** | Lines 85, 99, 184 | ⚠️ High | Scripts say `-p`, explanation says `-q` |
| **Missing transition** | Line 60 → 61 | ℹ️ Low | Jump from concept to execution needs bridge sentence |
| **Incomplete reference** | Line 166 | ℹ️ Low | Says `--broad-cutoff 0.1` but scripts don't show it |
| **No quality checkpoint** | After Section 3 | ⚠️ Medium | Should mention "How many peaks expected?" |

### **Recommended Flow Improvements:**

**Add after line 60 (before Section 3):**

```markdown
---

## Before You Start: Quality Validation

Confirm that your BAM files have passed QC checks from Tutorial 10:

✅ **Fingerprint:** Sharp enrichment curves (not flat)  
✅ **Coverage:** Adequate sequencing depth  
✅ **Correlation:** Replicates cluster together  

If any QC failed, peak calling may produce unreliable results.

---
```

**Add after Section 3 (after line 202):**

```markdown
### Quality Checkpoint: Expected Peak Counts

After running MACS3, check if peak counts are reasonable:

| Sample Type | Expected Peak Count | Your Results |
|-------------|---------------------|--------------|
| **TF (CEBPA)** | 10,000 - 50,000 peaks | `wc -l *ceb*narrowPeak` |
| **H3K9ac** | 30,000 - 80,000 peaks | `wc -l *H3K9ac*narrowPeak` |
| **H3K27me3** | 5,000 - 20,000 broad domains | `wc -l *H3K27me3*broadPeak` |

**Too few peaks?** Poor IP enrichment or overly stringent threshold.  
**Too many peaks?** High background or lenient threshold - consider using `-q` instead of `-p`.
```

---

## ✅ 6. SYNCHRONIZATION WITH FILE 10

### **File 10 → File 11a Transition:**

**File 10 ending (line 298-300):**

```markdown
> [!NOTE]
> **Up Next:** With comprehensive QC validation complete, we're ready 
> to call peaks with MACS2 and identify protein-DNA binding sites.
```

**File 11a opening (line 5-7):**

```markdown
## 1. Basic Concept (The "Heap" Hunt)

### What is a Peak?
```

**Issues:**

| Issue | Severity | Description |
|-------|----------|-------------|
| **Says MACS2, uses MACS3** | ⚠️ Low | File 10 says "MACS2" but 11a is MACS3 |
| **No reference back to QC** | ⚠️ Medium | Doesn't remind user to check Tutorial 10 results |
| **Abrupt start** | ℹ️ Low | Jumps straight to concept without transition |

### **Recommendations:**

**Fix File 10 line 300:**

```markdown
> [!NOTE]
> **Up Next:** With comprehensive QC validation complete, we're ready 
> to call peaks with MACS3 in Tutorial 11 and identify protein-DNA binding sites.
```

**Add to File 11a after line 4 (before Section 1):**

```markdown
## Prerequisites

Before calling peaks, ensure you've completed **[Tutorial 10: deepTools QC](./10_bam_summary_fingerprint.md)** and confirmed:

- ✅ Fingerprint plots show strong IP enrichment
- ✅ Coverage is adequate (>10M usable reads per sample)
- ✅ Biological replicates cluster together in PCA/correlation

These QC metrics determine whether peak calling will be successful.

---
```

---

## ✅ 7. MISSING ENCODE/nf-core BEST PRACTICES

Based on your nf-core comparison document, file 11a is missing:

### **High-Priority Additions:**

1. **Genome Size Discussion** (currently just says `hs`)

   ```markdown
   **Genome Size Options:**
   - `-g hs` or `-g 2.7e9` (human, conservative)
   - `-g 2.9e9` (human hg38, more accurate)
   - `-g mm` or `-g 1.87e9` (mouse)
   
   **Impact:** Slightly affects peak calling stringency (~1-2% difference in peak count)
   ```

2. **--SPMR for Normalized Tracks**

   ```markdown
   **Optional: Generate normalized bedGraph:**
   Add `--SPMR --bdg` to create signal tracks for visualization
   ```

3. **Peak Capping** (before IDR)

   ```markdown
   **Optional: Cap to top 300K peaks (ENCODE standard):**
   ```bash
   sort -k7,7nr peaks.narrowPeak | head -300000 > peaks.capped.narrowPeak
   ```

   This prevents IDR memory issues with very large peak sets.

   ```

4. **Blacklist Filtering Reference**

   ```markdown
   > [!TIP]
   > **Advanced:** Before downstream analysis, consider filtering peaks 
   > overlapping ENCODE blacklist regions (repetitive/problematic genomic areas). 
   > See Tutorial 11c for details.
   ```

---

## ✅ 8. TECHNICAL ERRORS

### **Error 1: Inconsistent -p vs -q**

**Lines affected:** 85, 99, 124, 144, 158, 184, 198

**Current state:** Scripts use `-p 0.01`, explanation says `-q 0.01`

**Fix:** Change all scripts to `-q 0.01` (FDR is standard for ChIP-seq)

---

### **Error 2: Missing --call-summits**

**Lines affected:** 79-87, 93-101, 178-186, 192-200

**Current state:** Narrow peak calls don't have `--call-summits`

**Impact:** No summit positions in column 10 of narrowPeak file

**Fix:** Add `--call-summits \` to all narrow peak commands

---

### **Error 3: Missing --broad-cutoff**

**Lines affected:** 144, 158

**Current state:** `--broad` without explicit cutoff

**Impact:** Uses default 0.1, but not transparent

**Fix:** Add `--broad-cutoff 0.1 \` explicitly

---

### **Error 4: model.r Output Redirect**

**Line 272:**

```bash
Rscript macs3_results/H3K9ac_ENCFF534IPX_model.r > macs3_results/H3K9ac_ENCFF534IPX_mode.pdf
```

**Issue:** R scripts don't work with stdout redirection to PDF!

**Fix:**

```bash
# The script creates the PDF automatically in the same directory
Rscript macs3_results/H3K9ac_ENCFF534IPX_model.r

# Output will be: H3K9ac_ENCFF534IPX_model.pdf (auto-generated)
```

---

## 📋 SUMMARY & ACTION ITEMS

### **Critical Fixes (Must Do):**

1. ✅ **Change all `-p 0.01` to `-q 0.01`** (Lines 85, 99, 144, 158, 184, 198)
2. ✅ **Add `--call-summits`** to all narrow peak commands (H3K9ac, CEBPA)
3. ✅ **Add `--broad-cutoff 0.1`** to broad peak commands (H3K27me3)
4. ✅ **Fix model.r execution** (Line 272)
5. ✅ **Update File 10 reference** from "MACS2" to "MACS3"

### **High-Priority Improvements:**

6. ⚠️ **Add prerequisites section** referencing Tutorial 10 QC
7. ⚠️ **Add quality checkpoint** with expected peak counts
8. ⚠️ **Update keywords** to include technical terms
9. ⚠️ **Add genome size discussion** (-g hs vs 2.7e9 vs 2.9e9)

### **Medium-Priority Enhancements:**

10. ℹ️ **Add --SPMR/--bdg option** for bedGraph generation
11. ℹ️ **Add peak capping section** (ENCODE standard)
12. ℹ️ **Reference blacklist filtering** (to be added as 11c)
13. ℹ️ **Rename file** from 11a to 11 (if no 11b planned)

### **Low-Priority Additions:**

14. ℹ️ Add fragment length discussion (auto vs manual --extsize)
15. ℹ️ Add troubleshooting section for failed model building
16. ℹ️ Add example peak count outputs for validation

---

## ✅ VERDICT

### **Overall Quality:** ⚠️ **GOOD with Critical Issues**

**Strengths:**

- ✅ Excellent pedagogical approach (heap analogy)
- ✅ Clear structure and progression
- ✅ Good visual aids and examples
- ✅ Correct use of `--keep-dup all` for deduplicated BAMs

**Critical Issues:**

- ❌ Scripts use `-p` but explanation says `-q` (major inconsistency)
- ❌ Missing `--call-summits` for narrow peaks (loses summit info)
- ❌ Missing explicit `--broad-cutoff`
- ❌ model.r command won't work as written

**Recommendation:** **Fix critical issues before next release.** The tutorial is well-written but has parameter inconsistencies that will confuse users and produce suboptimal results.

---

**Priority Order for Fixes:**

1. Parameter consistency (-p → -q everywhere)
2. Add --call-summits to narrow peaks
3. Add --broad-cutoff to broad peaks
4. Fix model.r command
5. Add QC checkpoint section
6. Update keywords
7. Add genome size discussion

---

**End of Review**
