# ChIP-seq Pipeline: Basic Commands (No Loops)

This guide shows the **exact commands** for every step of the pipeline. These are run for a single sample (e.g., `Sample1`), except for the "Compare Samples" section which is run once for all files.

### 1. Alignment (Bowtie2)

```bash
mkdir -p bowtie_align

# Align to reference
bowtie2 -x trim/ssindex \
    -1 trim/Sample1_R1_val_1.fq.gz \
    -2 trim/Sample1_R2_val_2.fq.gz \
    -p 6 --no-unal \
    -S bowtie_align/Sample1.sam
```

### 2. Sorting and Indexing (Samtools)

```bash
# Sort and convert to BAM
samtools sort -@ 6 -o bowtie_align/Sample1.sorted.bam bowtie_align/Sample1.sam

# Index
samtools index bowtie_align/Sample1.sorted.bam
```

### 3. Mark Duplicates (Picard)
*Mark duplicates but keep them in the file.*

```bash
mkdir -p Marked_duplicate3

picard MarkDuplicates \
    I=bowtie_align/Sample1.sorted.bam \
    O=Marked_duplicate3/Sample1.marked.bam \
    M=Marked_duplicate3/Sample1.marked_metrics.txt \
    REMOVE_DUPLICATES=false \
    VALIDATION_STRINGENCY=SILENT

samtools index Marked_duplicate3/Sample1.marked.bam
```

### 4. Remove Duplicates (Picard)
*Remove duplicates to create the final filtered file.*

```bash
mkdir -p de_duplicate4

picard MarkDuplicates \
    I=bowtie_align/Sample1.sorted.bam \
    O=de_duplicate4/Sample1.dedup.bam \
    M=de_duplicate4/Sample1.dedup_metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=SILENT

samtools index de_duplicate4/Sample1.dedup.bam
```

### 5. QC Metrics
*Generate mapping statistics, PBC complexity metrics, and Cross-Correlation.*

```bash
mkdir -p bam_QC_PBC

# 1. General Stats
samtools flagstat de_duplicate4/Sample1.dedup.bam > bam_QC_PBC/Sample1.flagstat.txt
samtools stats de_duplicate4/Sample1.dedup.bam > bam_QC_PBC/Sample1.stats.txt

# 2. Library Complexity (PBC/NRF)
#    a. Convert BAM to BED
bedtools bamtobed -i de_duplicate4/Sample1.dedup.bam \
    | awk 'BEGIN{OFS="\t"} ($6=="+"){print $1,$2,$2+1} ($6=="-"){print $1,$3-1,$3}' \
    | sort -k1,1 -k2,2n \
    > bam_QC_PBC/Sample1.read5.bed

#    b. Calculate Metrics
uniq -c bam_QC_PBC/Sample1.read5.bed \
    | awk '{
        c=$1; total+=c; uniq++; if(c==1) single++; if(c==2) double++;
      }
      END{
        if(total>0) {
            NRF=uniq/total; 
            PBC1=single/uniq; 
            PBC2=(double? single/double:"Inf");
        } else { NRF=0; PBC1=0; PBC2=0; }
        printf("NRF=%.3f\tPBC1=%.3f\tPBC2=%s\n", NRF, PBC1, PBC2);
      }' > bam_QC_PBC/Sample1.pbc.txt

# 3. Cross-Correlation (SPP/Phantompeakqualtools)
Rscript /opt/anaconda3/envs/chip/bin/run_spp.R \
    -c=de_duplicate4/Sample1.dedup.bam \
    -savp=bam_QC_PBC/Sample1_avp.pdf \
    -out=bam_QC_PBC/Sample1_spp.qc.txt
```

### 6. BigWig Generation (bamCoverage)

```bash
mkdir -p deeptools_QC

bamCoverage \
    -b de_duplicate4/Sample1.dedup.bam \
    -o deeptools_QC/Sample1.bw \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2652783500 \
    --binSize 10 \
    --ignoreDuplicates
```

### 7. Global Analysis (Compare Samples)
*Run these commands once you have bam files for ALL samples.*

```bash
# 1. Create a summary of all samples
#    (List all your bam files after -b)
multiBamSummary bins \
    -b de_duplicate4/Sample1.dedup.bam de_duplicate4/Sample2.dedup.bam \
    --labels Sample1 Sample2 \
    -o deeptools_QC/all_samples_bins.npz \
    --outRawCounts deeptools_QC/all_samples_bins.tab

# 2. Correlation Heatmap
plotCorrelation \
    -in deeptools_QC/all_samples_bins.npz \
    --corMethod spearman \
    --skipZeros \
    --whatToPlot heatmap \
    --colorMap RdYlBu \
    --plotNumbers \
    -o deeptools_QC/heatmap_SpearmanCorr.pdf

# 3. PCA Plot
plotPCA \
    -in deeptools_QC/all_samples_bins.npz \
    -o deeptools_QC/PCA.pdf \
    -T "PCA of Binned Coverage" \
    --transpose

# 4. Fingerprint Plot
plotFingerprint \
    -b de_duplicate4/Sample1.dedup.bam de_duplicate4/Sample2.dedup.bam \
    --labels Sample1 Sample2 \
    --skipZeros \
    --plotFile deeptools_QC/fingerprints.pdf
```

### 8. Profile Plotting (DeepTools)

```bash
# 1. Prepare Region Files from GTF (Run once)
#    TSS (Transcription Start Sites)
awk 'BEGIN{OFS="\t"} $3=="transcript" { 
      if($7=="+") print $1, $4-1, $4, $12, ".", $7; 
      else        print $1, $5-1, $5, $12, ".", $7 
    }' gencode.v49.annotation.gtf | tr -d '";' | sort -k1,1V -k2,2n > deeptools_QC/tss.bed

#    Gene Bodies
awk 'BEGIN{OFS="\t"} $3=="gene" { 
      print $1, $4-1, $5, $10, ".", $7 
    }' gencode.v49.annotation.gtf | tr -d '";' | sort -k1,1V -k2,2n > deeptools_QC/genes.bed


# 2. TSS Profile
computeMatrix reference-point \
    --referencePoint TSS \
    -b 3000 -a 10000 \
    -R deeptools_QC/tss.bed \
    -S deeptools_QC/Sample1.bw \
    -o deeptools_QC/matrix_all_bw_TSS.gz \
    --skipZeros

plotProfile \
    --matrixFile deeptools_QC/matrix_all_bw_TSS.gz \
    --outFileName deeptools_QC/all_samples_TSS_profile.pdf \
    --perGroup --refPointLabel TSS --plotType se


# 3. Gene Body Profile
computeMatrix scale-regions \
    -R deeptools_QC/genes.bed \
    -S deeptools_QC/Sample1.bw \
    --regionBodyLength 5000 \
    -b 1000 -a 1000 \
    -o deeptools_QC/matrix_all_bw_scalar.gz \
    --skipZeros

plotProfile \
    --matrixFile deeptools_QC/matrix_all_bw_scalar.gz \
    --outFileName deeptools_QC/all_samples_scalar_profile.pdf \
    --perGroup --plotType se
```
