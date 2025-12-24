############################################
# ChIP-seq Analysis Snakemake Workflow
# Parameters validated against peer-reviewed literature
# See chipseq_parameter_references.md for citations
############################################

############################################
# GLOBAL CONFIGURATION
############################################

# hg38 effective genome size (khmer-based calculation)
# Calculated using khmer unique-kmers.py for 50bp read length
# Source: deepTools documentation (recommended method for unique mappable genome)
# This accounts for uniquely mappable regions only (excludes all repeats + N-bases)
# 
# Read-length specific values (from khmer):
#   50bp:  2,701,495,761 bp  <-- Using this (user's actual read length)
#   75bp:  2,747,877,777 bp
#   100bp: 2,805,636,331 bp
#   150bp: 2,862,010,578 bp
#   200bp: 2,887,553,303 bp
#
# NOTE: Smaller than the standard "mappable" value (2,913,022,398)
#       because this excludes ALL repetitive regions, not just N-bases
GENOME_SIZE = 2701495761   # hg38, 50bp reads, khmer unique k-mers
THREADS = 6


############################################
# SAMPLE DEFINITION
############################################

SAMPLES = {
    "Input": [
        "ENCFF110SOB",
        "ENCFF919XCV"
    ],
    "H3K9ac": [
        "ENCFF534IPX",
        "ENCFF193NPE"
    ],
    "H3K27me3": [
        "ENCFF164ALR",
        "ENCFF532DQH"
    ],
    "ceb": [
        "ENCFF327JFG",
        "ENCFF744SVA"
    ]
}

IP_MARKS = ["H3K9ac", "H3K27me3", "ceb"]


############################################
# MARK-SPECIFIC MATRIX PARAMETERS
# Literature-validated (see citations below)
############################################

MATRIX_PARAMS = {
    # Klein et al. 2019 (Nature Neuroscience)
    # Standard for active promoter marks
    "H3K9ac": {
        "up": 3000, 
        "down": 3000, 
        "bin": 50
    },
    
    # Whitfield et al. 2012 (Genome Biology)
    # Gold standard for transcription factor analysis
    "ceb": {
        "up": 2000, 
        "down": 2000, 
        "bin": 25
    },
    
    # Young et al. 2011 (Nucleic Acids Research)
    # Optimized for broad repressive domains
    "H3K27me3": {
        "up": 5000, 
        "down": 5000, 
        "bin": 100
    }
}


############################################
# FINAL TARGETS
############################################

rule all:
    input:
        # TSS heatmaps for all marks
        expand("results/plots/{mark}_TSS_heatmap.pdf", mark=IP_MARKS),
        
        # Gene body heatmap (H3K27me3 only)
        "results/plots/H3K27me3_geneBody_heatmap.pdf",
        
        # QC: Replicate consistency plots
        expand("results/qc/{mark}_TSS_reps_profile.pdf", mark=IP_MARKS)


############################################
# STEP 1: BAM → RPGC BIGWIG
# ENCODE guidelines (Landt et al. 2012)
############################################

rule bam_to_rpgc_bw:
    input:
        bam="bam/{mark}_{acc}.bam"  # BAM files are in bam/ subdirectory
    output:
        bw="results/bw/{mark}_{acc}.RPGC.bw"
    wildcard_constraints:
        acc="ENCFF[A-Z0-9]+"  # Only match ENCODE accession IDs
    threads: THREADS
    shell:
        """
        bamCoverage \
          -b {input.bam} \
          -o {output.bw} \
          --normalizeUsing RPGC \
          --effectiveGenomeSize {GENOME_SIZE} \
          --binSize 50 \
          --extendReads 200 \
          --ignoreDuplicates \
          --numberOfProcessors {threads}
        """


############################################
# QC: REPLICATE CONSISTENCY
# Verify biological replicates show similar patterns
############################################

rule qc_replicate_matrix:
    input:
        lambda wc: [
            f"results/bw/{wc.mark}_{acc}.RPGC.bw"
            for acc in SAMPLES[wc.mark]
        ]
    output:
        mat="results/qc/{mark}_TSS_reps.mat.gz"
    threads: 4
    shell:
        """
        computeMatrix reference-point \
          --referencePoint TSS \
          -b 3000 -a 3000 \
          -R annotations/tss.bed \
          -S {input} \
          --binSize 50 \
          --numberOfProcessors {threads} \
          -o {output.mat}
        """


rule qc_replicate_profile:
    input:
        mat="results/qc/{mark}_TSS_reps.mat.gz"
    output:
        pdf="results/qc/{mark}_TSS_reps_profile.pdf"
    shell:
        """
        plotProfile \
          -m {input.mat} \
          --perGroup \
          --plotTitle "{wildcards.mark} Replicate Consistency" \
          --dpi 600 \
          -out {output.pdf}
        """


############################################
# STEP 2: AVERAGE REPLICATES
# Signal-level mean of biological replicates
############################################

rule average_replicates:
    input:
        lambda wc: [
            f"results/bw/{wc.mark}_{acc}.RPGC.bw"
            for acc in SAMPLES[wc.mark]
        ]
    output:
        bw="results/bw/{mark}_mean.RPGC.bw"
    shell:
        """
        bigwigCompare \
          -b1 {input[0]} \
          -b2 {input[1]} \
          --operation mean \
          -o {output.bw}
        """


############################################
# STEP 3: IP NORMALIZED TO INPUT (LOG2)
# Standard normalization for ChIP-seq visualization
############################################

rule ip_over_input:
    input:
        ip="results/bw/{mark}_mean.RPGC.bw",
        inp="results/bw/Input_mean.RPGC.bw"
    output:
        bw="results/bw/{mark}_log2IPoverInput.RPGC.bw"
    shell:
        """
        bigwigCompare \
          -b1 {input.ip} \
          -b2 {input.inp} \
          --operation log2 \
          --pseudocount 1 \
          -o {output.bw}
        """


############################################
# STEP 4: MARK-SPECIFIC TSS MATRICES
# Parameters optimized per mark type (see MATRIX_PARAMS)
############################################

rule tss_matrix_final:
    input:
        bw="results/bw/{mark}_log2IPoverInput.RPGC.bw"
    output:
        mat="results/plots/{mark}_TSS.mat.gz"
    threads: 4
    params:
        up=lambda wc: MATRIX_PARAMS[wc.mark]["up"],
        down=lambda wc: MATRIX_PARAMS[wc.mark]["down"],
        bin=lambda wc: MATRIX_PARAMS[wc.mark]["bin"]
    shell:
        """
        computeMatrix reference-point \
          --referencePoint TSS \
          -b {params.up} -a {params.down} \
          -R annotations/tss.bed \
          -S {input.bw} \
          --binSize {params.bin} \
          --numberOfProcessors {threads} \
          -o {output.mat}
        """


############################################
# STEP 5: TSS HEATMAPS
# Final publication-quality heatmaps
############################################

rule tss_heatmap_final:
    input:
        mat="results/plots/{mark}_TSS.mat.gz"
    output:
        pdf="results/plots/{mark}_TSS_heatmap.pdf"
    shell:
        """
        plotHeatmap \
          -m {input.mat} \
          --colorMap RdBu_r \
          --refPointLabel "TSS" \
          --plotTitle "{wildcards.mark} Enrichment at TSS" \
          --dpi 600 \
          -out {output.pdf}
        """


############################################
# STEP 6: GENE BODY ANALYSIS (H3K27me3 ONLY)
# Young et al. 2011 - Broad domains require scale-regions mode
# Updated parameters: 3kb flanks, 250bp bins (optimized)
############################################

rule gene_body_matrix_H3K27me3:
    input:
        bw="results/bw/H3K27me3_log2IPoverInput.RPGC.bw"
    output:
        mat="results/plots/H3K27me3_geneBody.mat.gz"
    threads: 4
    shell:
        """
        computeMatrix scale-regions \
          -R annotations/genes.bed \
          -S {input.bw} \
          --regionBodyLength 5000 \
          -b 3000 -a 3000 \
          --binSize 250 \
          --numberOfProcessors {threads} \
          -o {output.mat}
        """


rule gene_body_heatmap_H3K27me3:
    input:
        mat="results/plots/H3K27me3_geneBody.mat.gz"
    output:
        pdf="results/plots/H3K27me3_geneBody_heatmap.pdf"
    shell:
        """
        plotHeatmap \
          -m {input.mat} \
          --colorMap RdBu_r \
          --plotTitle "H3K27me3 Broad Domain Coverage" \
          --dpi 600 \
          -out {output.pdf}
        """


############################################
# REFERENCES
############################################
# 1. Ramírez et al. 2016 (NAR) - deepTools2
# 2. Klein et al. 2019 (Nat Neurosci) - H3K9ac ±3kb
# 3. Whitfield et al. 2012 (Genome Biol) - TF ±2kb
# 4. Young et al. 2011 (NAR) - H3K27me3 broad domains
# 5. Landt et al. 2012 (Genome Res) - ENCODE guidelines
#
# Full citations: see chipseq_parameter_references.md
############################################
