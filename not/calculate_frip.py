#!/usr/bin/env python3
"""
Compute FRiP (Fraction of Reads in Peaks) score for ChIP-seq BAM files.

Usage:
    # Single replicate
    python calculate_frip.py -b sample.bam -p peaks.narrowPeak
    
    # Multiple replicates
    python calculate_frip.py -b rep1.bam rep2.bam -p peaks.narrowPeak
    
    # Verbose output with interpretation
    python calculate_frip.py -b sample.bam -p peaks.narrowPeak --verbose
    
    # Save to file
    python calculate_frip.py -b sample.bam -p peaks.narrowPeak -o frip_results.tsv

Examples:
    python calculate_frip.py \\
        -b encode_bam/H3K9ac_ENCFF534IPX.bam \\
        -p macs3_results/H3K9ac_ENCFF534IPX_peaks.narrowPeak \\
        --verbose
"""

import argparse
import sys
import pysam
import numpy as np
from deeptools import countReadsPerBin


def get_mapped_reads(bam_path):
    """Get total number of mapped reads from BAM file."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    mapped = bam.mapped
    bam.close()
    if mapped == 0:
        raise ValueError(f"No mapped reads found in {bam_path}")
    return mapped


def interpret_frip(frip_score):
    """Interpret FRiP score based on ENCODE guidelines."""
    if frip_score >= 0.20:
        return "✅ EXCELLENT - High quality ChIP-seq (FRiP ≥ 20%)"
    elif frip_score >= 0.05:
        return "✓ GOOD - Acceptable quality (FRiP ≥ 5%)"
    elif frip_score >= 0.01:
        return "⚠ MARGINAL - Low but potentially valid (FRiP ≥ 1%)"
    else:
        return "❌ POOR - Very low enrichment (FRiP < 1%)"


def compute_frip(bam_files, bed_file, threads, verbose=False):
    """
    Compute FRiP scores for one or more BAM files.
    
    Parameters
    ----------
    bam_files : list of str
        Paths to BAM files
    bed_file : str
        Path to peak file (BED or narrowPeak format)
    threads : int
        Number of threads for processing
    verbose : bool
        Print progress messages
        
    Returns
    -------
    dict
        Dictionary mapping BAM path to FRiP score
    """
    if verbose:
        print(f"Processing {len(bam_files)} BAM file(s)...", file=sys.stderr)
        print(f"Peak file: {bed_file}", file=sys.stderr)
        print(f"Threads: {threads}\n", file=sys.stderr)
    
    # Count reads in peaks
    if verbose:
        print("Counting reads in peaks using deepTools...", file=sys.stderr)
    
    cr = countReadsPerBin.CountReadsPerBin(
        bam_files,
        bedFile=[bed_file],
        numberOfProcessors=threads
    )
    
    reads_at_peaks = cr.run()
    reads_at_peaks = np.asarray(reads_at_peaks)
    
    # Sum reads across all peaks for each BAM
    total_reads_in_peaks = reads_at_peaks.sum(axis=0)
    
    # Calculate FRiP for each BAM
    frip_scores = {}
    for i, bam in enumerate(bam_files):
        if verbose:
            print(f"\nProcessing: {bam}", file=sys.stderr)
        
        mapped_reads = get_mapped_reads(bam)
        frip = float(total_reads_in_peaks[i]) / mapped_reads
        frip_scores[bam] = frip
        
        if verbose:
            print(f"  Reads in peaks: {int(total_reads_in_peaks[i]):,}", file=sys.stderr)
            print(f"  Total mapped reads: {mapped_reads:,}", file=sys.stderr)
            print(f"  FRiP: {frip:.6f} ({frip*100:.2f}%)", file=sys.stderr)
            print(f"  {interpret_frip(frip)}", file=sys.stderr)
    
    return frip_scores


def main():
    parser = argparse.ArgumentParser(
        description="Compute FRiP score for ChIP-seq BAM files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single replicate
  %(prog)s -b sample.bam -p peaks.narrowPeak
  
  # Multiple replicates with verbose output
  %(prog)s -b rep1.bam rep2.bam -p peaks.narrowPeak --verbose
  
  # Save results to file
  %(prog)s -b sample.bam -p peaks.narrowPeak -o frip.tsv
        """
    )
    
    parser.add_argument(
        "-b", "--bam",
        nargs="+",
        required=True,
        help="Input BAM file(s) - can specify multiple replicates"
    )
    
    parser.add_argument(
        "-p", "--peaks",
        required=True,
        help="Peak regions in BED/narrowPeak format"
    )
    
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=4,
        help="Number of threads (default: 4)"
    )
    
    parser.add_argument(
        "-o", "--out",
        default=None,
        help="Optional output TSV file"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print progress messages and quality interpretation"
    )
    
    args = parser.parse_args()
    
    try:
        # Compute FRiP scores
        frip_scores = compute_frip(
            args.bam,
            args.peaks,
            args.threads,
            verbose=args.verbose
        )
        
        # Print results to stdout (TSV format)
        if args.verbose:
            print("\n" + "="*60, file=sys.stderr)
            print("SUMMARY", file=sys.stderr)
            print("="*60 + "\n", file=sys.stderr)
        
        print("BAM\tFRiP")
        for bam, score in frip_scores.items():
            print(f"{bam}\t{score:.6f}")
        
        # Save to file if requested
        if args.out:
            with open(args.out, "w") as fh:
                fh.write("bam\tfrip\n")
                for bam, score in frip_scores.items():
                    fh.write(f"{bam}\t{score:.6f}\n")
            
            if args.verbose:
                print(f"\nResults saved to: {args.out}", file=sys.stderr)
        
        return 0
        
    except Exception as e:
        print(f"\nERROR: {str(e)}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
