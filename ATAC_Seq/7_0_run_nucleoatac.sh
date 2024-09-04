#!/bin/bash

bedtools merge ./ATAC_Seq/bwa/merged_replicate/hmmratac_peaks/WT_ATAC_peaks.bed ./ATAC_Seq/bwa/merged_replicate/hmmratac_peaks/C5_ATAC_peaks.bed ./ATAC_Seq/bwa/merged_replicate/hmmratac_peaks/C9_ATAC_peaks.bed | sort -k1,1 -k2,2n > ./ATAC_Seq/bwa/merged_replicate/hmmratac_peaks/C5_C9_WT_merged.bed

mkdir ./ATAC_Seq/nucleoatac
nucleoatac run --bed ./ATAC_Seq/bwa/merged_replicate/hmmratac_peaks/C5_C9_WT_merged.bed --bam ./ATAC_Seq/bwa/merged_replicate/WT.mRp.clN.sorted.bam --fasta ../genome/genome.fa --out ./ATAC_Seq/nucleoatac/WT_nucleoatac --cores 28
nucleoatac run --bed ./ATAC_Seq/bwa/merged_replicate/hmmratac_peaks/C5_C9_WT_merged.bed --bam ./ATAC_Seq/bwa/merged_replicate/C5.mRp.clN.sorted.bam --fasta ../genome/genome.fa --out ./ATAC_Seq/nucleoatac/C5_nucleoatac --cores 28
nucleoatac run --bed ./ATAC_Seq/bwa/merged_replicate/hmmratac_peaks/C5_C9_WT_merged.bed --bam ./ATAC_Seq/bwa/merged_replicate/C9.mRp.clN.sorted.bam --fasta ../genome/genome.fa --out ./ATAC_Seq/nucleoatac/C9_nucleoatac --cores 28
