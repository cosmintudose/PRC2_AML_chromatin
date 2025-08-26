#!/bin/bash

cd ./1_0_Snakemake_run_sep_replicates/bigwigs/
multiBigwigSummary bins -bs 1000 -b WT_ac_1_deduped.bigwig WT_ac_2_deduped.bigwig C5_ac_1_deduped.bigwig C5_ac_2_deduped.bigwig C9_ac_1_deduped.bigwig C9_ac_2_deduped.bigwig -bl ../../../genome/hg38-blacklist.v2.bed -o ac_samples_pca.npz

multiBigwigSummary bins -bs 1000 -b WT_me3_1_deduped.bigwig WT_me3_2_deduped.bigwig C5_me3_1_deduped.bigwig C5_me3_2_deduped.bigwig C9_me3_1_deduped.bigwig C9_me3_2_deduped.bigwig -bl ../../../genome/hg38-blacklist.v2.bed -o me3_samples_pca.npz

plotPCA -in ac_samples_pca.npz --outFileNameData ac_pca.tab -o pca_ac.png
plotPCA -in me3_samples_pca.npz --outFileNameData me3_pca.tab -o pca_me3.png