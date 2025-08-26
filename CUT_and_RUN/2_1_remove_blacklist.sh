#!/bin/bash


mkdir ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/

bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/WT_ac_downsampled.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/WT_ac_filtered.stringent.bed
bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/C5_ac_downsampled.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/C5_ac_filtered.stringent.bed
bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/C9_ac_downsampled.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/C9_ac_filtered.stringent.bed

bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/WT_me3.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/WT_me3_filtered.stringent.bed
bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/C5_me3.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/C5_me3_filtered.stringent.bed
bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/C9_me3.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/C9_me3_filtered.stringent.bed

bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/WT_ac_1.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/WT_ac_1_filtered.stringent.bed
bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/WT_ac_2.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/WT_ac_2_filtered.stringent.bed
bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/C5_ac_1.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/C5_ac_1_filtered.stringent.bed
bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/C5_ac_2.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/C5_ac_2_filtered.stringent.bed
bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/C9_ac_1.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/C9_ac_1_filtered.stringent.bed
bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/C9_ac_2.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/C9_ac_2_filtered.stringent.bed


bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/WT_me3_1.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/WT_me3_1_filtered.stringent.bed
bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/WT_me3_2.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/WT_me3_2_filtered.stringent.bed
bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/C5_me3_1.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/C5_me3_1_filtered.stringent.bed
bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/C5_me3_2.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/C5_me3_2_filtered.stringent.bed
bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/C9_me3_1.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/C9_me3_1_filtered.stringent.bed
bedtools subtract -a ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/C9_me3_2.stringent.bed.stringent.bed -b ../genome/hg38-blacklist.v2.bed -A > ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/filtered/C9_me3_2_filtered.stringent.bed