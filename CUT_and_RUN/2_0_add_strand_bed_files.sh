#!/bin/bash
#Running this file is essential to create the upset plots from scripts 3 and 3_1

for file in ./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/seacr/*.stringent.bed; do
    new_file="${file%.stringent.bed}_stranded.bed"
    awk 'BEGIN{OFS="\t"} {$6="."}1' "$file" > "$new_file"
done


for file in ./CUT_and_RUN/1_2_Snakemake_run_merged_replicates/seacr/*.stringent.bed; do
    new_file="${file%.stringent.bed}_stranded.bed"
    awk 'BEGIN{OFS="\t"} {$6="."}1' "$file" > "$new_file"
done