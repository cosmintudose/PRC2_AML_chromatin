#!/bin/bash
#written with ChatGPT
# Loop over *_1_deduped.bam files in subdirectories
for bam1 in */*_1_deduped.bam; do
    # Extract base sample name: WT_me3 from WT_me3_1_deduped.bam
    filename=$(basename "$bam1")
    base=$(echo "$filename" | sed 's/_1_deduped\.bam//')

    # Build expected _2 file path using the same base
    bam2_dir="${base}_2"
    bam2="${bam2_dir}/${base}_2_deduped.bam"

    # Output directory and filename
    output_dir="${base}"
    output_file="${output_dir}/${base}_deduped.bam"

    # Check if bam2 exists
    if [[ -f "$bam2" ]]; then
        echo "Merging:"
        echo " - $bam1"
        echo " - $bam2"
        echo " => $output_file"

        mkdir -p "$output_dir"
        samtools merge "$output_file" "$bam1" "$bam2"
    else
        echo "No matching _2 file found for $base — skipping."
    fi
done
