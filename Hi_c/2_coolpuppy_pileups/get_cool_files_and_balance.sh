#!/bin/bash

input_dir="/mnt/data/pairtools_run/cooler_files/"

resolutions=("15000" "30000" "120000")

#for mcool_file in "$input_dir"*.mcool; do
#    # Extracting filename without extension
#    filename=$(basename "$mcool_file" .mcool)
#    # Extracting filename part until third underscore
#    output_prefix=$(echo "$filename" | awk -F'_' '{print $1"_"$2"_"$3}')
#    
#    for resolution in "${resolutions[@]}"; do
#        output_file="${input_dir}${output_prefix}_${resolution}.cool"
#        cooler dump "${mcool_file}::/resolutions/${resolution}" > "$output_file"
#    done
#done

# Iterate over each .cool file in the directory
for cool_file in "$input_dir"*.cool; do
    # Check if the file exists and is a regular file
    if [[ -f "$cool_file" ]]; then
        # Run cooler balance command on the current .cool file
        cooler balance "$cool_file"
    fi
done