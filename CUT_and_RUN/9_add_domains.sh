#!/bin/bash

cd ../seacr/contrasts

#mkdir domains

for d in 5000 10000 15000 30000 60000 120000; do
  for file in *_peaks.bed; do
    base=$(basename "$file" _peaks.bed)
    bedtools sort -i "$file" | bedtools merge -d "$d" -i - > "domains/${base}_${d}_bp_domains.bed"
  done
done