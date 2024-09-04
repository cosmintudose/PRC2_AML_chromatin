#!/bin/bash

cd ./ATAC_Seq/bwa/merged_replicate

samtools view -H WT.mRp.clN.sorted.bam | perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > WT_genome.info 
samtools view -H C5.mRp.clN.sorted.bam | perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > C5_genome.info 
samtools view -H C9.mRp.clN.sorted.bam | perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > C9_genome.info 

echo "Finished genome info files"

mkdir hmmratac_peaks

java -jar HMMRATAC_V1.2.10_exe.jar -b WT.mRp.clN.sorted.bam -i WT.mRp.clN.sorted.bam.bai -g WT_genome.info --bedgraph True --blacklist ../../../genome/hg38-blacklist.v2.bed -o ./hmmratac_peaks/WT_peaks
echo "Finished WT peak calling"

java -jar HMMRATAC_V1.2.10_exe.jar -b C5.mRp.clN.sorted.bam -i C5.mRp.clN.sorted.bam.bai -g C5_genome.info --bedgraph True --blacklist ../../../genome/hg38-blacklist.v2.bed -o ./hmmratac_peaks/C5_peaks
echo "Finished C5 peak calling"

java -jar HMMRATAC_V1.2.10_exe.jar -b C9.mRp.clN.sorted.bam -i C9.mRp.clN.sorted.bam.bai -g C9_genome.info --bedgraph True --blacklist ../../../genome/hg38-blacklist.v2.bed -o ./hmmratac_peaks/C9_peaks
echo "Finished C9 peak calling"