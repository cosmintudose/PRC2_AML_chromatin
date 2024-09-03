#!/bin/bash

mkdir heatmaps
mkdir ./heatmaps/plots
##########TSS h3k27me3
computeMatrix reference-point -S ./CUT_and_RUN/1_Snakemake_run/bigwigs/WT_h3k27me3_1_deduped.bigwig ./CUT_and_RUN/1_Snakemake_run/bigwigs/C9_h3k27me3_1_deduped.bigwig \
    --referencePoint TSS -R /mnt/data/hg38.knownGene.gtf \
    -p max/2 \
    -bl ./genome/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_tss_me3_cut_run.mat.gz

plotHeatmap -m ./heatmaps/heatmap_tss_me3_cut_run.mat.gz \
    -out ./heatmaps/plots/heatmap_tss_me3_cut_run.pdf \
    --refPointLabel=TSS \
    --samplesLabel "WT H3K27me3" "C9 H3K27me3"  \
    --regionsLabel Genes \
    --missingDataColor=1 \
    --legendLocation none \
    --xAxisLabel "Distance (bp)" \
    --colorMap Greens Greens \
    --zMax 0.3 0.3 \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=10


#####TSS h2ak119ub
computeMatrix reference-point -S ./CUT_and_RUN/1_Snakemake_run/bigwigs/WT_h2ak119ub_1_deduped.bigwig ./CUT_and_RUN/1_Snakemake_run/bigwigs/C9_h2ak119ub_1_deduped.bigwig \
    --referencePoint TSS -R /mnt/data/hg38.knownGene.gtf \
    -p max/2 \
    -bl ./genome/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_tss_ub_cut_run.mat.gz


plotHeatmap -m ./heatmaps/heatmap_tss_ub_cut_run.mat.gz \
    -out ./heatmaps/plots/heatmap_tss_ub_cut_run.pdf \
    --refPointLabel=TSS \
    --samplesLabel "WT H2AK119Ub" "C9 H2AK119Ub" \
    --regionsLabel Genes \
    --missingDataColor=1 \
    --legendLocation none \
    --xAxisLabel "Distance (bp)" \
    --colorMap Purples Purples \
    --zMax 0.3 0.3 \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=10


##########Called peaks
computeMatrix reference-point -S ./CUT_and_RUN/1_Snakemake_run/bigwigs/WT_h3k27me3_1_deduped.bigwig ./CUT_and_RUN/1_Snakemake_run/bigwigs/C9_h3k27me3_1_deduped.bigwig \
    --referencePoint center -R ./CUT_and_RUN/1_Snakemake_run/seacr/WT_h3k27me3.stringent.bed ./CUT_and_RUN/1_Snakemake_run/seacr/C9_h3k27me3.stringent.bed \
    -p max/2 \
    -bl ./genome/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_me3_at_called_peak_sites.mat.gz

plotHeatmap -m ./heatmaps/heatmap_me3_at_called_peak_sites.mat.gz \
    -out ./heatmaps/plots/heatmap_me3_at_called_peak_sites.pdf \
    --refPointLabel=Center \
    --samplesLabel WT_me3 C9_me3 \
    --regionsLabel WT_peaks C9_peaks \
    --missingDataColor=1 \
    --legendLocation none \
    --colorMap=Greens \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=18


##########Called peaks h2ak119ub
computeMatrix reference-point -S ./CUT_and_RUN/1_Snakemake_run/bigwigs/WT_h2ak119ub_1_deduped.bigwig ./CUT_and_RUN/1_Snakemake_run/bigwigs/C9_h2ak119ub_1_deduped.bigwig \
    --referencePoint center -R ./CUT_and_RUN/1_Snakemake_run/seacr/WT_h2ak119ub.stringent.bed ./CUT_and_RUN/1_Snakemake_run/seacr/C9_h2ak119ub.stringent.bed \
    -p max/2 \
    -bl ./genome/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_ub_at_called_peak_sites.mat.gz

plotHeatmap -m ./heatmaps/heatmap_ub_at_called_peak_sites.mat.gz \
    -out ./heatmaps/plots/heatmap_ub_at_called_peak_sites.pdf \
    --refPointLabel=Center \
    --samplesLabel WT_ub C9_ub \
    --regionsLabel WT_peaks C9_peaks \
    --missingDataColor=1 \
    --legendLocation none \
    --colorMap=Greens \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=18


###########MOLM13 peaks
computeMatrix reference-point -S ./CUT_and_RUN/1_Snakemake_run/bigwigs/WT_h3k27me3_1_deduped.bigwig ./CUT_and_RUN/1_Snakemake_run/bigwigs/C9_h3k27me3_1_deduped.bigwig \
    --referencePoint center -R ./publicly_available_data/GSM6893214_H3K27Tri-9733S_Input-q05_peaks.broadPeak \
    -p max/2 \
    -bl ./genome/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_me3_molm13_peaks.mat.gz


plotHeatmap -m ./heatmaps/heatmap_me3_molm13_peaks.mat.gz \
    -out ./heatmaps/plots/heatmap_me3_molm13_peaks.pdf \
    --refPointLabel=Peak \
    --samplesLabel "WT H3K27me3" "C9 H3K27me3"  \
    --regionsLabel "MOLM13 GSM6893214 H3K27me3 peaks" \
    --xAxisLabel "Distance (bp)" \
    --missingDataColor=1 \
    --legendLocation none \
    --colorMap=Greens \
    --zMax 0.7 \
    --dpi=800 \
    --heatmapWidth=3 \
    --heatmapHeight=10


########HL60 peaks from ENCODE
computeMatrix reference-point -S ./CUT_and_RUN/1_Snakemake_run/bigwigs/WT_h3k27me3_1_deduped.bigwig ./CUT_and_RUN/1_Snakemake_run/bigwigs/C9_h3k27me3_1_deduped.bigwig  \
    --referencePoint center -R ./publicly_available_data/GSM5330834_ENCFF432FAX_pseudoreplicated_peaks_GRCh38.bed \
    -p max/2 \
    -bl /mnt/data/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_me3_hl60_encode_peaks.mat.gz


plotHeatmap -m ./heatmaps/heatmap_me3_hl60_encode_peaks.mat.gz \
    -out ./heatmaps/plots/heatmap_me3_hl60_encode_peaks.pdf \
    --refPointLabel=Peak \
    --samplesLabel "WT H3K27me3" "C9 H3K27me3" \
    --regionsLabel "HL60 ENCODE H3K27me3 peaks" \
    --missingDataColor=1 \
    --xAxisLabel "Distance (bp)" \
    --legendLocation none \
    --colorMap=Greens \
    --zMax 0.7 \
    --dpi=800 \
    --heatmapWidth=3 \
    --heatmapHeight=10


########ub and me3 signal at called me3 peaks
computeMatrix reference-point -S ./CUT_and_RUN/1_Snakemake_run/bigwigs/WT_h3k27me3_1_deduped.bigwig ./CUT_and_RUN/1_Snakemake_run/bigwigs/C9_h3k27me3_1_deduped.bigwig ./CUT_and_RUN/1_Snakemake_run/bigwigs/WT_h2ak119ub_1_deduped.bigwig ./CUT_and_RUN/1_Snakemake_run/bigwigs/C9_h2ak119ub_1_deduped.bigwig \
    --referencePoint center -R ./CUT_and_RUN/1_Snakemake_run/seacr/contrasts/WT_C9_me3_ub_overlap_peaks.bed ./CUT_and_RUN/1_Snakemake_run/seacr/contrasts/WT_me3_unique_peaks.bed /mnt/data/cut_and_run_backup/seacr/contrasts/C9_me3_unique_peaks.bed  \
    -p max/2 \
    -bl ./genome/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_me3_ub_and_wt_at_c9_me3_peaks.mat.gz

plotHeatmap -m ./heatmaps/heatmap_me3_ub_and_wt_at_c9_me3_peaks.mat.gz \
    -out ./heatmaps/plots/heatmap_me3_ub_and_wt_at_c9_me3_peaks.pdf \
    --refPointLabel=Center \
    --samplesLabel "WT H3K27me3" "C9 H3K27me3" "WT H2AK119Ub" "C9 H2AK119Ub" \
    --regionsLabel "WT and C9"  "WT-only" "C9-only" \
    --missingDataColor=1 \
    --xAxisLabel "Distance (bp)" \
    --colorMap Greens Greens Purples Purples \
    --zMax 0.8 0.8 0.3 0.3 \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=14

