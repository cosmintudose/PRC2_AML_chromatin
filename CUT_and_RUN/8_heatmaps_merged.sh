#!/bin/bash

mkdir heatmaps
mkdir ./heatmaps/plots


#This was run on a server - Replace paths as needed (e.g. "/home/administrator")
##########TSS h3k27me3
computeMatrix reference-point -S ./bigwigs/WT_me3_deduped.bigwig ./bigwigs/C5_me3_deduped.bigwig ./bigwigs/C9_me3_deduped.bigwig \
    --referencePoint TSS -R /home/administrator/hg38.knownGene.gtf \
    -p max/2 \
    -bl /home/administrator/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_tss_me3_cut_run.mat.gz


plotHeatmap -m ./heatmaps/heatmap_tss_me3_cut_run.mat.gz \
    -out ./heatmaps/plots/heatmap_tss_me3_cut_run.pdf \
    --refPointLabel=TSS \
    --samplesLabel "WT H3K27me3" "C5 H3K27me3" "C9 H3K27me3"  \
    --regionsLabel Genes \
    --missingDataColor=1 \
    --legendLocation none \
    --xAxisLabel "Distance (bp)" \
    --colorMap Greens Greens Greens \
    --zMax 0.1 0.1 0.1 \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=8

##########TSS h3k27ac
computeMatrix reference-point -S ./bigwigs/WT_ac_deduped.bigwig ./bigwigs/C5_ac_deduped.bigwig ./bigwigs/C9_ac_deduped.bigwig \
    --referencePoint TSS -R /home/administrator/hg38.knownGene.gtf \
    -p max/2 \
    -bl /home/administrator/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_tss_ac_cut_run.mat.gz

plotHeatmap -m ./heatmaps/heatmap_tss_ac_cut_run.mat.gz \
    -out ./heatmaps/plots/heatmap_tss_ac_cut_run.pdf \
    --refPointLabel=TSS \
    --samplesLabel "WT H3K27ac" "C5 H3K27ac" "C9 H3K27ac" \
    --regionsLabel Genes \
    --missingDataColor=1 \
    --legendLocation none \
    --xAxisLabel "Distance (bp)" \
    --colorMap Oranges \
    --zMax 0.1 0.1 0.1 \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=8


###########MOLM13 peaks
computeMatrix reference-point -S ./bigwigs/WT_me3_deduped.bigwig ./bigwigs/C5_me3_deduped.bigwig ./bigwigs/C9_me3_deduped.bigwig \
    --referencePoint center -R /home/administrator/publicly_available_data/GSM6893214_H3K27Tri-9733S_Input-q05_peaks.broadPeak \
    -p max/2 \
    -bl /home/administrator/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_me3_molm13_peaks.mat.gz


plotHeatmap -m ./heatmaps/heatmap_me3_molm13_peaks.mat.gz \
    -out ./heatmaps/plots/heatmap_me3_molm13_peaks.pdf \
    --refPointLabel=Peak \
    --samplesLabel "WT H3K27me3" "C5 H3K27me3" "C9 H3K27me3"  \
    --regionsLabel "MOLM13 GSM6893214 H3K27me3 peaks" \
    --xAxisLabel "Distance (bp)" \
    --missingDataColor=1 \
    --legendLocation none \
    --colorMap=Greens \
    --zMax 0.3 \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=12


########HL60 peaks from ENCODE
computeMatrix reference-point -S ./bigwigs/WT_me3_deduped.bigwig ./bigwigs/C5_me3_deduped.bigwig ./bigwigs/C9_me3_deduped.bigwig \
    --referencePoint center -R /home/administrator/publicly_available_data/GSM5330834_ENCFF432FAX_pseudoreplicated_peaks_GRCh38.bed \
    -p max/2 \
    -bl /home/administrator/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_me3_hl60_encode_peaks.mat.gz


plotHeatmap -m ./heatmaps/heatmap_me3_hl60_encode_peaks.mat.gz \
    -out ./heatmaps/plots/heatmap_me3_hl60_encode_peaks.pdf \
    --refPointLabel=Peak \
    --samplesLabel "WT H3K27me3" "C5 H3K27me3" "C9 H3K27me3" \
    --regionsLabel "HL60 ENCODE H3K27me3 peaks" \
    --missingDataColor=1 \
    --xAxisLabel "Distance (bp)" \
    --legendLocation none \
    --colorMap=Greens \
    --zMax 0.3 \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=12

####HL60 peaks acetyl
computeMatrix reference-point -S ./bigwigs/WT_ac_deduped.bigwig ./bigwigs/C5_ac_deduped.bigwig ./bigwigs/C9_ac_deduped.bigwig \
    --referencePoint center -R /home/administrator/publicly_available_data/GSE167788_ENCFF763UAG_replicated_peaks_GRCh38.bed \
    -p max/2 \
    -bl /home/administrator/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_ac_hl60_peaks.mat.gz


plotHeatmap -m ./heatmaps/heatmap_ac_hl60_peaks.mat.gz \
    -out ./heatmaps/plots/heatmap_ac_hl60_peaks.pdf \
    --refPointLabel=Peak \
    --samplesLabel "WT H3K27Ac" "C5 H3K27Ac" "C9 H3K27Ac"  \
    --regionsLabel "HL60 GSE167788 H3K27ac peaks" \
    --xAxisLabel "Distance (bp)" \
    --missingDataColor=1 \
    --legendLocation none \
    --colorMap=Oranges \
    --zMax 0.1 \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=12

########ac and me3 signal at called me3 peaks
computeMatrix reference-point -S ./bigwigs/WT_me3_deduped.bigwig ./bigwigs/C5_me3_deduped.bigwig ./bigwigs/C9_me3_deduped.bigwig \
    --referencePoint center -R ./seacr/contrasts/WT_me3_unique_peaks.bed ./seacr/contrasts/C5_me3_unique_peaks.bed ./seacr/contrasts/C9_me3_unique_peaks.bed ./seacr/contrasts/c5_c9_me3_overlap_peaks.bed ./seacr/contrasts/wt_c5_c9_me3_overlap_peaks.bed \
    -p max/2 \
    -bl /home/administrator/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_me3_signal_at_me3_peaks.mat.gz

plotHeatmap -m ./heatmaps/heatmap_me3_signal_at_me3_peaks.mat.gz \
    -out ./heatmaps/plots/heatmap_me3_signal_at_me3_peaks.pdf \
    --refPointLabel=Center \
    --samplesLabel "WT H3K27me3" "C5 H3K27me3" "C9 H3K27me3" \
    --regionsLabel "WT-only" "C5-only" "C9-only" "C5-and-C9" "WT-C5-C9" \
    --missingDataColor=1 \
    --xAxisLabel "Distance (bp)" \
    --colorMap Greens Greens Greens \
    --zMax 0.3 0.3 0.3 \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=20


#########me3 extended#########
computeMatrix reference-point -S ./bigwigs/WT_me3_deduped.bigwig ./bigwigs/C5_me3_deduped.bigwig ./bigwigs/C9_me3_deduped.bigwig \
    --referencePoint center -R ./seacr/contrasts/WT_me3_unique_peaks.bed ./seacr/contrasts/C5_me3_unique_peaks.bed ./seacr/contrasts/C9_me3_unique_peaks.bed ./seacr/contrasts/c5_c9_me3_overlap_peaks.bed ./seacr/contrasts/wt_c5_c9_me3_overlap_peaks.bed \
    -p max/2 \
    -bl /home/administrator/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 10000 \
    --afterRegionStartLength 10000 \
    -o ./heatmaps/heatmap_me3_signal_at_me3_peaks_extended.mat.gz

plotHeatmap -m ./heatmaps/heatmap_me3_signal_at_me3_peaks_extended.mat.gz \
    -out ./heatmaps/plots/heatmap_me3_signal_at_me3_peaks_extended.pdf \
    --refPointLabel=Center \
    --samplesLabel "WT H3K27me3" "C5 H3K27me3" "C9 H3K27me3" \
    --regionsLabel "WT-only" "C5-only" "C9-only" "C5-and-C9" "WT-C5-C9" \
    --missingDataColor=1 \
    --xAxisLabel "Distance (bp)" \
    --colorMap Greens Greens Greens \
    --zMax 0.3 0.3 0.3 \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=20



############ac#########
computeMatrix reference-point -S ./bigwigs/WT_ac_deduped.bigwig ./bigwigs/C5_ac_deduped.bigwig ./bigwigs/C9_ac_deduped.bigwig \
    --referencePoint center -R ./seacr/WT_me3.stringent.bed.stringent.bed ./seacr/C5_me3.stringent.bed.stringent.bed ./seacr/C9_me3.stringent.bed.stringent.bed \
    -p max/2 \
    -bl /home/administrator/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_ac_signal_at_me3_full_peaks.mat.gz

plotHeatmap -m ./heatmaps/heatmap_ac_signal_at_me3_full_peaks.mat.gz \
    -out ./heatmaps/plots/heatmap_ac_signal_at_me3_full_peaks.pdf \
    --refPointLabel=Center \
    --samplesLabel "WT H3K27Ac" "C5 H3K27Ac" "C9 H3K27Ac" \
    --regionsLabel "WT H3K27me3" "C5 H3K27me3" "C9 H3K27me3" \
    --missingDataColor=1 \
    --xAxisLabel "Distance (bp)" \
    --colorMap Oranges Oranges Oranges \
    --zMax 0.05 0.05 0.05 \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=20


computeMatrix reference-point -S ./bigwigs/WT_ac_deduped.bigwig ./bigwigs/C5_ac_deduped.bigwig ./bigwigs/C9_ac_deduped.bigwig \
    --referencePoint center -R ./seacr/contrasts/WT_me3_unique_peaks.bed ./seacr/contrasts/C5_me3_unique_peaks.bed ./seacr/contrasts/C9_me3_unique_peaks.bed ./seacr/contrasts/c5_c9_me3_overlap_peaks.bed ./seacr/contrasts/wt_c5_c9_me3_overlap_peaks.bed \
    -p max/2 \
    -bl /home/administrator/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_ac_signal_at_me3_peaks.mat.gz

plotHeatmap -m ./heatmaps/heatmap_ac_signal_at_me3_peaks.mat.gz \
    -out ./heatmaps/plots/heatmap_ac_signal_at_me3_peaks.pdf \
    --refPointLabel=Center \
    --samplesLabel "WT H3K27Ac" "C5 H3K27Ac" "C9 H3K27Ac" \
    --regionsLabel "WT-only" "C5-only" "C9-only" "C5-and-C9" "WT-C5-C9" \
    --missingDataColor=1 \
    --xAxisLabel "Distance (bp)" \
    --colorMap Oranges Oranges Oranges \
    --zMax 0.3 0.3 0.3 \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=20





######me3 and ac signal at ac peaks
######## ac and me3 signal at called ac peaks
computeMatrix reference-point -S ./bigwigs/WT_me3_deduped.bigwig ./bigwigs/C5_me3_deduped.bigwig ./bigwigs/C9_me3_deduped.bigwig \
    --referencePoint center -R ./seacr/contrasts/WT_ac_unique_peaks.bed ./seacr/contrasts/C5_ac_unique_peaks.bed ./seacr/contrasts/C9_ac_unique_peaks.bed ./seacr/contrasts/c5_c9_ac_overlap_peaks.bed ./seacr/contrasts/wt_c5_c9_ac_overlap_peaks.bed \
    -p max/2 \
    -bl /home/administrator/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_me3_signal_at_ac_peaks.mat.gz

plotHeatmap -m ./heatmaps/heatmap_me3_signal_at_ac_peaks.mat.gz \
    -out ./heatmaps/plots/heatmap_me3_signal_at_ac_peaks.pdf \
    --refPointLabel=Center \
    --samplesLabel "WT H3K27me3" "C5 H3K27me3" "C9 H3K27me3" \
    --regionsLabel "WT-only" "C5-only" "C9-only" "C5-and-C9" "WT-C5-C9" \
    --missingDataColor=1 \
    --xAxisLabel "Distance (bp)" \
    --colorMap Greens Greens Greens \
    --zMax 0.3 0.3 0.3 \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=20

############ac#########
computeMatrix reference-point -S ./bigwigs/WT_ac_deduped.bigwig ./bigwigs/C5_ac_deduped.bigwig ./bigwigs/C9_ac_deduped.bigwig \
    --referencePoint center -R ./seacr/contrasts/WT_ac_unique_peaks.bed ./seacr/contrasts/C5_ac_unique_peaks.bed ./seacr/contrasts/C9_ac_unique_peaks.bed ./seacr/contrasts/c5_c9_ac_overlap_peaks.bed ./seacr/contrasts/wt_c5_c9_ac_overlap_peaks.bed \
    -p max/2 \
    -bl /home/administrator/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_ac_signal_at_ac_peaks.mat.gz

plotHeatmap -m ./heatmaps/heatmap_ac_signal_at_ac_peaks.mat.gz \
    -out ./heatmaps/plots/heatmap_ac_signal_at_ac_peaks.pdf \
    --refPointLabel=Center \
    --samplesLabel "WT H3K27Ac" "C5 H3K27Ac" "C9 H3K27Ac" \
    --regionsLabel "WT-only" "C5-only" "C9-only" "C5-and-C9" "WT-C5-C9" \
    --missingDataColor=1 \
    --xAxisLabel "Distance (bp)" \
    --colorMap Oranges Oranges Oranges \
    --zMax 0.3 0.3 0.3 \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=20





#######atac seq#########
computeMatrix reference-point -S ./bigwigs/WT_me3_deduped.bigwig ./bigwigs/C5_me3_deduped.bigwig ./bigwigs/C9_me3_deduped.bigwig \
    --referencePoint center -R /home/administrator/atac_seq/hmmratac_peaks/atac_peak_calling/WT_ATAC_peaks.bed /home/administrator/atac_seq/hmmratac_peaks/atac_peak_calling/C5_ATAC_peaks.bed /home/administrator/atac_seq/hmmratac_peaks/atac_peak_calling/C9_ATAC_peaks.bed  \
    -p max/2 \
    -bl /home/administrator/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_me3_signal_at_atac_peaks.mat.gz

plotHeatmap -m ./heatmaps/heatmap_me3_signal_at_atac_peaks.mat.gz \
    -out ./heatmaps/plots/heatmap_me3_signal_at_atac_peaks.pdf \
    --refPointLabel=Center \
    --samplesLabel "WT H3K27me3" "C5 H3K27me3" "C9 H3K27me3" \
    --regionsLabel "WT open" "C5 open" "C9 open" \
    --missingDataColor=1 \
    --xAxisLabel "Distance (bp)" \
    --colorMap Greens Greens Greens \
    --zMax 0.1 0.1 0.1 \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=10



#ac
computeMatrix reference-point -S ./bigwigs/WT_ac_deduped.bigwig ./bigwigs/C5_ac_deduped.bigwig ./bigwigs/C9_ac_deduped.bigwig \
    --referencePoint center -R /home/administrator/atac_seq/hmmratac_peaks/atac_peak_calling/WT_ATAC_peaks.bed /home/administrator/atac_seq/hmmratac_peaks/atac_peak_calling/C5_ATAC_peaks.bed /home/administrator/atac_seq/hmmratac_peaks/atac_peak_calling/C9_ATAC_peaks.bed  \
    -p max/2 \
    -bl /home/administrator/hg38-blacklist.v2.bed \
    -bs=100 \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    -o ./heatmaps/heatmap_ac_signal_at_atac_peaks.mat.gz

plotHeatmap -m ./heatmaps/heatmap_ac_signal_at_atac_peaks.mat.gz \
    -out ./heatmaps/plots/heatmap_ac_signal_at_atac_peaks.pdf \
    --refPointLabel=Center \
    --samplesLabel "WT H3K27Ac" "C5 H3K27Ac" "C9 H3K27Ac" \
    --regionsLabel "WT open" "C5 open" "C9 open" \
    --missingDataColor=1 \
    --xAxisLabel "Distance (bp)" \
    --colorMap Oranges Oranges Oranges \
    --zMax 0.1 0.1 0.1 \
    --dpi=800 \
    --heatmapWidth=4 \
    --heatmapHeight=10