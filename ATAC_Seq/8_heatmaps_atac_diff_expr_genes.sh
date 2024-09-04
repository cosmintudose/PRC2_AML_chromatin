#create bigwigs
bamCoverage -b WT.mRp.clN.sorted.bam -o WT.mRp.clN.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2
bamCoverage -b C5.mRp.clN.sorted.bam -o C5.mRp.clN.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2
bamCoverage -b C9.mRp.clN.sorted.bam -o C9.mRp.clN.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2

bamCoverage -b WT.mononucleosomal.sorted.bam -o WT.mononucleosomal.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2
bamCoverage -b C5.mononucleosomal.sorted.bam -o C5.mononucleosomal.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2
bamCoverage -b C9.mononucleosomal.sorted.bam -o C9.mononucleosomal.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2

bamCoverage -b WT.open.fragments.sorted.bam -o WT.open.fragments.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2
bamCoverage -b C5.open.fragments.sorted.bam -o C5.open.fragments.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2
bamCoverage -b C9.open.fragments.sorted.bam -o C9.open.fragments.sorted.bam.bw --binSize 1 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -p=max/2

#map the transcript IDs of differentially expressed genes to GTF
grep -f ./RNA_Seq/results_files/transcript_ids_upreg_genes_c5.txt ./genome/hg38.knownGene.gtf > ./genome/c5_upreg_genes.knownGene.gtf
grep -f ./RNA_Seq/results_files/transcript_ids_upreg_genes_c9.txt ./genome/hg38.knownGene.gtf > ./genome/c9_upreg_genes.knownGene.gtf
grep -f ./RNA_Seq/results_files/transcript_ids_downreg_genes_c5.txt ./genome/hg38.knownGene.gtf > ./genome/c5_downreg_genes.knownGene.gtf
grep -f ./RNA_Seq/results_files/transcript_ids_downreg_genes_c9.txt ./genome/hg38.knownGene.gtf > ./genome/c9_downreg_genes.knownGene.gtf


#open fragments

computeMatrix reference-point -S WT.open.fragments.sorted.bam.bw C5.open.fragments.sorted.bam.bw \
    --referencePoint TSS -R ./genome/c5_downreg_genes.knownGene.gtf \
    -p max/2 \
    -bs=10 \
    --beforeRegionStartLength 1000 \
    --afterRegionStartLength 1000 \
    -o heatmap_open_fragments_c5_downreg_genes_atac.mat.gz

plotHeatmap -m heatmap_open_fragments_c5_downreg_genes_atac.mat.gz \
    -out ./ATAC_Seq/plots/heatmap_open_fragments_c5_downreg_genes_atac.pdf \
    --refPointLabel=TSS \
    --samplesLabel WT C5 \
    --regionsLabel "Genes downregulated in C5" \
    --legendLocation none \
    --averageTypeSummaryPlot mean \
    --colorMap=Blues \
    --missingDataColor=1 \
    --dpi=1200 \
    --heatmapWidth=3 \
    --heatmapHeight=10

##############################################


#open fragments

computeMatrix reference-point -S WT.open.fragments.sorted.bam.bw C5.open.fragments.sorted.bam.bw \
    --referencePoint TSS -R ./genome/c5_upreg_genes.knownGene.gtf \
    -p max/2 \
    -bs=10 \
    --beforeRegionStartLength 1000 \
    --afterRegionStartLength 1000 \
    -o heatmap_open_fragments_c5_upreg_genes_atac.mat.gz

plotHeatmap -m heatmap_open_fragments_c5_upreg_genes_atac.mat.gz \
    -out ./ATAC_Seq/plots/heatmap_open_fragments_c5_upreg_genes_atac.pdf \
    --refPointLabel=TSS \
    --samplesLabel WT C5 \
    --regionsLabel "Genes upregulated in C5" \
    --legendLocation none \
    --averageTypeSummaryPlot mean \
    --colorMap=Reds \
    --missingDataColor=1 \
    --dpi=1200 \
    --heatmapWidth=3 \
    --heatmapHeight=10

#############################################



#open fragments

computeMatrix reference-point -S WT.open.fragments.sorted.bam.bw C9.open.fragments.sorted.bam.bw \
    --referencePoint TSS -R ./genome/c9_downreg_genes.knownGene.gtf \
    -p max/2 \
    -bs=10 \
    --beforeRegionStartLength 1000 \
    --afterRegionStartLength 1000 \
    -o heatmap_open_fragments_c9_downreg_genes_atac.mat.gz

plotHeatmap -m heatmap_open_fragments_c9_downreg_genes_atac.mat.gz \
    -out ./ATAC_Seq/plots/heatmap_open_fragments_c9_downreg_genes_atac.pdf \
    --refPointLabel=TSS \
    --samplesLabel WT C9 \
    --regionsLabel "Genes downregulated in C9" \
    --legendLocation none \
    --averageTypeSummaryPlot mean \
    --colorMap=Blues \
    --missingDataColor=1 \
    --dpi=1200 \
    --heatmapWidth=3 \
    --heatmapHeight=10





#open fragments

computeMatrix reference-point -S WT.open.fragments.sorted.bam.bw C9.open.fragments.sorted.bam.bw \
    --referencePoint TSS -R ./genome/c9_upreg_genes.knownGene.gtf \
    -p max/2 \
    -bs=10 \
    --beforeRegionStartLength 1000 \
    --afterRegionStartLength 1000 \
    -o heatmap_open_fragments_c9_upreg_genes_atac.mat.gz

plotHeatmap -m heatmap_open_fragments_c9_upreg_genes_atac.mat.gz \
    -out ./ATAC_Seq/plots/heatmap_open_fragments_c9_upreg_genes_atac.pdf \
    --refPointLabel=TSS \
    --samplesLabel WT C9 \
    --regionsLabel "Genes upregulated in C9" \
    --legendLocation none \
    --averageTypeSummaryPlot mean \
    --colorMap=Reds \
    --missingDataColor=1 \
    --dpi=1200 \
    --heatmapWidth=3 \
    --heatmapHeight=10