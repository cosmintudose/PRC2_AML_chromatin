

For reproducibility of the R scripts an environment file is provided in the form of renv.lock 


/envs_configs directory contains .yml files of [conda](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) environments for different tasks as following:
  * hic_and_cnr_env.yml: CUT_and_RUN processing and Hi_C processing
  * deeptools.yml: heatmaps on CUT&RUN and ATAC-Seq data
  * homer.yml: Homer motif enrichment on ATAC-Seq peaks
  * nucleoatac.yml: nucleosome positioning inference using ATAC-Seq peaks
  * coolpuppy.yml: pileups on Hi-C maps
  * mustache.yml: call peaks in Hi-C data
  * pygenometracks.yml: create figures of Hi-C, ATAC-Seq and CUT&RUN tracks
       

## Overview figures from manuscript mapped to scripts

| Script                                                   | Figures                        | Brief description                                        |
|:---------------------------------------------------------|:-------------------------------|:---------------------------------------------------------|
| RNA_Seq/2_QC_and_PCA.R                                   | Fig. S1B                       | PCA RNA-Seq OCI-AML2 |
| RNA_Seq/3_0_WT_vs_clones_diff_expr.R                     | Fig. S1C                       | GSEA Hallmarks clones vs WT  |
| RNA_Seq/3_1_WT_vs_clones_sep_diff_expr.R                 | Fig. 1C                        | Venn diagrams overlaps C5 and C9 vs WT |
| RNA_Seq/4_0_GSEA_atlas_blood cells.R                     | Fig. 1D, E, S1D                | GSEA atlas blood cells on OCI-AML2 RNA-Seq |
| RNA_Seq/4_1_LIN28B_FL_GSEA.R                             | Fig. 5C                        | LIN28B transcriptional signature GSEA on OCI-AML2 RNA-Seq |
| RNA_Seq/5_0_boxplots_select_genes.R                      | Fig. 5D, S4B                   | Boxlots RNA-Seq per gene |
| RNA_Seq/6_LIN28B_expression_CCLE.R                       | Fig. S4C                       | LIN28B expression in the CCLE |
| CUT_and_RUN/3_1_pca_figures.R                            | Fig. S2A, B                    | PCA H3K27me3 and H3K27ac C&R |
| CUT_and_RUN/4_0_upset_overlaps_me3.R                     | Fig. 2A                        | Overlaps of called H3K27me3 CUT&RUN peaks |
| CUT_and_RUN/4_1_upset_overlaps_ac.R                      | Fig. S2F                       | Overlaps of called H3K27ac CUT&RUN peaks |
| CUT_and_RUN/5_0_annotate_peaks_me3_peak_width.R          | Fig. 2F                        | Peak width boxplot  |
| CUT_and_RUN/5_1_annotate_peaks_contrasts_me3.R           | Fig. 2B                        | Annotated me3 peaks  |
| CUT_and_RUN/5_2_annotate_peaks_contrasts_ac.R            | Fig. S2G                       | Annotated ac peaks  |
| CUT_and_RUN/6_0_wt_c9_corr_lost_me3_rna.R                | Fig. 2D                        | Boxplot correlation lost Me3+Ub with RNA |
| CUT_and_RUN/6_1_wt_c9_corr_gained_me3_rna.R              | Fig. 2E                        | Boxplot correlation gained Me3+Ub with RNA |
| CUT_and_RUN/7_h3k27ac_quantification.R                   | Fig. S2I                       | Venn diagram me3 ac RNA |
| CUT_and_RUN/8_heatmaps.sh                                | Fig. 2C, 3B, S2C, D, E         | CUT&RUN signal heatmaps and profile plots |
| ATAC_Seq/1_atac_pca_macs2_peaks.R                        | Fig. S3A                       | ATAC-Seq replicates PCA |
| ATAC_Seq/4_upset_atac_called_peaks.R                     | Fig. 3A, S3B                   | ATAC-Seq summaries across WT C5 and C9 |
| ATAC_Seq/5_overlaps_inter_intra_rna_atac.R               | Fig. S3E, F, G, H, I           | ATAC-Seq and RNA-Seq overlaps |
| ATAC_Seq/6_annotated_peaks_contrasts_upset.R             | Fig. 3C                        | Annotate ATAC peaks to genomic regions |
| ATAC_Seq/7_1_nucs_positioning.R                          | Fig. 3E, F, G, S3J             | Nucleosome positioning analysis |
| ATAC_Seq/8_heatmaps_atac_diff_expr_genes                 | Fig. S3C, D                    | Heatmaps of NFRs at TSS of up- and down-regulated genes |
| Hi_C/2_coolpuppy_pileups/Snakefile                       | Fig. 4B, C, D, S4D             | Pileups of loops/H3K27me3/H3K27ac regions |
| Hi_C/3_upset_loops.R                                     | Fig. 4A                        | Overlaps of called loops |
| pygenometracks_figs/make_pygenometracks_figs.sh          | Fig. 5A, B, E, S2H, S4A        | All tracks plots for specific genomic regions |









