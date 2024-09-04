

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
| RNA_Seq/2_QC_and_PCA.R                                   | Fig. S2A                       | PCA RNA-Seq OCI-AML2 |
| RNA_Seq/3_0_WT_vs_clones_diff_expr.R                     | Fig. 2B                        | GSEA Hallmarks clones vs WT  |
| RNA_Seq/3_1_WT_vs_clones_sep_diff_expr.R                 | Fig. 2A                        | Venn diagrams overlaps C5 and C9 vs WT |
| RNA_Seq/4_0_GSEA_atlas_blood cells.R                     | Fig. 2C, S2B, E, F             | GSEA atlas blood cells on OCI-AML2 RNA-Seq |
| RNA_Seq/4_1_LIN28B_FL_GSEA.R                             | Fig. 6B                        | LIN28B transcriptional signature GSEA on OCI-AML2 RNA-Seq |
| RNA_Seq/5_0_boxplots_select_genes.R                      | Fig. 2E, 6C, S2E, F, S8C, D    | Boxlots RNA-Seq per gene |
| RNA_Seq/6_LIN28B_expression_CCLE.R                       | Fig. S7E                       | LIN28B expression in the CCLE |
| RNA_Seq/7_0_TARGET_cohort_analysis.R                     | Fig. S2C                       | Diff expr TARGET cohort & GSEA atlas blood cells |
| RNA_Seq/7_1_TARGET_cohort_analysis_KMT2Ar_restricted.R   | Fig. 2D, F, S2D, G, H          | Like 7_0, but restricted to KMT2Ar samples |
| CUT_and_RUN/3_0_upset_overlaps.R                         | Fig. 3A                        | Overlaps of called CUT&RUN peaks |
| CUT_and_RUN/4_0_annotate_peaks.R                         | Fig. 3B, C, F, S3B             | Annotated CUT&RUN peaks and peak width boxplot  |
| CUT_and_RUN/5_0_wt_c9_corr_lost_me3_ub_rna.R             | Fig. 3G                        | Boxplot correlation lost Me3+Ub with RNA |
| CUT_and_RUN/5_1_wt_c9_corr_gained_me3_ub_rna.R           | Fig. 3G                        | Boxplot correlation gained Me3+Ub with RNA |
| CUT_and_RUN/6_heatmaps.sh                                | Fig. 3D, E, S3A                | CUT&RUN signal heatmaps and profile plots |
| ATAC_Seq/1_atac_pca_macs2_peaks.R                        | Fig. S4A                       | ATAC-Seq replicates PCA |
| ATAC_Seq/4_upset_atac_called_peaks.R                     | Fig. 4A, S4B                   | ATAC-Seq summaries across WT C5 and C9 |
| ATAC_Seq/5_overlaps_inter_intra_rna_atac.R               | Fig. S5C, D, E, F, G           | ATAC-Seq and RNA-Seq overlaps |
| ATAC_Seq/6_annotate_peaks.R                              | Fig. 4B                        | Annotate ATAC peaks to genomic regions |
| ATAC_Seq/7_1_nucs_positioning.R                          | Fig. 4E, F, G, S6              | Nucleosome positioning analysis |
| ATAC_Seq/8_heatmaps_atac_diff_expr_genes                 | Fig. S5A, B                    | Heatmaps of NFRs at TSS of up- and down-regulated genes |
| Hi_C/2_coolpuppy_pileups/Snakefile                       | Fig. 5B, C, S7D                | Pileups of loops/H3K27me3 regions |
| Hi_C/3_upset_loops.R                                     | Fig. 5A                        | Overlaps of called loops |
| pygenometracks_figs/make_pygenometracks_figs.sh          | Fig. 6A, S3C, D, E, F, S8A, B  | All tracks plots for specific genomic regions |









