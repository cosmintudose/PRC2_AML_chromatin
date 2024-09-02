
/envs_configs directory contains .yml files of [conda](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) environments for different tasks as following:
  * hic_and_cnr_env.yml: CUT_and_RUN processing and Hi_C processing
  * deeptools.yml: heatmaps on CUT&RUN and ATAC-Seq data
  * homer.yml: Homer motif enrichment on ATAC-Seq peaks
  * nucleoatac.yml: nucleosome positioning inference using ATAC-Seq peaks
  * coolpuppy.yml: pileups on Hi-C maps
  * mustache.yml: call peaks in Hi-C data
  * pygenometracks.yml: create figures of Hi-C, ATAC-Seq and CUT&RUN tracks
       

## Overview figures from manuscript mapped to scripts

| Script                                                   | Figures                      | Brief description                                        |
|:---------------------------------------------------------|:-----------------------------|:---------------------------------------------------------|
| RNA_Seq/2_QC_and_PCA.R                                   | Fig. S2A                     | PCA RNA-Seq OCI-AML2 |
| RNA_Seq/3_0_WT_vs_clones_diff_expr.R                     | Fig. 2B                      | GSEA Hallmarks clones vs WT  |
| RNA_Seq/3_1_WT_vs_clones_sep_diff_expr.R                 | Fig. 2A                      | Venn diagrams overlaps C5 and C9 vs WT |
| RNA_Seq/4_0_GSEA_atlas_blood cells.R                     | Fig. 2C, S2B, E, F           | GSEA atlas blood cells on OCI-AML2 RNA-Seq |
| RNA_Seq/4_1_LIN28B_FL_GSEA.R                             | Fig. 6B                      | LIN28B transcriptional signature GSEA on OCI-AML2 RNA-Seq |
| RNA_Seq/5_0_boxplots_select_genes.R                      | Fig. 2E, 6C, S2E, F, S8C, D  | Boxlots RNA-Seq per gene |
| RNA_Seq/6_LIN28B_expression_CCLE.R                       | Fig. S7E                     | LIN28B expression in the CCLE |
| RNA_Seq/7_0_TARGET_cohort_analysis.R                     | Fig. S2C                     | Diff expr TARGET cohort & GSEA atlas blood cells |
| RNA_Seq/7_1_TARGET_cohort_analysis_KMT2Ar_restricted.R   | Fig. 2D, F, S2D, G, H        | Like 7_0, but restricted to KMT2Ar samples |
