library(biomaRt)
library(tidyverse)

#The outputs of this script are necessary for ATAC-Seq heatmaps of upregulated genes promoters

wt_vs_c5_diff_expr <- read.csv("./RNA_Seq/results_files/aml2_wt_v_c5_12k_genes.csv") %>%
  dplyr::filter(abs(logFC) >= 0.5) %>%
  dplyr::filter(adj.P.Val  < 0.1)
c5_down <- wt_vs_c5_diff_expr %>%
  dplyr::filter(logFC < -0.5) %>%
  dplyr::pull(geneID)
c5_up <- wt_vs_c5_diff_expr %>%
  dplyr::filter(logFC > 0.5) %>%
  dplyr::pull(geneID)

wt_vs_c9_diff_expr <- read.csv("./RNA_Seq/results_files/aml2_wt_v_c9_12k_genes.csv") %>%
  dplyr::filter(abs(logFC) > 0.5) %>%
  dplyr::filter(adj.P.Val  < 0.1)
c9_down <- wt_vs_c9_diff_expr %>%
  dplyr::filter(logFC < -0.5) %>%
  dplyr::pull(geneID)
c9_up <- wt_vs_c9_diff_expr %>%
  dplyr::filter(logFC > 0.5) %>%
  dplyr::pull(geneID)


mart <- useEnsembl("ensembl","hsapiens_gene_ensembl")

upreg_transcripts_c5 <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id", "transcript_is_canonical"),
                           filters = c("hgnc_symbol"),
                           values = list(c5_up),
                           mart = mart) %>%
  dplyr::filter(transcript_is_canonical == 1) %>%
  dplyr::pull("ensembl_transcript_id")

upreg_transcripts_c9 <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id", "transcript_is_canonical"),
                              filters = c("hgnc_symbol"),
                              values = list(c9_up),
                              mart = mart) %>%
  dplyr::filter(transcript_is_canonical == 1) %>%
  dplyr::pull("ensembl_transcript_id")

write.table(unique(upreg_transcripts_c5), col.names = F, row.names = F, quote = F, sep = "\t",
            file = "./RNA_Seq/results_files/transcript_ids_upreg_genes_c5.txt")

write.table(unique(upreg_transcripts_c9), col.names = F, row.names = F, quote = F, sep = "\t",
            file = "./RNA_Seq/results_files/transcript_ids_upreg_genes_c9.txt")


downreg_transcripts_c5 <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id", "transcript_is_canonical"),
                             filters = c("hgnc_symbol"),
                             values = list(c5_down),
                             mart = mart) %>%
  dplyr::filter(transcript_is_canonical == 1) %>%
  dplyr::pull("ensembl_transcript_id")

downreg_transcripts_c9 <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id", "transcript_is_canonical"),
                                filters = c("hgnc_symbol"),
                                values = list(c9_down),
                                mart = mart) %>%
  dplyr::filter(transcript_is_canonical == 1) %>%
  dplyr::pull("ensembl_transcript_id")

write.table(unique(downreg_transcripts_c5), col.names = F, row.names = F, quote = F, sep = "\t",
            file = "./RNA_Seq/results_files/transcript_ids_downreg_genes_c5.txt")
write.table(unique(downreg_transcripts_c9), col.names = F, row.names = F, quote = F, sep = "\t",
            file = "./RNA_Seq/results_files/transcript_ids_downreg_genes_c9.txt")