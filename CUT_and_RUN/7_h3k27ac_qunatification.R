library(Rsubread)
library(tidyverse)
library(edgeR)
library(reshape2)
library(VennDiagram)

#this script uses output from 5_0_annotate_peaks_peak_width.R

annotated_peaks_wt_ac_saf <- annotated_peaks_wt_ac %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()

annotated_peaks_c5_ac_saf <- annotated_peaks_c5_ac %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()

annotated_peaks_c9_ac_saf <- annotated_peaks_c9_ac %>%
  dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "seqnames", "start", "end", "strand")) %>%
  dplyr::rename(c("GeneID" = "SYMBOL", 
                  "Chr" = "seqnames",
                  "Start" = "start", 
                  "End" = "end", 
                  "Strand" = "strand")) %>%
  na.omit()


WT_counts <- featureCounts("./CUT_and_RUN/1_2_Snakemake_run_merged_replicates/deduped/WT_ac/WT_ac_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_wt_ac_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  rename(c("GeneID" = "Var1", "WT_counts" = "value")) %>%
  na.omit()


C5_counts <- featureCounts("./CUT_and_RUN/1_2_Snakemake_run_merged_replicates/deduped/C5_ac/C5_ac_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_c5_ac_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  rename(c("GeneID" = "Var1", "C5_counts" = "value")) %>%
  na.omit()

C9_counts <- featureCounts("./CUT_and_RUN/1_2_Snakemake_run_merged_replicates/deduped/C9_ac/C9_ac_deduped.bam", # count reads
                           annot.inbuilt = "hg38",
                           annot.ext = annotated_peaks_c9_ac_saf,
                           isGTFAnnotationFile = F,
                           isPairedEnd = T,
                           nthreads = 8) %>%
  melt() %>%
  dplyr::select(c("Var1", "value")) %>%
  rename(c("GeneID" = "Var1", "C9_counts" = "value")) %>%
  na.omit()


C9_counts$GeneID <- as.character(C9_counts$GeneID)
C5_counts$GeneID <- as.character(C5_counts$GeneID)
WT_counts$GeneID <- as.character(WT_counts$GeneID)


gene_lengths <- annotated_peaks_wt_ac %>%
  rbind(annotated_peaks_c9_ac) %>%
  #dplyr::filter(transcriptBiotype == "protein_coding") %>%
  dplyr::select(c("SYMBOL", "geneLength")) %>%
  slice_max(geneLength, by = SYMBOL) %>%
  na.omit() %>%
  unique() %>%
  rename("GeneID" = "SYMBOL")

gene_lengths$geneLength <- gene_lengths$geneLength + 3000

counts_matrix <- full_join(WT_counts, C9_counts, by = "GeneID") %>%
  full_join(C5_counts, by = "GeneID") %>%
  dplyr::left_join(gene_lengths, by = "GeneID")

counts_matrix$WT_counts <- as.numeric(counts_matrix$WT_counts)
counts_matrix$C9_counts <- as.numeric(counts_matrix$C9_counts)
counts_matrix$C5_counts <- as.numeric(counts_matrix$C5_counts)

counts_matrix <- counts_matrix %>%
  mutate_if(is.numeric, ~replace_na(., 0))

rownames(counts_matrix) <- counts_matrix$GeneID

counts_matrix <- counts_matrix %>%
  dplyr::select(-GeneID) %>%
  as.matrix()


counts_list <- DGEList(counts=counts_matrix[,c("WT_counts", "C5_counts", "C9_counts")], genes=data.frame(Length=counts_matrix[,"geneLength"]))

counts_list <- calcNormFactors(counts_list)

RPKM <- log2(rpkm(counts_list)+1) %>%
  as.data.frame()

RPKM$logFC_C5_WT <- RPKM$C5_counts - RPKM$WT_counts
RPKM$logFC_C9_WT <- RPKM$C9_counts - RPKM$WT_counts

RPKM$geneID <- rownames(RPKM)

write.csv(RPKM, file = "./ac_rpkm.csv", quote = FALSE, row.names = FALSE)

wt_vs_clones_diff_expr <- read.csv("./RNA_Seq/results_files/aml2_wt_v_clones_12k_genes.csv") %>%
  mutate(mRNA_status = case_when(logFC<0 ~ "Downregulated", 
                                 logFC>0 ~ "Upregulated"))

ac_and_rna <- left_join(RPKM, wt_vs_clones_diff_expr, by = "geneID") %>%
  na.omit()

ac_and_rna %>%
  dplyr::filter(logFC > 0.5) %>%
  dplyr::filter(adj.P.Val < 0.1) %>%
  dplyr::filter(logFC_C5_WT > 0.5) %>%
  dplyr::filter(logFC_C9_WT > 0.5) %>%
  pull(geneID)

rna_up_genes <- ac_and_rna %>%
  dplyr::filter(logFC > 0.5) %>%
  dplyr::filter(adj.P.Val < 0.1) %>%
  pull(geneID)

ac_up_c5 <- ac_and_rna %>%
  dplyr::filter(logFC_C5_WT > 0.5) %>%
  pull(geneID)

ac_up_c9 <- ac_and_rna %>%
  dplyr::filter(logFC_C9_WT > 0.5) %>%
  pull(geneID)



venn.diagram(list("C5 H3K27ac \u2191" = ac_up_c5, 
                  "C9 H3K27ac \u2191" = ac_up_c9,
                  "C5&C9 RNA \u2191" = rna_up_genes),
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5, 0.5), fill = c("#ffccd5", "#c9184a", "#E16036"),
             "./CUT_and_RUN/plots/rna_ac_overlaps.pdf", disable.logging = TRUE)

ac_and_rna %>%
  dplyr::filter(logFC < -0.5) %>%
  dplyr::filter(adj.P.Val < 0.1) %>%
  dplyr::filter(logFC_C5_WT < 0) %>%
  dplyr::filter(logFC_C9_WT < 0) %>%
  pull(geneID)


