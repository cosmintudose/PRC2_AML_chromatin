library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot) 
library(cowplot)

#These data were kindly provided by Rebecca Ling and Anindita Roy, Department of Paediatrics, University of Oxford, Oxford, United Kingdom
#script also needs outputs of 3.WT_vs_clones_diff_expr.R to run 
lin28b_fl_kds <- read.csv("./LIN28B_KO_DEGs/FL_LIN28BKD_allDEG_17MArch24.csv") %>%
  dplyr::filter(log2FoldChange > 1 | log2FoldChange < 1) %>%
  mutate(gs_name = case_when(log2FoldChange > 1 ~ "Genes upregulated upon LIN28B KD in FL CD34+",
                             log2FoldChange < -1 ~ "Genes downregulated upon LIN28B KD in FL CD34+")) %>%
  dplyr::select(c("gs_name", "gene")) #preparing data to test as gene sets

# columns corresponding to gene symbols and LogFC for one pairwise comparison for the enrichment analysis
aml2.df.sub <- dplyr::select(aml2.df, geneID, LogFC) #select genes and their logFCs from AML2s RNA-seq 
aml2.gsea <- aml2.df.sub$LogFC
names(aml2.gsea) <- as.character(aml2.df.sub$geneID)
aml2.gsea <- sort(aml2.gsea, decreasing = TRUE) #genes are sorted from most upregulated to most downregulated

# GSEA function from clusterProfiler
set.seed(1234)
GSEA.aml2.res <- GSEA(aml2.gsea, TERM2GENE=lin28b_fl_kds, 
                        verbose=FALSE, seed = TRUE, nPermSimple = 10000, eps = 0, pvalueCutoff = 0.2, pAdjustMethod = "bonferroni")
GSEA.aml2.df <- as_tibble(GSEA.aml2.res@result)

# NES graph for any each signature
gseaplot2(GSEA.aml2.res, 
          geneSetID = c(1,2), #can choose multiple signatures to overlay in this plot
          pvalue_table = TRUE, 
          base_size = 20, 
          color = c("#805D93", "#DF9A57"),
          title = "OCI-AML2 EZH2+/- vs EZH2+/+\nGene sets determined from LIN28B KD in FL CD34+") #can turn off 


anno <- GSEA.aml2.res[1, c("NES", "pvalue", "p.adjust")]
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
p1 <- gsearank(GSEA.aml2.res, 1, title = paste0("Gene set - ", GSEA.aml2.res[1, "Description"])) +
  annotate("text", 400, GSEA.aml2.res[1, "enrichmentScore"] * .75, label = lab, hjust=0, vjust=0.7) + 
  annotate("text", 10200, GSEA.aml2.res[1, "enrichmentScore"]  * .75, label = "Downregulated in EZH2+/-", colour = "#2C467A", vjust=-12) + 
  annotate("text", 2000, GSEA.aml2.res[1, "enrichmentScore"]  * .75, label = "Upregulated in EZH2+/-", colour = "#BE684D", vjust=-8.5)


anno <- GSEA.aml2.res[2, c("NES", "pvalue", "p.adjust")]
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
p2 <- gsearank(GSEA.aml2.res, 2, title = paste0("Gene set - ", GSEA.aml2.res[2, "Description"])) +
  annotate("text", 400, GSEA.aml2.res[2, "enrichmentScore"] * .75, label = lab, hjust=0, vjust=3.2) + 
  annotate("text", 10200, GSEA.aml2.res[2, "enrichmentScore"]  * .75, label = "Downregulated in EZH2+/-", colour = "#2C467A", vjust=6.8) + 
  annotate("text", 2200, GSEA.aml2.res[2, "enrichmentScore"]  * .75, label = "Upregulated in EZH2+/-", colour = "#BE684D", vjust=9.3) + 
  geom_segment(aes(x = which(rownames(as.data.frame(aml2.gsea)) == "CDK6"), 
                   xend = which(rownames(as.data.frame(aml2.gsea)) == "CDK6"),
                   y = 0, yend = 0.265), colour = "#F51AA4", linetype = "dashed", size = 0.4) +
  annotate("text", 200, 0.32, label = "CDK6", colour = "#BD4089")


pdf(file="./plots/running_enrichment_scores_with_CDK6.pdf", width = 6, height = 6)
plot_grid(p1, p2, nrow=2)
dev.off()

