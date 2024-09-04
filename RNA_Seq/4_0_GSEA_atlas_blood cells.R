#takes the output from 3.WT_vs_clones_diff_expr.R as input

library(tidyverse)
library(ensembldb) 
library(gplots)
library(RColorBrewer)
library(GSEABase) 
library(gprofiler2) 
library(clusterProfiler) 
library(enrichplot) 
library(colorspace)
library(Cairo)

#downloaded from http://scrna.sklehabc.com/
atlas_blood_cells <- read.csv("./publicly_available_data/signatures_atlas_human_blood_cells.txt", sep = "\t") %>%
  dplyr::select(RNA_Cluster, Gene) %>%
  rename("gs_name" = "RNA_Cluster", "gene_symbol" = "Gene")


expression.aml2.df <- read.csv("./RNA_Seq/results_files/normalised_expression_aml2_wt_vs_clones.csv") %>%
 dplyr::select(geneID, LogFC)

data.aml2.gsea <- expression.aml2.df$LogFC
names(data.aml2.gsea) <- as.character(expression.aml2.df$geneID)
data.aml2.gsea <- sort(data.aml2.gsea, decreasing = TRUE)

# GSEA function from clusterProfiler
set.seed(1234)
GSEA.aml2.res <- GSEA(data.aml2.gsea, TERM2GENE=atlas_blood_cells, 
                        verbose=FALSE, seed = TRUE, pvalueCutoff = 1, minGSSize = 1)
GSEA.aml2.df.clones <- as_tibble(GSEA.aml2.res@result)

write.csv(GSEA.aml2.df.clones, file = "./RNA_Seq/results_files/GSEA_atlas_blood_cells_oci_aml2_clones_vs_wt.csv", row.names = FALSE, quote = FALSE)

# NES graph for any each signature
gseaplot2(GSEA.aml2.res, 
          geneSetID = c(4, 3, 2), #can choose multiple signatures to overlay in this plot
          color = c("#2C467A", "#BE684D", "#5448C8"), 
          base_size = 18, #pvalue_table = TRUE,
          rel_heights = c(1.8, 0.6, 0.6))


plot_gsea_abc <- GSEA.aml2.df.clones %>% 
  dplyr::filter(p.adjust < 0.1) %>%
  ggplot(aes(x = reorder(ID, NES), y = NES, fill = factor(sign(NES)))) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("#2C467A", "#BE684D")) +
  coord_flip() +
  theme_classic(base_size = 22) + 
  theme(legend.position = "none") +
  labs(title = "GSEA cell WT vs clones", x = element_blank(), y = "NES", legend = NA, 
       caption = "FDR < 0.1") +
  annotate(geom = "text", x = 14.5, y = -0.25, size = 5.5, fontface = "bold",
           label = "Enriched in\nPRC2-depleted", colour = "#BE684D", angle = 90) +
  annotate(geom = "text", x = 6, y = 0.25, size = 5.5, fontface = "bold",
           label = "Enriched in\nPRC2-WT", colour = "#2C467A", angle = 270)


ggsave(filename = "./RNA_Seq/plots/GSEA_ABS_clones_vs_wt.pdf", device = cairo_pdf, plot = plot_gsea_abc, 
       width = 10, height = 8, dpi = 1000)