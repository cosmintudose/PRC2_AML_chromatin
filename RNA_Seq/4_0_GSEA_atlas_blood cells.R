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
library(ggrepel)

#downloaded from http://scrna.sklehabc.com/
atlas_blood_cells <- read.csv("./publicly_available_data/signatures_atlas_human_blood_cells.txt", sep = "\t") %>%
  dplyr::select(RNA_Cluster, Gene) %>%
  rename("gs_name" = "RNA_Cluster", "gene_symbol" = "Gene")


expression.aml2.df <- read.csv("./RNA_Seq/results_files/normalised_expression_aml2_wt_vs_clones.csv") %>%
 dplyr::select(geneID, LogFC, adj.P.Val)

data.aml2.gsea <- expression.aml2.df$LogFC
names(data.aml2.gsea) <- as.character(expression.aml2.df$geneID)
data.aml2.gsea <- sort(data.aml2.gsea, decreasing = TRUE)

# GSEA function from clusterProfiler
set.seed(1234)
GSEA.aml2.res <- GSEA(data.aml2.gsea, TERM2GENE=atlas_blood_cells, 
                        verbose=FALSE, seed = TRUE, pvalueCutoff = 0.05, minGSSize = 20)
GSEA.aml2.df.clones <- as_tibble(GSEA.aml2.res@result)

write.csv(GSEA.aml2.df.clones, file = "./RNA_Seq/results_files/GSEA_atlas_blood_cells_oci_aml2_clones_vs_wt.csv", row.names = FALSE, quote = FALSE)

# NES graph for each signature
gseaplot2(GSEA.aml2.res, 
          geneSetID = c(4, 3, 2), #can choose multiple signatures to overlay in this plot
          color = c("#2C467A", "#BE684D", "#5448C8"), 
          base_size = 18, #pvalue_table = TRUE,
          rel_heights = c(1.8, 0.6, 0.6))


#GSEA barplot
plot_gsea_abc <- GSEA.aml2.df.clones %>% 
  dplyr::filter(p.adjust < 0.05) %>%
  ggplot(aes(x = reorder(ID, NES), y = NES, fill = factor(sign(NES)))) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("#2C467A", "#BE684D")) +
  coord_flip() +
  theme_classic(base_size = 22) + 
  theme(legend.position = "none") +
  labs(title = "GSEA cell WT vs clones", x = element_blank(), y = "NES", legend = NA, 
       caption = "FDR < 0.05") +
  annotate(geom = "text", x = 8, y = -0.25, size = 5.5, fontface = "bold",
           label = "Enriched in\nPRC2-depleted", colour = "#BE684D", angle = 90) +
  annotate(geom = "text", x = 3, y = 0.25, size = 5.5, fontface = "bold",
           label = "Enriched in\nPRC2-WT", colour = "#2C467A", angle = 270)


#save GSEA barplot
ggsave(filename = "./RNA_Seq/plots/GSEA_ABC_clones_vs_wt.pdf", device = cairo_pdf, plot = plot_gsea_abc, 
       width = 10, height = 6, dpi = 1000)



#pull out core enchiment genes from GSEA
blood_cell_atlas_genes <- GSEA.aml2.df.clones %>%
  dplyr::filter(ID == "hMDP/cMoP" | ID == "Classical monocyte" | ID == "Non-classical monocyte") %>%
  dplyr::select(ID, core_enrichment) %>%
  separate_rows(core_enrichment, sep = "/")


#load differential exppression results and merge with core enrichment genes from GSEA
to_vplot <- read.csv("./RNA_Seq/results_files/aml2_wt_v_clones_12k_genes.csv") %>%
  left_join(blood_cell_atlas_genes, by = c("geneID" = "core_enrichment")) %>%
  mutate(ID = replace_na(ID, "none"))


#volcano plot differential expression analysis with GSEA core enrichment genes labelled
vplot <- to_vplot %>% 
  ggplot(aes(y=-log10(adj.P.Val), x = logFC, colour = ID, alpha = ID)) +
  geom_point(data = to_vplot %>% 
               filter(ID == "none"), size=2) +
  geom_hline(yintercept = -log10(0.1), linetype="longdash", colour="grey", size=0.5) +
  geom_vline(xintercept = 0.5, linetype="longdash", colour="grey", size=0.5) +
  geom_vline(xintercept = -0.5, linetype="longdash", colour="grey", size=0.5) +
  geom_point(data = to_vplot %>% 
               filter(ID != "none"), size=2) + #adding another layer so datapoints of interest are plotted on top
  geom_label_repel(data = to_vplot %>% 
                     filter(ID != "none" & abs(logFC) > 0.5 & adj.P.Val < 0.1),
                   aes(label = geneID, x = logFC, y = -log10(adj.P.Val)), force = 10,
                   hjust = "outward", box.padding = 0.1, max.overlaps = 8, min.segment.length = 0.1,
                   show.legend = FALSE, nudge_y = 3, alpha = 0.7, label.size = 0.1, direction = "y", size = 6) +
  scale_colour_manual(values = c("#2C467A", "#BE684D", "#5448C8", "grey40")) + 
  scale_alpha_manual(values = c(1, 1, 1, 0.5)) +
  xlim(c(-10, 10)) +
  theme_classic(base_size = 30) +
  labs(y = expression(-log[10]("FDR")))

#save vplot
ggsave(plot = vplot, file = "./RNA_Seq/plots/vplot_ABC_genes.pdf", device = cairo_pdf,
       dpi = 1000, height = 6, width = 12)