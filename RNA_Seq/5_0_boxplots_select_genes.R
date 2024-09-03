library(tidyverse)

#select any gene to plot boxplot WT, C5, and C9 - example given is with LIN28B and CDK6 (Figure 6C)
selected_genes <- c("LIN28B", "CDK6")

selected_genes_expression <- read.csv("./RNA_Seq/results_files/normalised_counts_aml2_wt_vs_clones.tsv", sep = "\t") %>%
  dplyr::filter(geneID %in% selected_genes)


boxplots <- selected_genes_expression %>%
  pivot_longer(-geneID, names_to = "Sample", values_to = "mRNA") %>%
  mutate(Condition = substr(Sample, 1, 2)) 

boxplots$Condition <- factor(boxplots$Condition, levels = c("WT", "C5", "C9"))
boxplots$geneID <- factor(boxplots$geneID, selected_genes)

plot_boxplots <- boxplots %>%
  ggplot(aes(x = Condition, y = mRNA, fill = Condition)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = c("black", "#E69F00", "#56B4E9")) +
  geom_jitter(width = 0.25, size = 2.8, fill = "black") +
  facet_wrap(~geneID, scales = "free") + #this line is not necessary if plotting just one gene
  theme_bw(base_size = 24) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(y = expression("log"[2] *"(tpm+1)"), legend = NA, title = "RNA-Seq")

ggsave("./RNA_Seq/plots/selected_genes_boxplots.pdf", 
       plot = p1, width = 16, height = 12, 
       dpi = 1000, units = "cm", device = cairo_pdf) 