library(tidyverse)
library(ggpubr)

#this takes annotated peaks outputs from script 5_0_annotate_peaks_me3_peak_width.R
annotated_peaks_wt_me3_genes <- annotated_peaks_WT_unique_me3 %>%
  dplyr::filter(SYMBOL != "NA")

expression_wt_vs_c5 <- read.csv("./RNA_Seq/results_files/normalised_expression_aml2_wt_vs_c9.csv") %>%
  dplyr::filter(geneID %in% annotated_peaks_wt_me3_genes$SYMBOL)

expression_wt_vs_c9 <- read.csv("./RNA_Seq/results_files/normalised_expression_aml2_wt_vs_c9.csv") %>%
  dplyr::filter(geneID %in% annotated_peaks_wt_me3_genes$SYMBOL)

corr_wt_vs_c9_wt_me3_peaks <- annotated_peaks_wt_me3_genes %>%
  left_join(expression_wt_vs_c9, by = c("SYMBOL" = "geneID")) %>%
  pivot_longer(c("aml2.wt.AVG", "aml2.c9.AVG"), names_to = "Condition", 
                      values_to = "RNA_expression") %>%
  dplyr::filter(RNA_expression != "NA") %>%
  dplyr::filter(annotation!= "Distal Intergenic") %>%
  dplyr::select("SYMBOL", "Condition", "RNA_expression") %>%
  unique()


plot_regions_losing_me3_c9 <- corr_wt_vs_c9_wt_me3_peaks %>%
  ggplot(aes(x = reorder(Condition, RNA_expression), y = RNA_expression, fill = Condition, colour = Condition)) +
  geom_boxplot(outlier.shape = NA, size = 0.7, alpha = 0.1) +
  geom_point(aes(fill=Condition, group = SYMBOL), size=1.8, shape=21, position = position_dodge(0.6), alpha = 0.3) +
  scale_fill_manual(values = c("#56B4E9", "black")) +
  scale_colour_manual(values = c("#56B4E9", "black")) +
  #geom_line(aes(group = SYMBOL), colour = "grey50", alpha = 0.3, position = position_dodge(0.6)) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none", plot.title = element_text(size = 15)) +
  stat_compare_means(paired = TRUE, label = "p.format", method = "wilcox.test", 
                     label.x.npc = "centre", size = 5) +
  # scale_y_continuous(limits = c(-2.5, 13), breaks = seq(-2, 12, by = 2)) +
  scale_x_discrete(labels = c("WT\nMe3", "C9\nNo Me3")) +
  labs(x = "Condition", y = expression(atop("RNA expression", log[2]~"(tpm+1)")), 
       title = "Regions losing me3 in C9")

  
ggsave(plot_regions_losing_me3_c9, file = "./CUT_and_RUN/plots/regions_losing_me3_c9.pdf",
       width = 12, height = 18, dpi = 800, units = "cm", device = cairo_pdf)
  




corr_wt_vs_c5_wt_me3_peaks <- annotated_peaks_wt_me3_genes %>%
  left_join(expression_wt_vs_c5, by = c("SYMBOL" = "geneID")) %>%
  pivot_longer(c("aml2.wt.AVG", "aml2.c5.AVG"), names_to = "Condition", 
               values_to = "RNA_expression") %>%
  dplyr::filter(RNA_expression != "NA") %>%
  dplyr::filter(annotation!= "Distal Intergenic") %>%
  dplyr::select("SYMBOL", "Condition", "RNA_expression") %>%
  unique()


plot_regions_losing_me3_c5 <- corr_wt_vs_c5_wt_me3_peaks %>%
  ggplot(aes(x = reorder(Condition, RNA_expression), y = RNA_expression, fill = Condition, colour = Condition)) +
  geom_boxplot(outlier.shape = NA, size = 0.7, alpha = 0.1) +
  geom_point(aes(fill=Condition, group = SYMBOL), size=1.8, shape=21, position = position_dodge(0.6), alpha = 0.3) +
  scale_fill_manual(values = c("#E69F00", "black")) +
  scale_colour_manual(values = c("#E69F00", "black")) +
  #geom_line(aes(group = SYMBOL), colour = "grey50", alpha = 0.3, position = position_dodge(0.6)) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none", plot.title = element_text(size = 15)) +
  stat_compare_means(paired = TRUE, label = "p.format", method = "wilcox.test", 
                     label.x.npc = "centre", size = 5) +
  # scale_y_continuous(limits = c(-2.5, 13), breaks = seq(-2, 12, by = 2)) +
  scale_x_discrete(labels = c("WT\nMe3", "C5\nNo Me3")) +
  labs(x = "Condition", y = expression(atop("RNA expression", log[2]~"(tpm+1)")), 
       title = "Regions losing me3 in C5")

ggsave(plot_regions_losing_me3_c5, file = "./CUT_and_RUN/plots/regions_losing_me3_c5.pdf",
       width = 12, height = 18, dpi = 800, units = "cm", device = cairo_pdf)

