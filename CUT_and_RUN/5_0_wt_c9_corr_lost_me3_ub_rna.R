library(tidyverse)
library(ggpubr)

#this takes annotated peaks outputs from script 4_annotate_peaks.R
annotated_peaks_wt_me3_ub_genes <- annotated_peaks_wt_me3_ub %>%
  dplyr::filter(SYMBOL != "NA")

expression_wt_vs_c9 <- read.csv("./RNA_Seq/results_files/normalised_expression_aml2_wt_vs_c9.csv") %>%
  dplyr::filter(geneID %in% annotated_peaks_wt_me3_ub_genes$SYMBOL)

corr_wt_vs_c9_wt_me3_ub_peaks <- annotated_peaks_wt_me3_ub_genes %>%
  left_join(expression_wt_vs_c9, by = c("SYMBOL" = "geneID")) %>%
  pivot_longer(c("aml2.wt.AVG", "aml2.c9.AVG"), names_to = "Condition", 
                      values_to = "RNA_expression") %>%
  dplyr::filter(RNA_expression != "NA")


plot_regions_losing_me3_ub_c9 <- corr_wt_vs_c9_wt_me3_ub_peaks %>%
  dplyr::filter(annotation!= "Distal Intergenic") %>%
  ggplot(aes(x = reorder(Condition, RNA_expression), y = RNA_expression, fill = Condition, colour = Condition)) +
  geom_boxplot(outlier.shape = NA, size = 1, alpha = 0.1) +
  geom_point(aes(fill=Condition, group = SYMBOL), size=1.8, shape=21, position = position_dodge(0.6), alpha = 0.3) +
  scale_fill_manual(values = c("#56B4E9", "black")) +
  scale_colour_manual(values = c("#56B4E9", "black")) +
  #geom_line(aes(group = SYMBOL), colour = "grey50", alpha = 0.3, position = position_dodge(0.6)) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none", plot.title = element_text(size = 15)) +
  stat_compare_means(paired = TRUE, label = "p.format", label.x.npc = "centre", size = 5) +
  # scale_y_continuous(limits = c(-2.5, 13), breaks = seq(-2, 12, by = 2)) +
  scale_x_discrete(labels = c("WT\nMe3+Ub", "C9\nNo Me3/Ub")) +
  labs(x = "Condition", y = "RNA expression\nlog2(tpm+1)", title = "Regions losing me3&ub in C9") #+
  facet_wrap(~annotation, scales = "free", nrow = 2)

  
ggsave(plot_regions_losing_me3_ub_c9, file = "./CUT_and_RUN/plots/regions_losing_me3_ub_c9.pdf",
       width = 11, height = 18, dpi = 800, units = "cm", device = cairo_pdf)
  