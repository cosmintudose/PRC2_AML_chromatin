library(tidyverse)
library(ggrepel)
library(ggbeeswarm)
library(scales)

#Downloaded from depmap.org > Custom Downloads > selected LIN28B and AML cell lines only version 24Q2
lin28b_expr_aml <- read.csv("./publicly_available_data/Batch_corrected_Expression_Public_24Q2_subsetted.csv")

plot_lin28b_ccle <- lin28b_expr_aml %>%
  ggplot(aes(x = lineage_1, y = LIN28B, label = cell_line_display_name)) +
  geom_quasirandom(dodge.width = 0.1, size = 4, cex = 4, color = "black", fill = "#F3DFA2", shape = 21, width = 0.5,
                mapping = aes(factor(lineage_1), LIN28B)) +
  geom_label_repel(mapping = aes(factor(lineage_1), LIN28B,
                                 label = ifelse(cell_line_display_name %in% c("MOLM13", "OCIAML2"), cell_line_display_name, "")),
                  position = position_quasirandom(width = 0.5), box.padding = 1.2, size = 6) +
  theme_classic(base_size = 28) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  scale_size_binned() +
  labs(x = "AML cell lines", y = expression("LIN28B log"[2] *"(tpm+1)"))
 
ggsave(filename = "./plots/LIN28B_expression_CCLE_AML.pdf", plot = plot_lin28b_ccle, 
              width = 12, height = 22, dpi = 800, units = "cm")
