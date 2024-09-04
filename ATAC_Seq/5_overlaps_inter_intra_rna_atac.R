library(VennDiagram)
library(WebGestaltR)
library(UpSetR)
library(tidyverse)
library(Cairo)

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

c5_open_promoters <- C5_minus_WT_annotated %>%
  dplyr::filter(annotation == "Promoter (<=1kb)") %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()

c9_open_promoters <- C9_minus_WT_annotated %>%
  dplyr::filter(annotation == "Promoter (<=1kb)") %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()

c5_closed_promoters <- WT_minus_C5_annotated %>%
  dplyr::filter(annotation == "Promoter (<=1kb)") %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()
#add to change thresholds # | annotation == "Promoter (1-2kb)" | annotation == "Promoter (2-3kb)" 

c9_closed_promoters <- WT_minus_C9_annotated %>%
  dplyr::filter(annotation == "Promoter (<=1kb)") %>%
  dplyr::pull("SYMBOL") %>%
  na.omit() %>%
  unique()

svg(file="./ATAC_Seq/plots/Upset_c5_RNA_ATAC.svg", width = 9, height = 6)

UpSetR::upset(fromList(list("C5 RNA \u2191" = c5_up, "C5 RNA \u2193" = c5_down, 
                    "C5 promoter accessibilty \u2191" = c5_open_promoters, "C5 promoter accessibilty \u2193" = c5_closed_promoters)),
      order.by = "freq", text.scale = 2, point.size = 3,
      queries = list(list(query = intersects, 
                          params = list("C5 RNA \u2191", "C5 promoter accessibilty \u2191"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C5 RNA \u2193", "C5 promoter accessibilty \u2193"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C5 RNA \u2193"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C5 promoter accessibilty \u2193"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C5 RNA \u2191"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C5 promoter accessibilty \u2191"), 
                          color = "#E68D74", active = T)),
      sets.bar.color = c("#E68D74", "#6679AB", "#E68D74", "#6679AB"))


dev.off()


svg(file="./ATAC_Seq/plots/Upset_c9_RNA_ATAC.svg", width = 9, height = 6)

UpSetR::upset(fromList(list("C9 RNA \u2191" = c9_up, "C9 RNA \u2193" = c9_down, 
                            "C9 promoter accessibilty \u2191" = c9_open_promoters, "C9 promoter accessibilty \u2193" = c9_closed_promoters)),
              order.by = "freq", text.scale = 2, point.size = 3,
              queries = list(list(query = intersects, 
                                  params = list("C9 RNA \u2191", "C9 promoter accessibilty \u2191"), 
                                  color = "#E68D74", active = T),
                             list(query = intersects, 
                                  params = list("C9 RNA \u2193", "C9 promoter accessibilty \u2193"), 
                                  color = "#6679AB", active = T),
                             list(query = intersects, 
                                  params = list("C9 RNA \u2193"), 
                                  color = "#6679AB", active = T),
                             list(query = intersects, 
                                  params = list("C9 promoter accessibilty \u2193"), 
                                  color = "#6679AB", active = T),
                             list(query = intersects, 
                                  params = list("C9 RNA \u2191"), 
                                  color = "#E68D74", active = T),
                             list(query = intersects, 
                                  params = list("C9 promoter accessibilty \u2191"), 
                                  color = "#E68D74", active = T)),
              sets.bar.color = c("#E68D74", "#6679AB", "#E68D74", "#6679AB"))

dev.off()


venn.diagram(list("C5 RNA \u2191\nC5 promoter accessibilty \u2191" = intersect(c5_up, c5_open_promoters), 
                  "C9 RNA \u2191\nC9 promoter accessibilty \u2191" = intersect(c9_up, c9_open_promoters)), 
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5), cat.pos = c(340, 20), fill = c("#ffccd5", "#c9184a"),
             "./upset_plots/atac_open_rna_up_c5_c9_overlaps.pdf", disable.logging = TRUE)

genes_to_plot <- intersect(intersect(c5_up, c5_open_promoters), intersect(c9_up, c9_open_promoters))


fisher.test(matrix(c(length(union(c5_up, c5_open_promoters))-length(union(intersect(c5_up, c5_open_promoters),intersect(c9_up, c9_open_promoters))), 
                     length(setdiff(intersect(c5_up, c5_open_promoters), intersect(c9_up, c9_open_promoters))), 
                     length(setdiff(intersect(c9_up, c9_open_promoters), intersect(c5_up, c5_open_promoters))), 
                     length(intersect(intersect(c5_up, c5_open_promoters), intersect(c9_up, c9_open_promoters)))), nrow = 2), alternative = "greater")


fisher.test(matrix(c(length(union(c5_down, c5_closed_promoters))-length(union(intersect(c5_down, c5_closed_promoters),intersect(c9_down, c9_closed_promoters))), 
                     length(setdiff(intersect(c5_down, c5_closed_promoters), intersect(c9_down, c9_closed_promoters))), 
                     length(setdiff(intersect(c9_down, c9_closed_promoters), intersect(c5_down, c5_closed_promoters))), 
                     length(intersect(intersect(c5_down, c5_closed_promoters), intersect(c9_down, c9_closed_promoters)))), nrow = 2), alternative = "greater")

library(pheatmap)
library(RColorBrewer)

normalised_expression_aml2 <- read.csv("/home/cosmin/rna_seq_clones/RNASeq_EZH2_KO/aml2_ezh2_ko/results_files/normalised_expression_aml2_wt_vs_clones.csv") %>%
  dplyr::filter(geneID %in% genes_to_plot)

rownames(normalised_expression_aml2) <- normalised_expression_aml2$geneID

normalised_expression_aml2 <- normalised_expression_aml2 %>%
  select(-c("geneID", "aml2.wt.AVG", "aml2.clones.AVG", "LogFC"))

#heatmap of selected genes
pheatmap(normalised_expression_aml2, scale = "row", cluster_cols = F, cluster_rows = F, angle_col = 0, fontsize = 20,
         color=colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))


venn.diagram(list("C5 RNA \u2193\nC5 promoter accessibilty \u2193" = intersect(c5_down, c5_closed_promoters), 
                  "C9 RNA \u2193\nC9 promoter accessibilty \u2193" = intersect(c9_down, c9_closed_promoters)), 
             lwd = 0, cex = 3, cat.cex = 1.6, print.mode = "raw", margin = 0.1,
             alpha = c(0.5, 0.5), cat.pos = c(340, 20), fill = c("#a9d6e5", "#014f86"),
             "./upset_plots/atac_closed_rna_down_c5_c9_overlaps.pdf", disable.logging = TRUE)

genes_to_plot <- intersect(intersect(c5_down, c5_closed_promoters), intersect(c9_down, c9_closed_promoters))

normalised_expression_aml2 <- read.csv("/home/cosmin/rna_seq_clones/RNASeq_EZH2_KO/aml2_ezh2_ko/results_files/normalised_expression_aml2_wt_vs_clones.csv") %>%
  dplyr::filter(geneID %in% genes_to_plot)

rownames(normalised_expression_aml2) <- normalised_expression_aml2$geneID

normalised_expression_aml2 <- normalised_expression_aml2 %>%
  select(-c("geneID", "aml2.wt.AVG", "aml2.clones.AVG", "LogFC"))

#heatmap of selected genes
pheatmap(normalised_expression_aml2, scale = "row", cluster_cols = F, cluster_rows = F, angle_col = 0, fontsize = 20,
         color=colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))



list_sets <- c("C5 RNA up", "C5 RNA down", 
                  "C5 promoter open", "C5 promoter closed",
                  "C9 RNA up", "C9 RNA down", 
                  "C9 promoter open", "C9 promoter closed")

list_sets <- factor(list_sets, levels=list_sets)


upset(fromList(list("C5 RNA up" = c5_up, "C5 RNA down" = c5_down, 
                    "C5 promoter open" = c5_open_promoters, "C5 promoter closed" = c5_closed_promoters,
                    "C9 RNA up" = c9_up, "C9 RNA down" = c9_down, 
                    "C9 promoter open" = c9_open_promoters, "C9 promoter closed" = c9_closed_promoters)),
      order.by = "freq", text.scale = 1.5, nsets = 8, nintersects = 55, sets = list_sets,
      queries = list(list(query = intersects, 
                          params = list("C5 RNA up", "C5 promoter open"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C5 RNA down", "C5 promoter closed"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C5 RNA down"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C5 promoter closed"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C5 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C5 promoter open"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C9 RNA up", "C9 promoter open"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C9 RNA down", "C9 promoter closed"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C9 RNA down"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C9 promoter closed"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C9 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C9 promoter open"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C5 promoter open", "C9 promoter open"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C5 promoter closed", "C9 promoter closed"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C5 RNA down", "C9 RNA down"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C5 promoter open", "C9 promoter open", "C5 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C9 promoter open", "C5 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C5 promoter open", "C9 promoter open", "C9 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C5 promoter closed", "C9 promoter closed", "C5 RNA down"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C5 promoter closed", "C9 promoter closed", "C5 RNA down", "C9 RNA down"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C5 promoter open", "C9 promoter open", "C5 RNA up", "C9 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C5 RNA up", "C9 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C5 promoter closed", "C9 RNA down"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C5 promoter closed", "C5 RNA down", "C9 RNA down"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C9 promoter closed", "C5 RNA down"), 
                          color = "#6679AB", active = T),
                     list(query = intersects, 
                          params = list("C9 promoter open", "C5 RNA up", "C9 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C5 promoter open", "C5 RNA up", "C9 RNA up"), 
                          color = "#E68D74", active = T),
                     list(query = intersects, 
                          params = list("C5 promoter open", "C9 RNA up"), 
                          color = "#E68D74", active = T)),
      sets.bar.color = c("#E68D74", "#E68D74", "#6679AB", "#6679AB",
                         "#E68D74", "#6679AB", "#6679AB", "#E68D74"))

venn.diagram(list("C5 promoter\nopen" = c5_open_promoters, "C9 promoter\nopen" = c9_open_promoters), 
             lwd = 0, cex = 2, cat.cex = 3, print.mode = "raw",
             alpha = c(0.5, 0.5), cat.pos = c(190, 152), fill = c("#ffccd5", "#c9184a"),
             "./ATAC_Seq/plots/open_atac_overlaps.svg", disable.logging = TRUE)

venn.diagram(list("C5 promoter\nclosed" = c5_closed_promoters, "C9 promoter\nclosed" = c9_closed_promoters), 
             lwd = 0, cex = 2, cat.cex = 3, print.mode = "raw",
             alpha = c(0.5, 0.5), cat.pos = c(200, 170), fill = c("#a9d6e5", "#014f86"),
             "./ATAC_Seq/plots/closed_atac_overlaps.svg", disable.logging = TRUE)

venn.diagram(list("C5 promoter\nopen" = c5_open_promoters, "C5 RNA up" = c5_up), 
             lwd = 0, cex = 2, cat.cex = 3, print.mode = "raw",
             alpha = c(0.5, 0.5), cat.pos = c(190, 152), fill = c("#ffccd5", "#c9184a"),
             "./ATAC_Seq/plots/atac_open_rna_up_c5_overlaps.svg", disable.logging = TRUE)

venn.diagram(list("C5 promoter\nclosed" = c5_closed_promoters, "C5 RNA down" = c5_down), 
             lwd = 0, cex = 2, cat.cex = 3, print.mode = "raw",
             alpha = c(0.5, 0.5), cat.pos = c(195, 175), fill = c("#a9d6e5", "#014f86"),
             "./ATAC_Seq/plots/atac_closed_rna_down_c5_overlaps.svg", disable.logging = TRUE)


venn.diagram(list("C9 promoter\nopen" = c9_open_promoters, "C9 RNA up" = c9_up), 
             lwd = 0, cex = 2, cat.cex = 3, print.mode = "raw",
             alpha = c(0.5, 0.5), cat.pos = c(190, 152), fill = c("#ffccd5", "#c9184a"),
             "./ATAC_Seq/plots/atac_open_rna_up_c9_overlaps.svg", disable.logging = TRUE)

venn.diagram(list("C9 promoter\nclosed" = c9_closed_promoters, "C9 RNA down" = c9_down), 
             lwd = 0, cex = 2, cat.cex = 3, print.mode = "raw",
             alpha = c(0.5, 0.5), cat.pos = c(195, 175), fill = c("#a9d6e5", "#014f86"),
             "./ATAC_Seq/plots/atac_closed_rna_down_c9_overlaps.svg", disable.logging = TRUE)