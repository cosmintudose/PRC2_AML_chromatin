library(tidyverse)
library(hicVennDiagram)
library(GenomicInteractions)
library(ComplexUpset)


dir <- "./ATAC_Seq/bwa/merged_replicate/hmmratac_peaks"

fs <- dir("./ATAC_Seq/bwa/merged_replicate/hmmratac_peaks/", 
          pattern = "_ATAC_peaks.bed", full.names = TRUE)
fs <- fs[!grepl("overlap", fs)]

venn <- vennCount(fs, maxgap = 50, FUN = max)


upset_themes_fix <- lapply(ComplexUpset::upset_themes, function(.ele){
  lapply(.ele, function(.e){
    do.call(theme, .e[names(.e) %in% names(formals(theme))])
  })
})

upsetPlot(venn,
          themes = list(default=theme_bw()))

combinations <- venn$combinations
expInput <- venn$counts

plotdata <- combinations[rep(rownames(combinations), expInput), ] %>%
  as.data.frame()

names(plotdata) <- c("OCI AML2 C5", "OCI AML2 C9", "OCI AML2 WT")

plotdata <- plotdata %>%
  dplyr::select("OCI AML2 C9", "OCI AML2 C5", "OCI AML2 WT")

 

upset_atac <- ComplexUpset::upset(
  data=as.data.frame(plotdata),
  intersect=colnames(plotdata),
  themes = upset_modify_themes(list('overall_sizes'=theme(axis.text.x=element_text(angle=90), 
                                                          panel.grid = element_blank()), 
                                    'Intersection size'=theme(panel.grid = element_blank(),
                                                              text=element_text(size=14)),
                                    'intersections_matrix'=theme(panel.grid = element_blank(),
                                                                 text=element_text(size=14)))),
  min_size = 50,
  sort_sets = FALSE,
  queries = list(
    upset_query(set="OCI AML2 WT",fill="black"),
    upset_query(set = "OCI AML2 C5", fill = "#E69F00"),
    upset_query(set = "OCI AML2 C9", fill = "#56B4E9"), 
    upset_query(intersect=c('OCI AML2 C5', 'OCI AML2 C9'), fill='#E68D74', color = "#E68D74"),
    upset_query(intersect=c('OCI AML2 C5'), fill='#E68D74', color = "#E68D74"),
    upset_query(intersect=c('OCI AML2 C9'), fill='#E68D74', color = "#E68D74"),
    upset_query(intersect=c('OCI AML2 WT', 'OCI AML2 C9'), fill='#6679AB', color = "#6679AB"),
    upset_query(intersect=c('OCI AML2 WT', 'OCI AML2 C5'), fill='#6679AB', color = "#6679AB"),
    upset_query(intersect=c('OCI AML2 WT'), fill='#6679AB', color = "#6679AB")),
  height_ratio=0.7,
  width_ratio=0.2,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE))) + 
  labs(x = "Called accessible peaks")

ggsave(filename = "./ATAC_Seq/plots/upset_peaks_overlaps.pdf", plot = upset_atac, 
       width = 18, height = 10, dpi = 800, units = "cm", device = cairo_pdf)

dir.create("./homer_input")
####printing files for homer analysis
OCI_AML2_WT_peaks_only <- venn$overlapList[["001"]][["WT_ATAC_peaks.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

write.table(OCI_AML2_WT_peaks_only, file = "./homer_input_after_upset/OCI_AML2_WT_peaks_only.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


OCI_AML2_C5_C9_overlap_peaks <- venn$overlapList[["110"]][["C5_ATAC_peaks.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

write.table(OCI_AML2_C5_C9_overlap_peaks, file = "./homer_input_/OCI_AML2_C5_C9_overlap_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

OCI_AML2_C5_peaks <- venn$overlapList[["100"]][["C5_ATAC_peaks.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end")) %>%
  rbind(OCI_AML2_C5_C9_overlap_peaks)

write.table(OCI_AML2_C5_peaks, file = "./homer_input/OCI_AML2_C5_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

OCI_AML2_C9_peaks <- venn$overlapList[["010"]][["C9_ATAC_peaks.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end")) %>%
  rbind(OCI_AML2_C5_C9_overlap_peaks)

write.table(OCI_AML2_C9_peaks, file = "./homer_input/OCI_AML2_C9_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

#For simplicity, these numbers are manually added in the plot, taken from the upSet plot created earlier - make sure they match
plot <- data.frame("C5 gained\naccessible regions" = venn@counts[["110"]]+venn@counts[["100"]], 
                   "C9 gained\naccessible regions" = venn@counts[["010"]]+venn@counts[["110"]], 
                   "C5 lost\naccessible regions" = venn@counts[["011"]]+venn@counts[["001"]], 
                   "C9 lost\naccessible regions" = venn@counts[["001"]] + venn@counts[["101"]], 
                   check.names = FALSE) %>%
  pivot_longer(cols = everything(), names_to = "Description", values_to = "Count")

plot$Type <- c("Gained vs WT", "Gained vs WT", "Lost vs WT", "Lost vs WT")

plot$Description <- as.factor(plot$Description) %>%
  fct_rev()


gained_lost_peaks_plot <- plot %>%
  ggplot(aes(x = Description, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label=scales::comma(Count)), hjust = -0.1, size = 7) +
  theme_classic(base_size = 26) + 
  theme(axis.text.y=element_text(hjust = 0.5)) +
  coord_flip() + 
  scale_y_continuous(labels = scales::comma, limits = c(0, 40000)) +
  scale_fill_manual(values = c("#BE684D", "#2C467A")) + 
  labs(x = "", y = "No. accessible regions") 

ggsave(filename = "./ATAC_Seq/plots/summary_gained_lost_peaks.pdf", plot = gained_lost_peaks_plot, 
       width = 29, height = 12, dpi = 800, units = "cm", device = cairo_pdf)
