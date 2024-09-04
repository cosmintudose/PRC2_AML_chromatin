library(ChIPseeker)
library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(Cairo)
library(colorspace)

edb <- EnsDb.Hsapiens.v86
seqlevelsStyle(edb) <- "UCSC"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

dir <- "./ATAC_Seq/bwa/merged_replicate/hmmratac_peaks/atac_peak_calling/"

dir.create("./ATAC_Seq/bwa/merged_replicate/hmmratac_peaks/atac_peak_calling/annotated_peaks")

peak_files <- list(WT = paste0(dir, "WT_ATAC_peaks.bed"), 
                   C5 = paste0(dir, "C5_ATAC_peaks.bed"), 
                   C9 = paste0(dir, "C9_ATAC_peaks.bed"))


annotated_peaks_wt <- annotatePeak(peak_files[["WT"]], tssRegion = c(-3000, 3000),
                                   TxDb = edb, annoDb = "org.Hs.eg.db", overlap = "all")
annotated_peaks_wt_df <- as.data.frame(annotated_peaks_wt)
write.csv(annotated_peaks_wt_df, file = paste0(dir, "./annotated_peaks/WT_ATAC_peaks.csv"), row.names = FALSE, quote = FALSE)

annotated_peaks_wt_df$annotation <- gsub("^Intron.*", "Intron", annotated_peaks_wt_df$annotation)
annotated_peaks_wt_df$annotation <- gsub("^Exon.*", "Exon", annotated_peaks_wt_df$annotation)


annotated_peaks_c5 <- annotatePeak(peak_files[["C5"]], tssRegion = c(-3000, 3000),
                                   TxDb = edb, annoDb = "org.Hs.eg.db", overlap = "all")
annotated_peaks_c5_df <- as.data.frame(annotated_peaks_c5)
write.csv(annotated_peaks_c5_df, file = paste0(dir, "./annotated_peaks/C5_ATAC_peaks.csv"), row.names = FALSE, quote = FALSE)

annotated_peaks_c5_df$annotation <- gsub("^Intron.*", "Intron", annotated_peaks_c5_df$annotation)
annotated_peaks_c5_df$annotation <- gsub("^Exon.*", "Exon", annotated_peaks_c5_df$annotation)



annotated_peaks_c9 <- annotatePeak(peak_files[["C9"]], tssRegion = c(-3000, 3000),
                                   TxDb = edb, annoDb = "org.Hs.eg.db", overlap = "all")
annotated_peaks_c9_df <- as.data.frame(annotated_peaks_c9)
write.csv(annotated_peaks_c9_df, file = paste0(dir, "./annotated_peaks/C9_ATAC_peaks.csv"), row.names = FALSE, quote = FALSE)


annotated_peaks_c9_df$annotation <- gsub("^Intron.*", "Intron", annotated_peaks_c9_df$annotation)
annotated_peaks_c9_df$annotation <- gsub("^Exon.*", "Exon", annotated_peaks_c9_df$annotation)

wt_summary <- as.data.frame(table(annotated_peaks_wt_df$annotation))
colnames(wt_summary) <- c("Annotation", "WT")
c5_summary <- as.data.frame(table(annotated_peaks_c5_df$annotation))
colnames(c5_summary) <- c("Annotation", "C5")
c9_summary <- as.data.frame(table(annotated_peaks_c9_df$annotation))
colnames(c9_summary) <- c("Annotation", "C9")


summary_annotation <- left_join(wt_summary, c5_summary, by = "Annotation") %>%
  left_join(c9_summary, by = "Annotation") %>%
  pivot_longer(cols = -Annotation, names_to = "Condition", values_to = "value") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation$Annotation <- factor(summary_annotation$Annotation, 
                                               levels = c("3' UTR", "Exon", 
                                                          "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", 
                                                          "Distal Intergenic", "Intron"))

legend_ordfactor(summary_annotation$Annotation, 
       levels = c("3' UTR", "Exon", 
                  "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", 
                  "Distal Intergenic", "Intron"))


annotated_peaks_plot <- summary_annotation %>%
  ggplot(aes(x = reorder(Condition, value, decreasing = TRUE), y = value, fill = Annotation)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_text(aes(label=scales::comma(value)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", 
                               lighten("#009E73", amount = 0.5), "#009E73", darken("#009E73", amount = 0.5),  
                               "#FEB95F", "#CC79A7"), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma, limits = c(0, 45000)) +
  labs(x = "Sample") +
  ylab("No. accessible regions")

ggsave(filename = "./ATAC_Seq/plots/summary_peaks_annotated.pdf", plot = annotated_peaks_plot, 
       width = 28, height = 18, dpi = 800, units = "cm", device = cairo_pdf)

wt_summary_perc <- wt_summary %>%
  mutate(WT_percentage = WT/sum(wt_summary$WT)*100) %>%
  select(-WT)

c5_summary_perc <- c5_summary %>%
  mutate(C5_percentage = C5/sum(c5_summary$C5)*100) %>%
  select(-C5)

c9_summary_perc <- c9_summary %>%
  mutate(C9_percentage = C9/sum(c9_summary$C9)*100) %>%
  select(-C9)


summary_annotation_perc <- left_join(wt_summary_perc, c5_summary_perc, by = "Annotation") %>%
  left_join(c9_summary_perc, by = "Annotation") %>%
  pivot_longer(cols = -Annotation, names_to = "Condition", values_to = "percentage") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_perc$Annotation <- factor(summary_annotation_perc$Annotation, 
                                        levels = c("3' UTR", "Exon", 
                                                   "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", 
                                                   "Distal Intergenic", "Intron"))

annotated_peaks_plot <- summary_annotation_perc %>%
  ggplot(aes(x = reorder(Condition, percentage), y = percentage, fill = Annotation)) +
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#0072B2", "#D55E00", 
                               lighten("#009E73", amount = 0.5), "#009E73", darken("#009E73", amount = 0.5),  
                               "#FEB95F", "#CC79A7"), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  scale_x_discrete(labels=c("WT_percentage" = "WT", "C5_percentage" = "C5",
                            "C9_percentage" = "C9")) +
  labs(x = "Sample") +
  ylab("Accessible regions (%)")

ggsave(filename = "./ATAC_Seq/plots/summary_peaks_annotated_percentage.pdf", 
       plot = annotated_peaks_plot, 
       width = 27, height = 10, dpi = 800, units = "cm", device = cairo_pdf)