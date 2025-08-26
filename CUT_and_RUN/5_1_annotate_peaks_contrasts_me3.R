library(ChIPseeker)
library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(Cairo)
library(colorspace)

#Function to annotate peak files
annotate_peaks <- function(peak_file_name, condition) {
  annotated_peaks <- annotatePeak(peak_file_name[[condition]], tssRegion = c(-3000, 3000),
                                  TxDb = edb, annoDb = "org.Hs.eg.db", overlap = "all")
  annotated_peaks_df <- as.data.frame(annotated_peaks)
  annotated_peaks_df$annotation <- gsub("^Intron.*", "Intron", annotated_peaks_df$annotation) #This merges all Intron/Exon as one annotation, independent of intron/exon number
  annotated_peaks_df$annotation <- gsub("^Exon.*", "Exon", annotated_peaks_df$annotation)
  return(annotated_peaks_df)
}


#load annotations
edb <- EnsDb.Hsapiens.v86
seqlevelsStyle(edb) <- "UCSC"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

dir <- "./CUT_and_RUN/1_2_Snakemake_run_merged_replicates/seacr/contrasts/"

#load peak files
peak_files <- list(C5_C9_common_me3 = paste0(dir, "c5_c9_me3_overlap_peaks.bed"),
                   WT_unique_me3 = paste0(dir, "WT_me3_unique_peaks.bed"),
                   C5_unique_me3 = paste0(dir, "C5_me3_unique_peaks.bed"),
                   C9_unique_me3 = paste0(dir, "C9_me3_unique_peaks.bed"), 
                   WT_C5_C9_overlap_me3 = paste0(dir, "wt_c5_c9_me3_overlap_peaks.bed"))


#Annotate peaks for each condition and save files
annotated_peaks_WT_unique_me3 <- annotate_peaks(peak_files, "WT_unique_me3") 
annotated_peaks_C5_unique_me3 <- annotate_peaks(peak_files, "C5_unique_me3")
annotated_peaks_C9_unique_me3 <- annotate_peaks(peak_files, "C9_unique_me3")
annotated_peaks_C5_C9_common_me3 <- annotate_peaks(peak_files, "C5_C9_common_me3") 
annotated_peaks_WT_C5_C9_overlap_me3 <- annotate_peaks(peak_files, "WT_C5_C9_overlap_me3")



write.csv(annotated_peaks_WT_unique_me3, file = paste0(dir, "./annotated_peaks/WT_unique_me3.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_C5_unique_me3, file = paste0(dir, "./annotated_peaks/C5_unique_me3.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_C9_unique_me3, file = paste0(dir, "./annotated_peaks/C9_unique_me3.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_C5_C9_common_me3, file = paste0(dir, "./annotated_peaks/C5_C9_common_me3.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_WT_C5_C9_overlap_me3, file = paste0(dir, "./annotated_peaks/WT_C5_C9_overlap_me3.csv"), quote = FALSE, row.names = FALSE)


#Count number of peaks in each type of annotation for H3K27me3 peaks
WT_unique_me3_summary <- as.data.frame(table(annotated_peaks_WT_unique_me3$annotation))
colnames(WT_unique_me3_summary) <- c("Annotation", "WT_unique_me3")

C5_unique_me3_summary <- as.data.frame(table(annotated_peaks_C5_unique_me3$annotation))
colnames(C5_unique_me3_summary) <- c("Annotation", "C5_unique_me3")

C9_unique_me3_summary <- as.data.frame(table(annotated_peaks_C9_unique_me3$annotation))
colnames(C9_unique_me3_summary) <- c("Annotation", "C9_unique_me3")

C5_C9_common_me3_summary <- as.data.frame(table(annotated_peaks_C5_C9_common_me3$annotation))
colnames(C5_C9_common_me3_summary) <- c("Annotation", "C5_C9_common_me3")

WT_C5_C9_overlap_me3_summary <- as.data.frame(table(annotated_peaks_WT_C5_C9_overlap_me3$annotation))
colnames(WT_C5_C9_overlap_me3_summary) <- c("Annotation", "WT_C5_C9_overlap_me3")



summary_annotation_me3 <- left_join(WT_unique_me3_summary, C5_unique_me3_summary, by = "Annotation") %>%
  left_join(C9_unique_me3_summary, by = "Annotation") %>%
  left_join(C5_C9_common_me3_summary, by = "Annotation") %>%
  left_join(WT_C5_C9_overlap_me3_summary, by = "Annotation") %>%
  pivot_longer(cols = -Annotation, names_to = "Condition", values_to = "value") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_me3$Annotation <- factor(summary_annotation_me3$Annotation, 
                                            levels = c("3' UTR", "Exon", 
                                                       "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", 
                                                       "Distal Intergenic", "Intron"))

summary_annotation_me3$Condition <- factor(summary_annotation_me3$Condition, 
                                           levels = c("WT_C5_C9_overlap_me3", "C5_C9_common_me3", "C9_unique_me3", "C5_unique_me3", "WT_unique_me3"))

#plot
annotated_me3_peaks_plot <- summary_annotation_me3 %>%
  ggplot(aes(x = Condition, y = value, fill = Annotation)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_text(aes(label=scales::comma(value)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", 
                               lighten("#009E73", amount = 0.5), "#009E73", darken("#009E73", amount = 0.5),  
                               "#FEB95F", "#CC79A7"), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma, limits = c(0, 50000)) +
  labs(x = "Sample") +
  ylab("No. H3K27me3 peaks")


ggsave(filename = "./CUT_and_RUN/plots/annotated_contrast_peaks_me3.pdf", plot = annotated_me3_peaks_plot, 
       width = 36, height = 24, dpi = 800, units = "cm", device = cairo_pdf)


WT_unique_me3_summary_perc <- WT_unique_me3_summary %>%
  mutate(WT_unique_me3_percentage = WT_unique_me3/sum(WT_unique_me3_summary$WT_unique_me3)*100) %>%
  select(-WT_unique_me3) 

C5_unique_me3_summary_perc <- C5_unique_me3_summary %>%
  mutate(C5_unique_me3_percentage = C5_unique_me3/sum(C5_unique_me3_summary$C5_unique_me3)*100) %>%
  select(-C5_unique_me3) 

C9_unique_me3_summary_perc <- C9_unique_me3_summary %>%
  mutate(C9_unique_me3_percentage = C9_unique_me3/sum(C9_unique_me3_summary$C9_unique_me3)*100) %>%
  select(-C9_unique_me3) 

C5_C9_common_me3_summary_perc <- C5_C9_common_me3_summary %>%
  mutate(C5_C9_common_me3_percentage = C5_C9_common_me3/sum(C5_C9_common_me3_summary$C5_C9_common_me3)*100) %>%
  select(-C5_C9_common_me3) 

WT_C5_C9_overlap_me3_summary_perc <- WT_C5_C9_overlap_me3_summary %>%
  mutate(WT_C5_C9_overlap_me3_percentage = WT_C5_C9_overlap_me3/sum(WT_C5_C9_overlap_me3_summary$WT_C5_C9_overlap_me3)*100) %>%
  select(-WT_C5_C9_overlap_me3) 

summary_annotation_perc <- left_join(WT_unique_me3_summary_perc, C5_unique_me3_summary_perc, by = "Annotation") %>%
  left_join(C9_unique_me3_summary_perc, by = "Annotation") %>%
  left_join(C5_C9_common_me3_summary_perc, by = "Annotation") %>%
  left_join(WT_C5_C9_overlap_me3_summary_perc, by = "Annotation") %>%
  pivot_longer(cols = -c("Annotation"), names_to = "Condition", values_to = "percentage") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_perc$Annotation <- factor(summary_annotation_perc$Annotation, 
                                             levels = c("3' UTR", "Exon",
                                                        "Distal Intergenic", "Intron",
                                                        "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)"))

summary_annotation_perc$Condition <- factor(summary_annotation_perc$Condition,
                                            levels = c("WT_C5_C9_overlap_me3_percentage",
                                                       "C5_C9_common_me3_percentage",
                                                       "C9_unique_me3_percentage",
                                                       "C5_unique_me3_percentage",
                                                       "WT_unique_me3_percentage"))

annotated_peaks_plot_perc <- summary_annotation_perc %>%
  ggplot(aes(x = Condition, y = percentage, fill = Annotation)) +
  geom_bar(stat = "identity", position = "stack") + 
  #geom_text(aes(label=scales::comma(value)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00",  
                               "#FEB95F", "#CC79A7", 
                               lighten("#009E73", amount = 0.5), "#009E73", darken("#009E73", amount = 0.5)), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  #scale_x_discrete(labels=c("WT_percentage" = "WT", "C9_percentage" = "C9")) +
  labs(x = "Sample") +
  ylab("H3K27me3 peaks (%)") + 
  theme(panel.grid = element_blank()) 

ggsave(filename = "./CUT_and_RUN/plots/annotated_constrast_me3_peaks_percentage.pdf", 
       plot = annotated_peaks_plot_perc, 
       width = 38, height = 13, dpi = 800, units = "cm", device = cairo_pdf)



peak_width_plot <- annotated_peaks_wt_me3 %>%
  mutate(Condition = "WT") %>%
  rbind(mutate(annotated_peaks_c5_me3, Condition = "C5")) %>%
  rbind(mutate(annotated_peaks_c9_me3, Condition = "C9")) %>%
  ggplot(aes(x = reorder(Condition, width), y = log10(width), colour = Condition)) +
  geom_boxplot(linewidth = 1.2, outlier.shape = 21, outlier.stroke = 0.4, outlier.fill = "white") +
  scale_colour_manual(values = c("#E69F00", "#56B4E9", "#000000")) +
  theme_classic(base_size = 26) +
  theme(legend.position = "none") +
  annotate("text", label = paste0(median(annotated_peaks_wt_me3$width), "bp"), x = 1, y = 3.12, size = 5.5) +
  annotate("text", label = paste0(median(annotated_peaks_c5_me3$width), "bp"), x = 2, y = 3.32, size = 5.5) +
  annotate("text", label = paste0(median(annotated_peaks_c9_me3$width), "bp"), x = 3, y = 3.6, size = 5.5) + 
  xlab("Sample") +
  ylab("H3K27me3 peaks width\n(log10bp)")

ggsave(filename = "./plots/peak_width.pdf", 
       plot = peak_width_plot, 
       width = 16, height = 14, dpi = 800, units = "cm", device = cairo_pdf)