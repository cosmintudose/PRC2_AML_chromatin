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

dir <- "./CUT_and_RUN/1_2_Snakemake_run_merged_replicates/seacr/filtered/"

#load peak files
peak_files <- list(WT_me3 = paste0(dir, "WT_me3_filtered.stringent.bed"), 
                   C5_me3 = paste0(dir, "C5_me3_filtered.stringent.bed"),
                   C9_me3 = paste0(dir, "C9_me3_filtered.stringent.bed"),
                   WT_ac = paste0(dir, "WT_ac_filtered.stringent.bed"), 
                   C5_ac = paste0(dir, "C5_ac_filtered.stringent.bed"),
                   C9_ac = paste0(dir, "C9_ac_filtered.stringent.bed"))

#Annotate peaks for each condition and save files
annotated_peaks_wt_me3 <- annotate_peaks(peak_files, "WT_me3") 
annotated_peaks_c5_me3 <- annotate_peaks(peak_files, "C5_me3")
annotated_peaks_c9_me3 <- annotate_peaks(peak_files, "C9_me3")
annotated_peaks_wt_ac <- annotate_peaks(peak_files, "WT_ac") 
annotated_peaks_c5_ac <- annotate_peaks(peak_files, "C5_ac")
annotated_peaks_c9_ac <- annotate_peaks(peak_files, "C9_ac")

dir.create("./CUT_and_RUN/1_2_Snakemake_run_merged_replicates/seacr/filtered/annotated_peaks")

write.csv(annotated_peaks_wt_me3, file = paste0(dir, "./annotated_peaks/WT_h3k27me3.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_c5_me3, file = paste0(dir, "./annotated_peaks/C5_h3k27me3.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_c9_me3, file = paste0(dir, "./annotated_peaks/C9_h3k27me3.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_wt_ac, file = paste0(dir, "./annotated_peaks/WT_h3k27ac.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_c5_ac, file = paste0(dir, "./annotated_peaks/C5_h3k27ac.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_c9_ac, file = paste0(dir, "./annotated_peaks/C9_h3k27ac.csv"), quote = FALSE, row.names = FALSE)

#Count number of peaks in each type of annotation for H3K27me3 peaks
wt_me3_summary <- as.data.frame(table(annotated_peaks_wt_me3$annotation))
colnames(wt_me3_summary) <- c("Annotation", "WT")

c5_me3_summary <- as.data.frame(table(annotated_peaks_c5_me3$annotation))
colnames(c5_me3_summary) <- c("Annotation", "C5")

c9_me3_summary <- as.data.frame(table(annotated_peaks_c9_me3$annotation))
colnames(c9_me3_summary) <- c("Annotation", "C9")


summary_annotation_me3 <- left_join(wt_me3_summary, c5_me3_summary, by = "Annotation") %>%
  left_join(c9_me3_summary, by = "Annotation") %>%
  pivot_longer(cols = -Annotation, names_to = "Condition", values_to = "value") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_me3$Annotation <- factor(summary_annotation_me3$Annotation, 
                                        levels = c("3' UTR", "Exon", 
                                                   "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", 
                                                   "Distal Intergenic", "Intron"))

summary_annotation_me3$Condition <- factor(summary_annotation_me3$Condition, 
                                           levels = c("C9", "C5", "WT"))

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
  scale_y_continuous(labels = scales::comma, limits = c(0, 65000)) +
  labs(x = "Sample") +
  ylab("No. H3K27me3 peaks")



#Count number of peaks in each type of annotation for H2AK119Ub peaks
wt_ub_summary <- as.data.frame(table(annotated_peaks_wt_ub$annotation))
colnames(wt_ub_summary) <- c("Annotation", "WT")

c5_ub_summary <- as.data.frame(table(annotated_peaks_c5_ub$annotation))
colnames(c5_ub_summary) <- c("Annotation", "C5")

c9_ub_summary <- as.data.frame(table(annotated_peaks_c9_ub$annotation))
colnames(c9_ub_summary) <- c("Annotation", "C9")


#Count number of peaks in each type of annotation for H3K27ac peaks
wt_ac_summary <- as.data.frame(table(annotated_peaks_wt_ac$annotation))
colnames(wt_ac_summary) <- c("Annotation", "WT")

c5_ac_summary <- as.data.frame(table(annotated_peaks_c5_ac$annotation))
colnames(c5_ac_summary) <- c("Annotation", "C5")

c9_ac_summary <- as.data.frame(table(annotated_peaks_c9_ac$annotation))
colnames(c9_ac_summary) <- c("Annotation", "C9")


summary_annotation_ac <- left_join(wt_ac_summary, c5_ac_summary, by = "Annotation") %>%
  left_join(c9_ac_summary, by = "Annotation") %>%
  pivot_longer(cols = -Annotation, names_to = "Condition", values_to = "value") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_ac$Annotation <- factor(summary_annotation_ac$Annotation, 
                                           levels = c("3' UTR", "Exon", 
                                                      "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", 
                                                      "Distal Intergenic", "Intron"))

summary_annotation_ac$Condition <- factor(summary_annotation_ac$Condition, 
                                          levels = c("C9", "C5", "WT"))


#plot
annotated_ac_peaks_plot <- summary_annotation_ac %>%
  ggplot(aes(x = Condition, y = value, fill = Annotation)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_text(aes(label=scales::comma(value)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", 
                               lighten("#009E73", amount = 0.5), "#009E73", darken("#009E73", amount = 0.5),  
                               "#FEB95F", "#CC79A7"), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma, limits = c(0, 190000)) +
  labs(x = "Sample") +
  ylab("No. H3K27ac peaks")


###Peak width
to_plot_peak_width <- annotated_peaks_wt_me3 %>%
  mutate(Condition = "WT") %>%
  rbind(mutate(annotated_peaks_c5_me3, Condition = "C5")) %>%
  rbind(mutate(annotated_peaks_c9_me3, Condition = "C9"))
  
  
to_plot_peak_width$Condition <- factor(to_plot_peak_width$Condition, levels = c("WT", "C5", "C9"))

peak_width_plot <- to_plot_peak_width %>%
  ggplot(aes(x = Condition, y = log10(width), colour = Condition)) +
  geom_boxplot(linewidth = 1.2, outlier.shape = 21, outlier.stroke = 0.4, outlier.fill = "white") +
  scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9")) +
  theme_classic(base_size = 26) +
  theme(legend.position = "none") +
  annotate("text", label = paste0(median(annotated_peaks_wt_me3$width), "bp"), x = 1, y = 3.28, size = 5.5) +
  annotate("text", label = paste0(median(annotated_peaks_c5_me3$width), "bp"), x = 2, y = 3.32, size = 5.5) +
  annotate("text", label = paste0(median(annotated_peaks_c9_me3$width), "bp"), x = 3, y = 3.6, size = 5.5) + 
  xlab("Sample") +
  ylab(expression(atop("H3K27me3 peaks width", (log[10]~"bp"))))


ggsave(filename = "./CUT_and_RUN/plots/peak_width.pdf", 
       plot = peak_width_plot, 
       width = 14, height = 15, dpi = 800, units = "cm", device = cairo_pdf)