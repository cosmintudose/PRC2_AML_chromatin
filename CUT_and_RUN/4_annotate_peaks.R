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

dir <- "./CUT_and_RUN/1_Snakemake_run/seacr/"

#load peak files
peak_files <- list(WT_me3 = paste0(dir, "WT_h3k27me3.stringent.bed"), 
                   C9_me3 = paste0(dir, "C9_h3k27me3.stringent.bed"),
                   WT_ub = paste0(dir, "WT_h2ak119ub.stringent.bed"), 
                   C9_ub = paste0(dir, "C9_h2ak119ub.stringent.bed"))

#Annotate peaks for each condition and save files
annotated_peaks_wt_me3 <- annotate_peaks(peak_files, "WT_me3") 
annotated_peaks_c9_me3 <- annotate_peaks(peak_files, "C9_me3")
annotated_peaks_wt_ub <- annotate_peaks(peak_files, "WT_ub") 
annotated_peaks_c9_ub <- annotate_peaks(peak_files, "C9_ub")


create.dir("/CUT_and_RUN/1_Snakemake_run/seacr/annotated_peaks")

write.csv(annotated_peaks_wt_me3, file = paste0(dir, "./annotated_peaks/WT_h3k27me3.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_c9_me3, file = paste0(dir, "./annotated_peaks/C9_h3k27me3.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_wt_ub, file = paste0(dir, "./annotated_peaks/WT_h2ak119ub.csv"), quote = FALSE, row.names = FALSE)
write.csv(annotated_peaks_c9_ub, file = paste0(dir, "./annotated_peaks/C9_h2ak119ub.csv"), quote = FALSE, row.names = FALSE)

#Count number of peaks in each type of annotation for H3K27me3 peaks
wt_me3_summary <- as.data.frame(table(annotated_peaks_wt_me3$annotation))
colnames(wt_me3_summary) <- c("Annotation", "WT")

c9_me3_summary <- as.data.frame(table(annotated_peaks_c9_me3$annotation))
colnames(c9_me3_summary) <- c("Annotation", "C9")


summary_annotation_me3 <- left_join(wt_me3_summary, c9_me3_summary, by = "Annotation") %>%
  pivot_longer(cols = -Annotation, names_to = "Condition", values_to = "value") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_me3$Annotation <- factor(summary_annotation_me3$Annotation, 
                                        levels = c("3' UTR", "Exon", 
                                                   "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", 
                                                   "Distal Intergenic", "Intron"))

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
  scale_y_continuous(labels = scales::comma, limits = c(0, 18000)) +
  labs(x = "Sample") +
  ylab("No. H3K27me3 peaks")

dir.create("./plots")

ggsave(filename = "./CUT_and_RUN/plots/summary_me3_peaks_annotated.pdf", plot = annotated_me3_peaks_plot, 
       width = 26, height = 14, dpi = 800, units = "cm", device = cairo_pdf)


#Count number of peaks in each type of annotation for H2AK119Ub peaks
wt_ub_summary <- as.data.frame(table(annotated_peaks_wt_ub$annotation))
colnames(wt_ub_summary) <- c("Annotation", "WT")

c9_ub_summary <- as.data.frame(table(annotated_peaks_c9_ub$annotation))
colnames(c9_ub_summary) <- c("Annotation", "C9")


summary_annotation_ub <- left_join(wt_ub_summary, c9_ub_summary, by = "Annotation") %>%
  pivot_longer(cols = -Annotation, names_to = "Condition", values_to = "value") %>%
  dplyr::filter(Annotation != "5' UTR") %>%
  dplyr::filter(Annotation != "Downstream (<=300bp)")

summary_annotation_ub$Annotation <- factor(summary_annotation_ub$Annotation, 
                                           levels = c("3' UTR", "Exon", 
                                                      "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", 
                                                      "Distal Intergenic", "Intron"))

#plot
annotated_ub_peaks_plot <- summary_annotation_ub %>%
  ggplot(aes(x = Condition, y = value, fill = Annotation)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_text(aes(label=scales::comma(value)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", 
                               lighten("#009E73", amount = 0.5), "#009E73", darken("#009E73", amount = 0.5),  
                               "#FEB95F", "#CC79A7"), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma, limits = c(0, 10000)) +
  labs(x = "Sample") +
  ylab("No. H2AK119ub peaks")

ggsave(filename = "./CUT_and_RUN/plots/summary_ub_peaks_annotated.pdf", plot = annotated_ub_peaks_plot, 
       width = 26, height = 14, dpi = 800, units = "cm", device = cairo_pdf)


wt_me3_summary_perc <- wt_me3_summary %>%
  mutate(WT_percentage = WT/sum(wt_me3_summary$WT)*100,
         PTM = "H3K27me3") %>%
  select(-WT) 

c9_me3_summary_perc <- c9_me3_summary %>%
  mutate(C9_percentage = C9/sum(c9_me3_summary$C9)*100) %>%
  select(-C9)

wt_ub_summary_perc <- wt_ub_summary %>%
  mutate(WT_percentage = WT/sum(wt_ub_summary$WT)*100,
         PTM = "H2AK119ub") %>%
  select(-WT) 

c9_ub_summary_perc <- c9_ub_summary %>%
  mutate(C9_percentage = C9/sum(c9_ub_summary$C9)*100) %>%
  select(-C9)



summary(annotated_peaks_wt_df$width)
summary(annotated_peaks_c9_df$width)


peak_width_plot <- annotated_peaks_wt_df %>%
  mutate(Condition = "WT") %>%
  rbind(mutate(annotated_peaks_c9_df, Condition = "C9")) %>%
  ggplot(aes(x = reorder(Condition, width), y = log10(width), colour = Condition)) +
  geom_boxplot(linewidth = 1.2, outlier.shape = 21, outlier.stroke = 0.4, outlier.fill = "white") +
  scale_colour_manual(values = c("#56B4E9", "#000000")) +
  theme_classic(base_size = 26) +
  theme(legend.position = "none") +
  annotate("text", label = paste0(median(annotated_peaks_wt_df$width), "bp"), x = 1, y = 3.42, size = 5.5) +
  annotate("text", label = paste0(median(annotated_peaks_c9_df$width), "bp"), x = 2, y = 3.75, size = 5.5) + 
  xlab("Sample") +
  ylab("H3K27me3 peaks width\n(log10bp)")

ggsave(filename = "./plots/peak_width.pdf", 
       plot = peak_width_plot, 
       width = 12, height = 14, dpi = 800, units = "cm", device = cairo_pdf)

genes_covered_by_peaks_wt <- annotated_peaks_wt_df %>%
  dplyr::filter(SYMBOL != "NA") %>%
  dplyr::filter(annotation != "Distal Intergenic") %>%
  dplyr::filter(annotation != "Downstream (<=300bp)") %>%
  dplyr::filter(annotation != "5' UTR") %>%
  group_by(annotation) %>%
  summarise(genes_covered = n_distinct(SYMBOL)) %>%
  mutate(Condition = "WT")

genes_covered_by_peaks_c9 <- annotated_peaks_c9_df %>%
  dplyr::filter(SYMBOL != "NA") %>%
  dplyr::filter(annotation != "Distal Intergenic") %>%
  dplyr::filter(annotation != "Downstream (<=300bp)") %>%
  dplyr::filter(annotation != "5' UTR") %>%
  group_by(annotation) %>%
  summarise(genes_covered = n_distinct(SYMBOL)) %>%
  mutate(Condition = "C9")

summary_annotation_genes_covered <- rbind(genes_covered_by_peaks_wt, genes_covered_by_peaks_c9) %>%
  rename("Annotation" = "annotation")

summary_annotation_genes_covered$Annotation <- factor(summary_annotation_genes_covered$Annotation, 
                                        levels = c("3' UTR", "Exon", 
                                                   "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", 
                                                   "Intron"))

annotation_genes_covered_plot <- summary_annotation_genes_covered %>%
  ggplot(aes(x = Condition, y = genes_covered, fill = Annotation)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  geom_text(aes(label=scales::comma(genes_covered)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", 
                               lighten("#009E73", amount = 0.5), "#009E73", darken("#009E73", amount = 0.5),  
                               "#CC79A7"), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma, limits = c(0, 3900)) +
  labs(x = "Sample") +
  ylab("No. genes covered by H3K27me3")

ggsave(filename = "./CUT_and_RUN/plots/annotation_genes_covered.pdf", plot = annotation_genes_covered_plot, 
       width = 26, height = 13, dpi = 800, units = "cm", device = cairo_pdf)


genes_covered_by_peaks_wt_perc <- genes_covered_by_peaks_wt %>%
  mutate(percentage = genes_covered/sum(genes_covered_by_peaks_wt$genes_covered)*100) %>%
  select(-genes_covered)

genes_covered_by_peaks_c9_perc <- genes_covered_by_peaks_c9 %>%
  mutate(percentage = genes_covered/sum(genes_covered_by_peaks_c9$genes_covered)*100) %>%
  select(-genes_covered)


summary_annotation_perc <- rbind(genes_covered_by_peaks_wt_perc, genes_covered_by_peaks_c9_perc) %>%
  dplyr::filter(annotation != "5' UTR") %>%
  dplyr::filter(annotation != "Downstream (<=300bp)")

summary_annotation_perc$annotation <- factor(summary_annotation_perc$annotation, 
                                             levels = c("3' UTR", "Exon", 
                                                        "Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", 
                                                        "Distal Intergenic", "Intron"))

annotated_peaks_plot <- summary_annotation_perc %>%
  ggplot(aes(x = reorder(Condition, percentage), y = percentage, fill = annotation)) +
  geom_bar(stat = "identity", position = "stack") + 
  #geom_text(aes(label=scales::comma(value)), hjust = -0.1, size = 6.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", 
                               lighten("#009E73", amount = 0.5), "#009E73", darken("#009E73", amount = 0.5),  
                               "#FEB95F", "#CC79A7"), guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 26) +
  coord_flip() +
  scale_x_discrete(labels=c("WT_percentage" = "WT", "C9_percentage" = "C9")) +
  labs(x = "Sample") +
  ylab("H3K27me3 peaks (%)")
