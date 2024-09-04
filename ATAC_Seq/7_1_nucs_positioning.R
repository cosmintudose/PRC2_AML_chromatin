library(tidyverse)
library(ChIPseeker)
library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(Cairo)
library(ggpubr)
library(mixtools)
library(plotGMM)
library(ClusterR)


#create function for mixed modelling
mixed_model <- function(data, sample_selection) {
  
  data_selected <- data %>%
    dplyr::filter(sample == sample_selection)
  
  opt_gmm <- Optimal_Clusters_GMM(as.data.frame(data_selected$distanceToTSS), max_clusters = 8, 
                                  criterion = "AIC")
  
  set.seed(12)
  mixres3 = normalmixEM(data_selected$distanceToTSS, k=5, maxit = 10000)
  
  peaks <- round(mixres3[["mu"]], digits = 0)
  
  post.df3 <- as.data.frame(cbind(score = mixres3$x, mixres3$posterior))
  post.df3 <- post.df3 %>%
    mutate(Component = case_when(comp.1 > 0.5 ~ 1,
                                 comp.2 > 0.5 ~ 2,
                                 comp.4 > 0.5 ~ 4, 
                                 comp.5 > 0.5 ~ 5))
  
  
  threshold1 <- post.df3 %>% 
    dplyr::filter(Component == 1) %>%
    pull("score") %>%
    min()
  
  threshold2 <- post.df3 %>% 
    dplyr::filter(Component == 2) %>%
    pull("score") %>%
    min()
  
  threshold3 <- post.df3 %>% 
    dplyr::filter(Component == 5) %>%
    pull("score") %>%
    min()
  
  threshold4 <- post.df3 %>% 
    dplyr::filter(Component == 5) %>%
    pull("score") %>%
    max()
  
  return(c(threshold1, threshold2, threshold3, threshold4, peaks))
}


nuc_dist <- function(df) {
  df %>%
    arrange(seqnames, start) %>%
    group_by(transcriptId) %>%
    mutate(
      difference = lead(start) - start
    ) %>%
    dplyr::filter(!is.na(difference))
}

edb <- EnsDb.Hsapiens.v86
seqlevelsStyle(edb) <- "UCSC"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

dir <- "./ATAC_Seq/nucleoatac/"

nucs_files <- list(WT = paste0(dir, "WT_nucleoatac.nucpos.bed"), 
                   C5 = paste0(dir, "C5_nucleoatac.nucpos.bed"), 
                   C9 = paste0(dir, "C9_nucleoatac.nucpos.bed"))




annotated_nucs_wt <- annotatePeak(nucs_files[["WT"]], tssRegion = c(-500, 500),
                                  TxDb = edb, annoDb = "org.Hs.eg.db", overlap = "all")
annotated_nucs_wt_df <- as.data.frame(annotated_nucs_wt)
annotated_nucs_wt_df <- annotated_nucs_wt_df %>%
  dplyr::filter(grepl("Promoter", annotation)) %>%
  dplyr::filter(SYMBOL != "NA")



annotated_nucs_c5 <- annotatePeak(nucs_files[["C5"]], tssRegion = c(-500, 500),
                                  TxDb = edb, annoDb = "org.Hs.eg.db", overlap = "all")
annotated_nucs_c5_df <- as.data.frame(annotated_nucs_c5)
annotated_nucs_c5_df <- annotated_nucs_c5_df %>%
  dplyr::filter(grepl("Promoter", annotation)) %>%
  dplyr::filter(SYMBOL != "NA")


annotated_nucs_c9 <- annotatePeak(nucs_files[["C9"]], tssRegion = c(-500, 500),
                                  TxDb = edb, annoDb = "org.Hs.eg.db", overlap = "all")
annotated_nucs_c9_df <- as.data.frame(annotated_nucs_c9)
annotated_nucs_c9_df <- annotated_nucs_c9_df %>%
  dplyr::filter(grepl("Promoter", annotation)) %>%
  dplyr::filter(SYMBOL != "NA")



# Print the resulting dataframe
WT_nuc_dist <- annotated_nucs_wt_df %>%
  nuc_dist()
C5_nuc_dist <- annotated_nucs_c5_df %>%
  nuc_dist()
C9_nuc_dist <- annotated_nucs_c9_df %>%
  nuc_dist()

median(WT_nuc_dist$difference)
median(C5_nuc_dist$difference)
median(C9_nuc_dist$difference)

median(WT_nuc_dist$V13)
median(C5_nuc_dist$V13)
median(C9_nuc_dist$V13)

WT_nuc_dist$sample <- "WT" 
C5_nuc_dist$sample <- "C5"
C9_nuc_dist$sample <- "C9"


nuc_dist_merged <- rbind(WT_nuc_dist, C5_nuc_dist) %>%
  rbind(C9_nuc_dist)

nuc_dist_merged$sample <- factor(nuc_dist_merged$sample,
                                 levels = c("WT", "C5", "C9"), ordered = TRUE)


WT_thresholds_mixed_model <- nuc_dist_merged %>%
  mixed_model("WT")
C5_thresholds_mixed_model <- nuc_dist_merged %>%
  mixed_model("C5")
C9_thresholds_mixed_model <- nuc_dist_merged %>%
  mixed_model("C9")


nuc_dist_merged <- nuc_dist_merged %>%
  dplyr::mutate(nuc_type = case_when(sample == "WT" & 
                                       distanceToTSS < WT_thresholds_mixed_model[3] &
                                       distanceToTSS >= 0 ~ "+1 nuc",
                                     sample == "WT" & 
                                       distanceToTSS >= WT_thresholds_mixed_model[3] ~ "+2 nuc",
                                     sample == "WT" & 
                                       distanceToTSS >= WT_thresholds_mixed_model[1] &
                                       distanceToTSS < 0 ~ "-1 nuc",
                                     sample == "WT" & 
                                       distanceToTSS < WT_thresholds_mixed_model[1] ~ "-2 nuc",
                                     sample == "C5" & 
                                       distanceToTSS < C5_thresholds_mixed_model[3] &
                                       distanceToTSS >= 0 ~ "+1 nuc",
                                     sample == "C5" & 
                                       distanceToTSS >= C5_thresholds_mixed_model[3] ~ "+2 nuc",
                                     sample == "C5" & 
                                       distanceToTSS >= C5_thresholds_mixed_model[1] &
                                       distanceToTSS < 0 ~ "-1 nuc",
                                     sample == "C5" & 
                                       distanceToTSS < C5_thresholds_mixed_model[1] ~ "-2 nuc",
                                     sample == "C9" & 
                                       distanceToTSS < C9_thresholds_mixed_model[3] &
                                       distanceToTSS >= 0 ~ "+1 nuc",
                                     sample == "C9" & 
                                       distanceToTSS >= C9_thresholds_mixed_model[3] ~ "+2 nuc",
                                     sample == "C9" & 
                                       distanceToTSS >= C9_thresholds_mixed_model[1] &
                                       distanceToTSS < 0 ~ "-1 nuc",
                                     sample == "C9" & 
                                       distanceToTSS < C9_thresholds_mixed_model[1] ~ "-2 nuc"
                                    ))



inter_dyad_dist_plot <- nuc_dist_merged %>%
  ggplot(aes(x = sample, y = difference, fill = sample)) +
  geom_boxplot(linewidth = 1, alpha = 0.8) +
  scale_fill_manual(values = c("grey40", "#E69F00", "#56B4E9")) +
  #coord_cartesian(y = c(100, 500)) +
  theme_classic(base_size = 26) +
  theme(legend.position = "none") +
  xlab("Sample") +
  ylab("Inter-dyad distance (bp)") + 
  stat_compare_means(method= "wilcox.test", comparisons = nucs_comparisons, label.y = c(950, 1000), size = 6,
                     bracket.size = 1, alpha = 0.2, tip.length = 0.01) + 
  annotate("text", label = median(WT_nuc_dist$difference), x = 1, y = 300, size = 5.5) +
  annotate("text", label = median(C5_nuc_dist$difference), x = 2, y = 290, size = 5.5) + 
  annotate("text", label = median(C9_nuc_dist$difference), x = 3, y = 290, size = 5.5) +
  scale_y_continuous(breaks = c(200, 400, 600, 800, 1000)) 


ggsave(filename = "./ATAC_Seq/plots/inter_dyad_dist.pdf", plot = inter_dyad_dist_plot,
       width = 4.5, height = 6.5, dpi = 1000, device = cairo_pdf())

median(WT_nuc_dist$difference)
median(C5_nuc_dist$difference)
median(C9_nuc_dist$difference)

nuc_dist_merged %>%
  ggplot(aes(x = expression_c5, y = V13, colour = sample)) +
  geom_boxplot(outlier.colour = NA, linewidth = 1.2) +
  scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9")) +
  #coord_cartesian(y = c(0, 500)) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none") +
  xlab("Sample") +
  ylab("Nucleosome occupancy score")


nuc_dist_merged %>%
  ggplot(aes(x = distanceToTSS, colour = sample, y = V5)) +
  #geom_line() +
  #geom_line(stat = "summary_bin", binwidth = 200) +
  geom_smooth(se = FALSE, method = "gam") +
  scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9")) +
  theme_classic() +
  xlim(c(-700, 700)) + 
  geom_vline(xintercept = WT_thresholds_mixed_model[2], 
             colour = "#000000", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = C5_thresholds_mixed_model[2], 
             colour = "#E69F00", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = C9_thresholds_mixed_model[2], 
             colour = "#56B4E9", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = WT_thresholds_mixed_model[3], 
             colour = "#000000", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = C5_thresholds_mixed_model[3], 
             colour = "#E69F00", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = C9_thresholds_mixed_model[3], 
             colour = "#56B4E9", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = WT_thresholds_mixed_model[3], 
             colour = "#000000", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = C5_thresholds_mixed_model[3], 
             colour = "#E69F00", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = C9_thresholds_mixed_model[3], 
             colour = "#56B4E9", linetype = "dashed", linewidth = 1)


nucs_density <- nuc_dist_merged %>%
  ggplot(aes(x = distanceToTSS, fill = sample, colour = sample)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9")) +
  scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9")) +
  theme_bw(base_size = 20) +
  theme(panel.grid = element_blank(), legend.position = "none", plot.margin = margin(r = 0, l = 0)) +
  labs(x = "Distance to TSS (bp)", y = "Nucleosome density") +
  facet_wrap(~sample, nrow = 3) +
  ylim(c(0, 0.0023)) +
  xlim(c(-500, 500)) +
  geom_vline(data=dplyr::filter(nuc_dist_merged, sample=="WT"), 
             aes(xintercept = WT_thresholds_mixed_model[5]), 
             colour = "#000000", linetype = "dashed", linewidth = 1) +
  geom_vline(data=dplyr::filter(nuc_dist_merged, sample=="C5"), 
             aes(xintercept = C5_thresholds_mixed_model[5]), 
             colour = "#E69F00", linetype = "dashed", linewidth = 1) +
  geom_vline(data=dplyr::filter(nuc_dist_merged, sample=="C9"), 
             aes(xintercept = C9_thresholds_mixed_model[5]), 
             colour = "#56B4E9", linetype = "dashed", linewidth = 1) +
  geom_vline(data=dplyr::filter(nuc_dist_merged, sample=="WT"), 
             aes(xintercept = WT_thresholds_mixed_model[6]), 
             colour = "#000000", linetype = "dashed", linewidth = 1) +
  geom_vline(data=dplyr::filter(nuc_dist_merged, sample=="C5"), 
             aes(xintercept = C5_thresholds_mixed_model[6]), 
             colour = "#E69F00", linetype = "dashed", linewidth = 1) +
  geom_vline(data=dplyr::filter(nuc_dist_merged, sample=="C9"), 
             aes(xintercept = C9_thresholds_mixed_model[6]), 
             colour = "#56B4E9", linetype = "dashed", linewidth = 1) +
  geom_vline(data=dplyr::filter(nuc_dist_merged, sample=="WT"), 
             aes(xintercept = WT_thresholds_mixed_model[7]), 
             colour = "#000000", linetype = "dashed", linewidth = 1) +
  geom_vline(data=dplyr::filter(nuc_dist_merged, sample=="C5"), 
             aes(xintercept = C5_thresholds_mixed_model[7]), 
             colour = "#E69F00", linetype = "dashed", linewidth = 1) +
  geom_vline(data=dplyr::filter(nuc_dist_merged, sample=="C9"), 
             aes(xintercept = C9_thresholds_mixed_model[7]), 
             colour = "#56B4E9", linetype = "dashed", linewidth = 1) +
  geom_vline(data=dplyr::filter(nuc_dist_merged, sample=="WT"), 
             aes(xintercept = WT_thresholds_mixed_model[9]), 
             colour = "#000000", linetype = "dashed", linewidth = 1) +
  geom_vline(data=dplyr::filter(nuc_dist_merged, sample=="C5"), 
             aes(xintercept = C5_thresholds_mixed_model[9]), 
             colour = "#E69F00", linetype = "dashed", linewidth = 1) +
  geom_vline(data=dplyr::filter(nuc_dist_merged, sample=="C9"), 
             aes(xintercept = C9_thresholds_mixed_model[9]), 
             colour = "#56B4E9", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = 0), 
             colour = "#BF4342", linetype = "dashed", linewidth = 1) +
  geom_text(data=dplyr::filter(nuc_dist_merged, sample=="WT"), #+1 nuc
            label = WT_thresholds_mixed_model[7], 
            colour = "#000000", x = mean(c(WT_thresholds_mixed_model[7], 0)), y = 0.0021, size = 5.5) +
  geom_text(data=dplyr::filter(nuc_dist_merged, sample=="C5"), 
            label = C5_thresholds_mixed_model[7], 
            colour = "#E69F00", x = mean(c(C5_thresholds_mixed_model[7], 0)), y = 0.0021, size = 5.5) +
  geom_text(data=dplyr::filter(nuc_dist_merged, sample=="C9"), 
            label = C9_thresholds_mixed_model[7], 
            colour = "#56B4E9", x = mean(c(C9_thresholds_mixed_model[7], 0)), y = 0.0021, size = 5.5) +
  geom_text(data=dplyr::filter(nuc_dist_merged, sample=="WT"), #+2 nuc
            label = WT_thresholds_mixed_model[9] - WT_thresholds_mixed_model[7], 
            colour = "#000000", x = mean(c(WT_thresholds_mixed_model[9], WT_thresholds_mixed_model[7])), y = 0.0021, size = 5.5) +
  geom_text(data=dplyr::filter(nuc_dist_merged, sample=="C5"), 
            label = C5_thresholds_mixed_model[9] - C5_thresholds_mixed_model[7], 
            colour = "#E69F00", x = mean(c(C5_thresholds_mixed_model[9], C5_thresholds_mixed_model[7])), y = 0.0021, size = 5.5) +
  geom_text(data=dplyr::filter(nuc_dist_merged, sample=="C9"), 
            label = C9_thresholds_mixed_model[9] - C9_thresholds_mixed_model[7], 
            colour = "#56B4E9", x = mean(c(C9_thresholds_mixed_model[9], C9_thresholds_mixed_model[7])), y = 0.0021, size = 5.5) +
  geom_text(data=dplyr::filter(nuc_dist_merged, sample=="WT"), #-2 nuc
            label = WT_thresholds_mixed_model[5] - WT_thresholds_mixed_model[6], 
            colour = "#000000", x = mean(c(WT_thresholds_mixed_model[5], WT_thresholds_mixed_model[6])), y = 0.0021, size = 5.5) +
  geom_text(data=dplyr::filter(nuc_dist_merged, sample=="C5"), 
            label = C5_thresholds_mixed_model[5] - C5_thresholds_mixed_model[6], 
            colour = "#E69F00", x = mean(c(mean(c(C5_thresholds_mixed_model[5], C5_thresholds_mixed_model[6])))), y = 0.0021, size = 5.5) +
  geom_text(data=dplyr::filter(nuc_dist_merged, sample=="C9"), 
            label = C9_thresholds_mixed_model[5] - C9_thresholds_mixed_model[6], 
            colour = "#56B4E9", x = mean(c(C9_thresholds_mixed_model[5], C9_thresholds_mixed_model[6])), y = 0.0021, size = 5.5) +
  geom_text(data=dplyr::filter(nuc_dist_merged, sample=="WT"), #-1 nuc
            label = WT_thresholds_mixed_model[5], 
            colour = "#000000", x = mean(c(0, WT_thresholds_mixed_model[5])), y = 0.0021, size = 5.5) +
  geom_text(data=dplyr::filter(nuc_dist_merged, sample=="C5"), 
            label = C5_thresholds_mixed_model[5], 
            colour = "#E69F00", x = mean(c(mean(c(0, C5_thresholds_mixed_model[5])))), y = 0.0021, size = 5.5) +
  geom_text(data=dplyr::filter(nuc_dist_merged, sample=="C9"), 
            label = C9_thresholds_mixed_model[5], 
            colour = "#56B4E9", x = mean(c(0, C9_thresholds_mixed_model[5])), y = 0.0021, size = 5.5, size = 5.5)


ggsave(filename = "./ATAC_Seq/plots/nucs_density.png", plot = nucs_density,
       width = 8.7, height = 6.5, dpi = 1000)

nucs_comparisons <- list(c("WT", "C5"), c("WT", "C9"))
nucleosome_fuziness_plot <- nuc_dist_merged %>%
  ggplot(aes(x = sample, y = V13, fill = sample)) +
  coord_cartesian(y = c(15, 60)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("grey40", "#E69F00", "#56B4E9")) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none") +
  xlab("Sample") +
  ylab("Nucleosome fuzziness") +
  facet_wrap(~nuc_type, nrow = 1, labeller = as_labeller(c("-1 nuc" = "-1 nucleosome", 
                                                 "-2 nuc" = "-2 nucleosome",
                                                 "+1 nuc" = "+1 nucleosome", 
                                                 "+2 nuc" = "+2 nucleosome"))) + 
  stat_compare_means(method= "wilcox.test", comparisons = nucs_comparisons) +
  theme(panel.grid = element_blank())

ggsave(filename = "./ATAC_Seq/plots/nucleosome_fuzzines_per_nuc_1row.pdf", plot = nucleosome_fuziness_plot,
       width = 10, height = 3.5, dpi = 1000, device = cairo_pdf())

nus_occupancy_score_plot <- nuc_dist_merged %>%
  ggplot(aes(x = sample, y = V5, fill = sample)) +
  geom_boxplot(outlier.colour = NA) +
  scale_fill_manual(values = c("grey40", "#E69F00", "#56B4E9")) +
  theme_bw(base_size = 20) +
  coord_cartesian(y = c(0.2, 1.25)) +
  xlab("Sample") +
  ylab("Nucleosome occupancy score") +
  facet_wrap(~nuc_type, labeller = as_labeller(c("-1 nuc" = "-1 nucleosome", 
                                                 "-2 nuc" = "-2 nucleosome",
                                                 "+1 nuc" = "+1 nucleosome", 
                                                 "+2 nuc" = "+2 nucleosome"))) +
  stat_compare_means(method= "wilcox.test", comparisons = nucs_comparisons) +
  theme(panel.grid = element_blank(), legend.position = "none") + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1.25))

ggsave(filename = "./ATAC_Seq/plots/nucleosome_ocuupancy_score_per_nuc.pdf", plot = nus_occupancy_score_plot,
        width = 4.8, height = 6, dpi = 1000, device = cairo_pdf())