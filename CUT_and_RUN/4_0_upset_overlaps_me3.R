library(hicVennDiagram)
library(GenomicInteractions)
library(ComplexUpset)
library(tidyverse)
library(ChIPseeker)

fs <- dir("./CUT_and_RUN/1_2_Snakemake_run_merged_replicates/seacr/filtered", pattern = "_me3.filtered_stranded.bed", full.names = TRUE)


venn <- vennCount(fs, maxgap = 1000, FUN = min)


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
colnames(plotdata)
names(plotdata) <- c("OCI AML2 C5 H3K27me3", "OCI AML2 C9 H3K27me3", "OCI AML2 WT H3K27me3")

plotdata <- plotdata %>%
  dplyr::select("OCI AML2 C5 H3K27me3", "OCI AML2 C9 H3K27me3", "OCI AML2 WT H3K27me3")



venn_cnr <- ComplexUpset::upset(
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
    upset_query(set="OCI AML2 WT H3K27me3",fill="black"),
    upset_query(set="OCI AML2 C5 H3K27me3",fill="#E69F00"),
    upset_query(set = "OCI AML2 C9 H3K27me3", fill = "#56B4E9")), 
  height_ratio=0.7,
  width_ratio=0.2,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE))) + 
  labs(x = "Peaks")


wt_me3_unique_peaks <- venn$overlapList[["001"]][["WT_me3_filtered_stranded.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

c9_me3_unique_peaks <- venn$overlapList[["010"]][["C9_me3_filtered_stranded.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

c5_me3_unique_peaks <- venn$overlapList[["100"]][["C5_me3_filtered_stranded.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

wt_c5_c9_me3_overlap_peaks <- venn$overlapList[["111"]][["WT_me3_filtered_stranded.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

c5_c9_me3_overlap_peaks <- venn$overlapList[["110"]][["C5_me3_filtered_stranded.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

create.dir("./CUT_and_RUN/1_2_Snakemake_run_merged_replicates/contrasts")
write.table(wt_me3_unique_peaks, file = "./CUT_and_RUN/1_2_Snakemake_run_merged_replicates/seacr/contrasts/WT_me3_unique_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(c5_me3_unique_peaks, file = "./CUT_and_RUN/1_2_Snakemake_run_merged_replicates/seacr/contrasts/C5_me3_unique_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(c9_me3_unique_peaks, file = "./CUT_and_RUN/1_2_Snakemake_run_merged_replicates/seacr/contrasts/C9_me3_unique_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(wt_c5_c9_me3_overlap_peaks, file = "./CUT_and_RUN/1_2_Snakemake_run_merged_replicates/seacr/contrasts/wt_c5_c9_me3_overlap_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(c5_c9_me3_overlap_peaks, file = "./CUT_and_RUN/1_2_Snakemake_run_merged_replicates/seacr/contrasts/c5_c9_me3_overlap_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

rbind(c5_me3_unique_peaks, c9_me3_unique_peaks) %>%
  rbind(c5_c9_me3_overlap_peaks) %>%
  write.table(file = "./CUT_and_RUN/1_2_Snakemake_run_merged_replicates/seacr/contrasts/c5_c9_me3_union_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


rbind(wt_me3_unique_peaks, c5_me3_unique_peaks) %>%
  rbind(c9_me3_unique_peaks) %>%
  rbind(c5_c9_me3_overlap_peaks) %>%
  rbind(wt_c5_c9_me3_overlap_peaks) %>%
  write.table(file = "./CUT_and_RUN/1_2_Snakemake_run_merged_replicates/seacr/contrasts/all_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


#rbind(wt_me3_unique_peaks, c9_me3_unique_peaks) %>%
 # rbind(wt_c9_me3_overlap_peaks) %>%
  #write.table(file = "./CUT_and_RUN/1_Snakemake_run/seacr/contrasts/me3_all.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)