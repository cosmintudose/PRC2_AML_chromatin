library(hicVennDiagram)
library(GenomicInteractions)
library(ComplexUpset)
library(tidyverse)
library(ChIPseeker)

fs <- dir("./CUT_and_RUN/1_Snakemake_run/seacr", 
          pattern = "_stranded.bed", full.names = TRUE)

#Overlap peak files
venn <- vennCount(fs, maxgap = 1000, FUN = min) #considering peaks as overlapping if they're <1kb apart

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
names(plotdata) <- c("OCI AML2 C9 H2AK119ub", "OCI AML2 C9 H3K27me3", "OCI AML2 WT H2AK119ub", "OCI AML2 WT H3K27me3")

plotdata <- plotdata %>%
  dplyr::select("OCI AML2 C9 H2AK119ub", "OCI AML2 WT H2AK119ub", "OCI AML2 C9 H3K27me3", "OCI AML2 WT H3K27me3")


#Make upset plot
venn_cnr <- ComplexUpset::upset(
  data=as.data.frame(plotdata),
  intersect=colnames(plotdata),
  themes = upset_modify_themes(list('overall_sizes'=theme(axis.text.x=element_text(angle=90), 
                                                          panel.grid = element_blank()), 
                                    'Intersection size'=theme(panel.grid = element_blank(),
                                                              text=element_text(size=14)),
                                    'intersections_matrix'=theme(panel.grid = element_blank(),
                                                                 text=element_text(size=14)))),
  min_size = 900,
  sort_sets = FALSE,
  queries = list(
    upset_query(set="OCI AML2 WT H2AK119ub",fill="black"),
    upset_query(set = "OCI AML2 C9 H2AK119ub", fill = "#56B4E9"),
    upset_query(set="OCI AML2 WT H3K27me3",fill="black"),
    upset_query(set = "OCI AML2 C9 H3K27me3", fill = "#56B4E9")), 
  height_ratio=0.5,
  width_ratio=0.2,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE))) + 
  labs(x = "Peaks")

ggsave(venn_cnr, file = "./CUT_and_RUN/plots/upset_overlaps.pdf",
       width = 26, height = 9, dpi = 800, units = "cm", device = cairo_pdf)


#extracting and saving peaks overlapping and condition-specific
wt_me3_ub_unique_peaks <- venn$overlapList[["0011"]][["WT_h2ak119ub_stranded.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

c9_me3_ub_unique_peaks <- venn$overlapList[["1100"]][["C9_h2ak119ub_stranded.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))
  
wt_c9_me3_ub_overlap_peaks <- venn$overlapList[["1111"]][["WT_h2ak119ub_stranded.bed"]] %>%
  as.data.frame() %>%
  dplyr::select(c("seqnames", "start", "end"))

dir.create("./CUT_and_RUN/1_Snakemake_run/seacr/contrasts")

write.table(wt_me3_ub_unique_peaks, file = "./CUT_and_RUN/1_Snakemake_run/seacr/contrasts/WT_me3_ub_unique_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(c9_me3_ub_unique_peaks, file = "./CUT_and_RUN/1_Snakemake_run/seacr/contrasts/C9_me3_ub_unique_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(wt_c9_me3_ub_overlap_peaks, file = "./CUT_and_RUN/1_Snakemake_run/seacr/contrasts/WT_C9_me3_ub_overlap_peaks.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

rbind(wt_me3_ub_unique_peaks, c9_me3_ub_unique_peaks) %>%
  rbind(wt_c9_me3_ub_overlap_peaks) %>%
  write.table(file = "./CUT_and_RUN/1_Snakemake_run/seacr/contrasts/me3_ub_all.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)





