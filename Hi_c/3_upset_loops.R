library(tidyverse)
library(hicVennDiagram)
library(GenomicInteractions)
library(ComplexUpset)
library(stringr)
library(tidyr)
library(ChIPseeker)
library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#transform loop calls into bedpe
dir.create("./Hi_C/1_mustache_loop_calling/bed/")
read.table("./Hi_C/1_mustache_loop_calling/OCI_AML2_WT_15000_1.6_mustache_loops.tsv", sep = "\t", 
                              header = TRUE)[1:6] %>%
  write.table("./Hi_C/1_mustache_loop_calling/bed/OCI_AML2_WT_15000_1.6_mustache_loops.bedpe", 
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

read.table("./Hi_C/1_mustache_loop_calling/OCI_AML3_WT_15000_1.6_mustache_loops.tsv", sep = "\t", 
           header = TRUE)[1:6] %>%
  write.table("./Hi_C/1_mustache_loop_calling/bed/OCI_AML3_WT_15000_1.6_mustache_loops.bedpe", 
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


read.table("./Hi_C/1_mustache_loop_calling/OCI_AML2_C9_15000_1.6_mustache_loops.tsv", sep = "\t", 
           header = TRUE)[1:6] %>% 
  write.table("./Hi_C/1_mustache_loop_calling/bed/OCI_AML2_C9_15000_1.6_mustache_loops.bedpe", 
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

all_files <- Sys.glob("./Hi_C/1_mustache_loop_calling/bed/*_15000_1.6_mustache_loops.bedpe")
interactions <- lapply(all_files, makeGenomicInteractionsFromFile, type="bedpe") %>% 
  magrittr::set_names(c("OCI_AML2_C9", "OCI_AML2_WT", "OCI_AML3_WT"))

fs <- dir("./Hi_C/1_mustache_loop_calling/bed/", pattern = "_1.6_mustache_loops.bedpe", full.names = TRUE)

venn <- vennCount(fs, maxgap = 20000, FUN = max)

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

names(plotdata) <- c("OCI AML2 C9", "OCI AML2 WT", "OCI AML3 WT")

plotdata <- plotdata %>%
  dplyr::select("OCI AML2 C9","OCI AML2 WT", "OCI AML3 WT")

upset_hic <- ComplexUpset::upset(
  data=as.data.frame(plotdata),
  intersect=colnames(plotdata),
  themes = upset_themes_fix,
  min_size = 50,
  sort_sets = FALSE,
  queries = list(
    upset_query(set="OCI AML2 WT",fill="black"),
    upset_query(set = "OCI AML3 WT", fill = "grey60"),
    upset_query(set = "OCI AML2 C9", fill = "#56B4E9")
  ),
  height_ratio=1,
  width_ratio=0.2,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=FALSE)))

ggsave(filename = "./Hi_c/upset_loops.pdf", plot = upset_hic, 
       width = 16, height = 10, dpi = 800, units = "cm", device = cairo_pdf)

OCI_AML2_WT_loops <- venn$overlapList[["010"]][["OCI_AML2_WT_15000_1.6_mustache_loops.bedpe"]]
OCI_AML3_WT_loops <- venn$overlapList[["001"]][["OCI_AML3_WT_15000_1.6_mustache_loops.bedpe"]]
OCI_AML2_WT_AML3_WT_loops <- venn$overlapList[["011"]][["OCI_AML2_WT_15000_1.6_mustache_loops.bedpe"]]
OCI_AML2_C9_loops <- venn$overlapList[["100"]][["OCI_AML2_C9_15000_1.6_mustache_loops.bedpe"]]
all_samples_loops <- venn$overlapList[["111"]][["OCI_AML2_WT_15000_1.6_mustache_loops.bedpe"]]

export.bedpe(OCI_AML2_WT_loops, fn = "./Hi_C/results_files/OCI_AML2_WT_loops.bedpe")
export.bedpe(OCI_AML3_WT_loops, fn = "./Hi_C/results_files/OCI_AML3_WT_loops.bedpe")
export.bedpe(OCI_AML2_WT_AML3_WT_loops, fn = "./Hi_C/results_files/OCI_AML2_WT_AML3_WT_loops.bedpe")
export.bedpe(OCI_AML2_C9_loops, fn = "./Hi_C/results_files/OCI_AML2_C9_loops.bedpe")
export.bedpe(all_samples_loops, fn = "./Hi_C/results_files/all_samples_loops.bedpe")


all_called_loops <- as.data.frame(venn$overlapList[["001"]][[1]]) %>%
  rbind(as.data.frame(venn$overlapList[["010"]][[1]])) %>%
  rbind(as.data.frame(venn$overlapList[["100"]][[1]])) %>%
  rbind(as.data.frame(venn$overlapList[["011"]][[1]])) %>%
  rbind(as.data.frame(venn$overlapList[["110"]][[1]])) %>%
  rbind(as.data.frame(venn$overlapList[["101"]][[1]])) %>%
  rbind(as.data.frame(venn$overlapList[["111"]][[1]])) %>%
  select(c("seqnames1", "start1", "end1", "seqnames2", "start2", "end2" ))

write.table(all_called_loops, file = "./Hi_C/results_files/all_called_loops.bedpe", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")