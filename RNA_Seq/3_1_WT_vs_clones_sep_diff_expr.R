###aml2 wt vs 
###load packages ----
library(tidyverse)
library(tximport)
library(ensembldb) 
library(EnsDb.Hsapiens.v86)
library(edgeR)
library(matrixStats)
library(limma)
library(gplots)
library(RColorBrewer)
library(GSEABase) 
library(Biobase) 
library(gprofiler2) 
library(clusterProfiler) 
library(msigdbr) 
library(enrichplot) 
library(colorspace)
library(Cairo)
library(VennDiagram)


### read data----
targets.aml2 <- read_tsv("./RNA_Seq/study_design.txt") %>%
  dplyr::select(-phenotype)
path <- file.path("./RNA_Seq/", targets.aml2$code, "abundance.tsv") # set file paths to mapped data

Tx.aml2 <- transcripts(EnsDb.Hsapiens.v86, columns = c("tx_id", "gene_name")) # annotations and gene symbols 
Tx.aml2 <- as_tibble(Tx.aml2)
Tx.aml2 <- dplyr::rename(Tx.aml2, target_id = tx_id)
Tx.aml2 <- dplyr::select(Tx.aml2, "target_id", "gene_name")
Txi_gene.aml2 <- tximport(path, 
                          type = "kallisto", 
                          tx2gene = Tx.aml2, 
                          txOut = FALSE, # determines whether data represented at transcript or gene level (here read at gene level)
                          countsFromAbundance = "lengthScaledTPM",
                          ignoreTxVersion = TRUE)

### preprocessing----
sample_labels <- targets.aml2$sample
DGEList <- DGEList(Txi_gene.aml2$counts)
log2.cpm <- cpm(DGEList, log=TRUE)  ###log2 normalised counts per million

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID") ###convert as data frame
colnames(log2.cpm.df) <- c("geneID", sample_labels)

###tidy data
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df,
                                  cols = c(WT_1, WT_2, WT_3, C5_1, C5_2, C5_3, C9_1, C9_2, C9_3),
                                  names_to = "samples", # name of new column
                                  values_to = "expression") # name of new column storing all the data



###filter data
cpm <- cpm(DGEList)
keepers <- rowSums(cpm>1)>=3 #user defined - depends on studydesign, eliminates lowly expressed genes
DGEList.filtered <- DGEList[keepers,]

###normalize data
DGEList.filtered.norm <- calcNormFactors(DGEList.filtered, method = "TMM") #TMM normalization
log2.cpm.filtered.norm <- DGEList.filtered.norm %>% 
  cpm(log=TRUE) %>% 
  as_tibble(rownames = "geneID")

colnames(log2.cpm.filtered.norm) <- c("geneID", sample_labels)


### multivariate analysis ----

mydata.aml2.df <- mutate(log2.cpm.filtered.norm,
                         aml2.wt.AVG = (WT_1 + WT_2 + WT_3)/3, 
                         aml2.c5.AVG = (C5_1 + C5_2 + C5_3)/3,
                         aml2.c9.AVG = (C9_1 + C9_2 + C9_3)/3,
                         #now make columns comparing each of the averages above 
                         LogFC_c5_wt = (aml2.c5.AVG - aml2.wt.AVG), 
                         LogFC_c9_wt = (aml2.c9.AVG - aml2.wt.AVG))%>% 
  mutate_if(is.numeric, round, 2)

mydata.aml2.df %>%
  dplyr::select(-c("C9_1", "C9_2", "C9_3", "aml2.c9.AVG", "LogFC_c9_wt")) %>%
  write.csv(file = "./RNA_Seq/results_files/normalised_expression_aml2_wt_vs_c5.csv", row.names = FALSE)

mydata.aml2.df %>%
  dplyr::select(-c("C5_1", "C5_2", "C5_3", "aml2.c5.AVG", "LogFC_c5_wt")) %>%
  write.csv(file = "./RNA_Seq/results_files/normalised_expression_aml2_wt_vs_c9.csv", row.names = FALSE)

groupaml2 <- targets.aml2$phenotype_sep
groupaml2 <- factor(groupaml2)


###diff genes - volcano plot ----

designaml2 <- model.matrix(~0 + groupaml2)
colnames(designaml2) <- levels(groupaml2)

v.DEGList.aml2.filtered.norm <- voom(DGEList.filtered.norm, designaml2, plot = TRUE) ###models mean-variance relationship
fit.aml2 <- lmFit(v.DEGList.aml2.filtered.norm, designaml2)
contrast.matrix.aml2 <- makeContrasts(AML2_c5 - AML2_wt, 
                                      AML2_c9 - AML2_wt,
                                      levels = designaml2)

fits.aml2 <- contrasts.fit(fit.aml2, contrast.matrix.aml2)
ebFit.aml2 <- eBayes(fits.aml2)
TopHits.aml2_c5 <- topTable(ebFit.aml2, adjust ="BH", coef=1, number=20000, sort.by="logFC")
TopHits.aml2_c9 <- topTable(ebFit.aml2, adjust ="BH", coef=2, number=20000, sort.by="logFC")

TopHits.aml2_c5.df <- TopHits.aml2_c5 %>%
  as_tibble(rownames = "geneID")

TopHits.aml2_c9.df <- TopHits.aml2_c9 %>%
  as_tibble(rownames = "geneID")

write.csv(TopHits.aml2_c5.df, file = "./RNA_Seq/results_files/aml2_wt_v_c5_12k_genes.csv", quote = FALSE, row.names = FALSE)
write.csv(TopHits.aml2_c9.df, file = "./RNA_Seq/results_files/aml2_wt_v_c9_12k_genes.csv", quote = FALSE, row.names = FALSE)

#venn diagrams
c5_signif_genes <- TopHits.aml2_c5.df %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::filter(logFC <= -0.5 | logFC >= 0.5)

c9_signif_genes <- TopHits.aml2_c9.df %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::filter(logFC <= -0.5 | logFC >= 0.5)

c5_up <- c5_signif_genes %>%
  dplyr::filter(logFC >= 0.5) %>%
  dplyr::pull("geneID")
c9_up <- c9_signif_genes %>%
  dplyr::filter(logFC >= 0.5) %>%
  dplyr::pull("geneID")

c5_down <- c5_signif_genes %>%
  dplyr::filter(logFC <= -0.5) %>%
  dplyr::pull("geneID")
c9_down <- c9_signif_genes %>%
  dplyr::filter(logFC <= -0.5) %>%
  dplyr::pull("geneID")


genes.down <- list(C5 = c5_down, C9 = c9_down)
genes.up <- list(C5 = c5_up, C9 = c9_up)

venn.diagram(genes.down, lwd = 0, cex = 2, cat.cex = 3, print.mode = c("raw", "percent"),
             alpha = c(0.5, 0.5), fill = c("#a9d6e5", "#014f86"), 
             "genes.down.venn.tiff")
venn.diagram(genes.up, lwd = 0, cex = 2, cat.cex = 3, print.mode = c("raw", "percent"),
             alpha = c(0.5, 0.5), fill = c("#ffccd5", "#c9184a"), 
             "genes.up.venn.lfc1.tiff")
