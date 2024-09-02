###aml2 wt vs 
###load packages ----
library(tidyverse)
library(tximport)
library(ensembldb) 
library(EnsDb.Hsapiens.v86)
library(edgeR)
library(matrixStats)
library(cowplot)
library(plotly)
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

### read data----
targets.aml2 <- read_tsv("./RNA_Seq/study_design.txt") %>%
  dplyr::select(-phenotype_sep)
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
                         aml2.clones.AVG = (C5_1 + C5_2 + C5_3 + C9_1 + C9_2 + C9_3)/6,
                         #now make columns comparing each of the averages above 
                         LogFC = (aml2.clones.AVG - aml2.wt.AVG))%>% 
  mutate_if(is.numeric, round, 2)


write.csv(mydata.aml2.df, file = "./RNA_Seq/results_files/normalised_expression_aml2_wt_vs_clones.csv", row.names = FALSE)


groupaml2 <- targets.aml2$phenotype
groupaml2 <- factor(groupaml2)


###diff genes - volcano plot ----

designaml2 <- model.matrix(~0 + groupaml2)
colnames(designaml2) <- levels(groupaml2)

v.DEGList.aml2.filtered.norm <- voom(DGEList.filtered.norm, designaml2, plot = TRUE) ###models mean-variance relationship
fit.aml2 <- lmFit(v.DEGList.aml2.filtered.norm, designaml2)
contrast.matrix.aml2 <- makeContrasts(AML2_c - AML2_wt,
                                      levels = designaml2)

fits.aml2 <- contrasts.fit(fit.aml2, contrast.matrix.aml2)
ebFit.aml2 <- eBayes(fits.aml2)
myTopHits.aml2 <- topTable(ebFit.aml2, adjust ="BH", coef=1, number=20000, sort.by="logFC")

myTopHits.aml2.df <- myTopHits.aml2 %>%
  as_tibble(rownames = "geneID")

write.csv(myTopHits.aml2.df, file = "./RNA_Seq/results_files/aml2_wt_v_clones_12k_genes.csv", quote = FALSE, row.names = FALSE)


signif.genes <- myTopHits.aml2.df %>%
  filter(adj.P.Val < 0.1)# & (logFC <= -0.5 | logFC >= 0.5))

results.aml2 <- decideTests(ebFit.aml2, method="global", adjust.method="BH", p.value=0.1, lfc=0.5)

summary(results.aml2)


vplot <- ggplot(myTopHits.aml2.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  annotate("rect", xmin = 1, xmax = 10, ymin = -log10(0.1), ymax = 6.35, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -1, xmax = -10, ymin = -log10(0.1), ymax = 6.35, alpha=.2, fill="#2C467A") +
  annotate(geom = "text", x = -8, y = 1.5, label = "175", size = 9, colour = "#2C467A") +
  annotate(geom = "text", x = 8, y = 1.5, label = "74", size = 9, colour = "#BE684D") +
  labs(title="AML2 WT & clones",
       subtitle = "Volcano plot") +
  theme_bw(base_size = 24) +
  theme(axis.text = element_text(size = 24)) + 
  labs(y = expression(-log[10]("FDR"))) + 
  theme(panel.grid.minor = element_blank())

vplot

ggsave(vplot, 
       filename = "./RNA_Seq/plots/vplot_aml2_wt_vs_clones.pdf",
       width = 18, height = 13, dpi = 1000, units = "cm", device = cairo_pdf())



# GSEA ----
hs_gsea_h <- msigdbr(species = "Homo sapiens",
                      category = "H") %>% #msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) #just columns corresponding to signature name and gene symbols in each signature 

# columns corresponding to gene symbols and LogFC for  one pairwise comparison for the enrichment analysis
mydata.aml2.df.sub <- dplyr::select(mydata.aml2.df, geneID, LogFC)
mydata.aml2.gsea <- mydata.aml2.df.sub$LogFC
names(mydata.aml2.gsea) <- as.character(mydata.aml2.df.sub$geneID)
mydata.aml2.gsea <- sort(mydata.aml2.gsea, decreasing = TRUE)

# GSEA function from clusterProfiler
set.seed(1234)
myGSEA.aml2.res <- GSEA(mydata.aml2.gsea, TERM2GENE=hs_gsea_h, 
                        verbose=FALSE, seed = TRUE, nPermSimple = 10000, eps = 0)
myGSEA.aml2.df <- as_tibble(myGSEA.aml2.res@result)


# NES graph for any each signature
gseaplot2(myGSEA.aml2.res, 
          geneSetID = 1, #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, 
          title = myGSEA.aml2.res$Description[1]) #can turn off


myGSEA.aml2.df <- myGSEA.aml2.df %>%
  mutate(phenotype = dplyr::case_when(
    NES < 0 ~ "EZH2+/+\n(WT)",
    NES > 0 ~ "EZH2+/-\n(C5&C9)")) %>%
  rename("FDR" = "p.adjust")

myGSEA.aml2.df$phenotype <- factor(myGSEA.aml2.df$phenotype, 
                                   levels = c("EZH2+/+\n(WT)", "EZH2+/-\n(C5&C9)"))

gsea_barplot <- myGSEA.aml2.df %>%
  ggplot(aes(x=NES, y= fct_reorder(ID, NES), fill = NES)) +
  geom_bar(stat = "identity") +
  scale_fill_distiller(palette = "RdBu") +
  labs(y = "ID") +
  theme_minimal(base_size = 16) + 
  scale_alpha_continuous(limits = c(0, 3)) +
  scale_x_continuous(limits = c(-2, 2)) +
  theme(axis.title = element_text(size = 22), 
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 18), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank(), legend.box = "horizontal")

ggsave(gsea_barplot, 
       filename = "./RNA_Seq/plots/gsea_barplot_clones_hallmarks.pdf",
       width = 27, height = 7, dpi = 500, units = "cm", device = cairo_pdf)