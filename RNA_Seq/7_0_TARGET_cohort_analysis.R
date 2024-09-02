library(tidyverse)
library(readxl)
library(limma)
library(edgeR)
library(RColorBrewer)
library(GSEABase) 
library(Biobase) 
library(gprofiler2) 
library(clusterProfiler) 
library(msigdbr) 
library(gplots)
library(enrichplot)
library(ggpubr)

#TARGET RNA seq data
rna_seq <- read.table("./publicly_available_data/aml_target_2018_pub/data_mrna_seq_rpkm.txt", sep = "\t", header = TRUE, check.names = FALSE) %>%
  dplyr::select(-Entrez_Gene_Id)

#removing duplicate genes
rna_seq <- rna_seq[!duplicated(rna_seq$Hugo_Symbol), ]
rownames(rna_seq) <- rna_seq$Hugo_Symbol
rna_seq <- rna_seq %>%
  dplyr::select(-Hugo_Symbol) %>%
  na.omit()
rna_seq <- rna_seq[as.logical(rowSums(rna_seq != 0)), ]
log_rna_seq <- log2(rna_seq + 0.1) #log2(RPKM+1)


t_rna_seq <- log_rna_seq %>%
  t() %>%
  as.data.frame()

EZH2_target <- log_rna_seq["EZH2", ] %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "EZH2_mRNA")

EZH2_target %>%
  ggplot(aes(x = EZH2_mRNA)) +
  geom_histogram(aes(y=..density..), colour="black", fill="grey80")+
  geom_density(alpha=.2, fill="#FF6666") + 
  theme_minimal(base_size = 20)

EZH2_target %>% 
  summary()

#take median EZH2 expression value from summary
low_EZH2_samples <- EZH2_target %>%
  dplyr::filter(EZH2_mRNA <= 3.981) %>%
  pull("Sample")

high_EZH2_samples <- EZH2_target %>%
  dplyr::filter(EZH2_mRNA > 3.981) %>%
  pull("Sample")

#creating a design file to separate the sample in high and low EZH2 expressing
study_design_EZH2_expr <- data.frame(Sample = c(low_EZH2_samples, high_EZH2_samples)) %>%
  mutate(Condition_EZH2 = case_when(Sample %in% low_EZH2_samples ~ "EZH2_low",
                                    Sample %in% high_EZH2_samples ~ "EZH2_high"))

#only keep samples in study design
log_rna_seq <- log_rna_seq %>%
  dplyr::select(c(low_EZH2_samples, high_EZH2_samples))

group_EZH2 <- study_design_EZH2_expr$Condition_EZH2 %>%
  factor()

#differential expression analysis
design <- model.matrix(~0 + group_EZH2)
colnames(design) <- sub("group_EZH2", "", colnames(design))

#remove lowly expressed genes in at least 10% of samples
keepers <- rowSums(log_rna_seq>1)>= 4 
log_rna_seq <- log_rna_seq[keepers, ]


fit <- lmFit(log_rna_seq, design)
contrast.matrix <- makeContrasts(differences = EZH2_low - EZH2_high,
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)
TopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="p")

TopHits.df <- TopHits %>%
  as_tibble(rownames = "geneID")

write.csv(TopHits.df, file = "./RNA_Seq/results_files/differential_expression_EZH2_low_vs_high_TARGET.csv", row.names = FALSE, quote = FALSE)

# GSEA ----
# the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
#mydata.aml.df.sub <- dplyr::select(mydata.aml.df, geneID, LogFC)
mydata.gsea <- TopHits$logFC
names(mydata.gsea) <- as.character(rownames(TopHits))
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

#load atlas of blood cells signatures
atlas_blood_cells <- read.csv("./publicly_available_data/signatures_atlas_human_blood_cells.txt", sep = "\t") %>%
  dplyr::select(RNA_Cluster, Gene) %>%
  rename("gs_name" = "RNA_Cluster", "gene_symbol" = "Gene")

# GSEA function from clusterProfiler
set.seed(1234)
GSEA.res <- GSEA(mydata.gsea, TERM2GENE=atlas_blood_cells, verbose=FALSE, seed = TRUE, nPermSimple = 1000, pvalueCutoff = 0.1)
GSEA.df <- as_tibble(GSEA.res@result)

write.csv(GSEA.df, file = "./RNA_Seq/results_files/GSEA_altas_blood_cells_ezh2_low_vs_ezh2_high_TARGET.csv", row.names = FALSE, quote = FALSE)

gsea_plot <- GSEA.df %>%
  filter(p.adjust < 0.1) %>% #uncomment this if you want entire df
  ggplot(aes(x=NES, y= fct_reorder(ID, NES), fill = factor(sign(NES)))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#2C467A", "#BE684D")) +
  labs(y = "ID", title = "Gene signatures from cell lines\nEZH2 low vs EZH2 high (TARGET)", 
       caption = "*Median EZH2 expression used as stratification point") +
  theme_classic(base_size = 22) + 
  theme(legend.position = "none") + 
  annotate(geom = "text", x = -0.3, y = 20.5, size = 5, fontface = "bold",
           label = "Enriched in\nEZH2-low", colour = "#BE684D", angle = 90) +
  annotate(geom = "text", x = 0.3, y = 9, size = 5.5, fontface = "bold",
           label = "Enriched in\nEZH2-high", colour = "#2C467A", angle = 270)


ggsave(filename = "./RNA_Seq/plots/GSEA_TARGET_EZH2_low_vs_high.pdf", device = cairo_pdf, plot = gsea_plot, 
       width = 10, height = 8, dpi = 1000)