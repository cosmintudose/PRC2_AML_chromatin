###aml2 wt vs 
###load packages ----
library(tidyverse)
library(tximport)
library(ensembldb) 
library(EnsDb.Hsapiens.v86)
library(edgeR)
library(colorspace)
library(Cairo)


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
sampleLabels.aml2 <- targets.aml2$sample
myDGEList.aml2 <- DGEList(Txi_gene.aml2$counts)
log2.cpm.aml2 <- cpm(myDGEList.aml2, log=TRUE)  ###log2 normalised counts per million

log2.cpm.aml2.df <- as_tibble(log2.cpm.aml2, rownames = "geneID") ###convert as data frame
colnames(log2.cpm.aml2.df) <- c("geneID", sampleLabels.aml2)

###tidy data
log2.cpm.aml2.df.pivot <- pivot_longer(log2.cpm.aml2.df, 
                                       cols = WT_1:C9_3,
                                       names_to = "samples", # name of new column
                                       values_to = "expression") # name of new column storing all the data

###plot tidy data of log2 expression

p1 <- ggplot(log2.cpm.aml2.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               colour = "black", 
               show.legend = FALSE) +
  scale_fill_brewer(palette = "Set3") +
  labs(y="log2 expression", x = "Sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw(base_size = 18)


###filter data
cpm.aml2 <- cpm(myDGEList.aml2)
keepers.aml2 <- rowSums(cpm.aml2>1)>=3 #user defined - depends on studydesign
myDGEList.aml2.filtered <- myDGEList.aml2[keepers.aml2,]

log2.cpm.aml2.filtered <- cpm(myDGEList.aml2.filtered, log=TRUE)
log2.cpm.aml2.filtered.df <- as_tibble(log2.cpm.aml2.filtered, rownames = "geneID")
colnames(log2.cpm.aml2.filtered.df) <- c("geneID", sampleLabels.aml2)
log2.cpm.aml2.filtered.df.pivot <- pivot_longer(log2.cpm.aml2.filtered.df,
                                                cols = WT_1:C9_3, 
                                                names_to = "samples", # name of new column
                                                values_to = "expression") # name of new column storing all the data


p2 <- ggplot(log2.cpm.aml2.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  scale_fill_brewer(palette = "Set3") +
  labs(y="log2 expression", x = "Sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  
  theme_bw(base_size = 18)
myDGEList.aml2.filtered.norm <- calcNormFactors(myDGEList.aml2.filtered, method = "TMM")
log2.cpm.aml2.filtered.norm <- cpm(myDGEList.aml2.filtered.norm, log=TRUE)
log2.cpm.aml2.filtered.norm.df <- as_tibble(log2.cpm.aml2.filtered.norm, rownames = "geneID")
colnames(log2.cpm.aml2.filtered.norm.df) <- c("geneID", sampleLabels.aml2)
log2.cpm.aml2.filtered.norm.df.pivot <- pivot_longer(log2.cpm.aml2.filtered.norm.df,
                                                     cols = WT_1:C9_3,
                                                     names_to = "samples", # name of new column
                                                     values_to = "expression") # name of new column storing all the data


p3 <- ggplot(log2.cpm.aml2.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  scale_fill_brewer(palette = "Set3") +
  labs(y="log2 expression", x = "Sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw(base_size = 18)

plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)


### multivariate analysis ----
mydata.aml2.df <- mutate(log2.cpm.aml2.filtered.norm.df,
                         aml2.wt.AVG = (WT_1 + WT_2 + WT_3)/3, 
                         aml2.c5.AVG = (C5_1 + C5_2 + C5_3)/3,
                         aml2.c9.AVG = (C9_1 + C9_2 + C9_3)/3,
                         #now make columns comparing each of the averages above 
                         LogFC_C5_WT = (aml2.c5.AVG - aml2.wt.AVG),
                         LogFC_C9_WT = (aml2.c9.AVG - aml2.wt.AVG))%>% 
  mutate_if(is.numeric, round, 2)

groupaml2 <- targets.aml2$phenotype_sep
groupaml2 <- factor(groupaml2)

### pca aml2 ----
dark_okabe <- darken(c("#E69F00", "#56B4E9", "grey60"), amount = 0.2) 

pca.res.aml2 <- prcomp(t(log2.cpm.aml2.filtered.norm), scale.=F, retx=T)
pc.var.aml2 <- pca.res.aml2$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per.aml2 <- round(pc.var.aml2/sum(pc.var.aml2)*100, 1) 
pca.res.aml2.df <- as_tibble(pca.res.aml2$x)
pca.plot.aml2 <- ggplot(pca.res.aml2.df) +
  aes(x=PC1, y=PC2, label=sampleLabels.aml2, fill = groupaml2, colour = groupaml2, scale = TRUE) +
  geom_point(size=6, stroke = 1, shape = 21) +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per.aml2[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per.aml2[2],"%",")")) +
  #labs(title="PCA plot") +
  coord_fixed() +
  theme_classic(base_size = 22) +
  scale_fill_manual(labels = c("C5", "C9", "WT"), values = c("#E69F00", "#56B4E9", "black")) +
  scale_colour_manual(labels = c("C5", "C9", "WT"), values = dark_okabe) +
  labs(fill = "Group", colour = "Group")

pca.plot.aml2

ggsave(filename = "./RNA_Seq/plots/rna_pca.pdf", plot = pca.plot.aml2, 
       width = 17, height = 14, dpi = 800, units = "cm")  