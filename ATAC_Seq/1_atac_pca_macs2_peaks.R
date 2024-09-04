library(tidyverse)
library(colorspace)
library(Cairo) 

pca <- read.csv("./ATAC_Seq/bwa/merged_library/macs2/narrow_peak/consensus/deseq2/consensus_peaks.mLb.clN.pca.vals.txt", sep = "\t", comment.char = "#", check.names = FALSE) %>%
  mutate(Group = gsub('_.*','', sample))

dark_okabe <- darken(c("#E69F00", "#56B4E9", "grey60"), amount = 0.2) 

atac_pca <- ggplot(pca) +
  aes(x=`PC1: 44% variance`, y=`PC2: 22% variance`, colour = Group, fill = Group) +
  geom_point(size=6, stroke = 1, shape = 21) +
  xlab("PC1 (44%)") + 
  ylab("PC2 (22%)") +
  coord_fixed() +
  expand_limits(y = c(-1.2, 1.8)) +
  theme_classic(base_size = 22) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#000000")) + 
  scale_colour_manual(values = dark_okabe)

ggsave(filename = "./plots/atac_pca.pdf", plot = atac_pca, 
       width = 17, height = 14, dpi = 800, units = "cm", device = cairo_pdf)