#takes the output from 3.WT_vs_clones_diff_expr.R as input

library(tidyverse)
library(ensembldb) 
library(gplots)
library(RColorBrewer)
library(GSEABase) 
library(gprofiler2) 
library(clusterProfiler) 
library(enrichplot) 
library(colorspace)
library(Cairo)

#downloaded from http://scrna.sklehabc.com/
atlas_blood_cells <- read.csv("./publicly_available_data/signatures_atlas_human_blood_cells.txt", sep = "\t") %>%
  dplyr::select(RNA_Cluster, Gene) %>%
  rename("gs_name" = "RNA_Cluster", "gene_symbol" = "Gene")