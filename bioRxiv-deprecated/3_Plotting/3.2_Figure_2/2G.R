###############################################################################
###############################################################################

# Program: 2G.R
# Author: Sergio Cámara Peña
# Date: 09/07/2025
# Version: FINAL

###############################################################################
###############################################################################

##### Load required libraries #####
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

##### Set seed #####
set.seed(2504)

##### Load data #####
setwd("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Dreamlet_to_ClusterProfiler")
res_all <- readRDS("dreamlet_DEG_IACs_Post_vs_Infusion.RDS") # This object comes from Dreamlet_IACs_clusters.R script
res_all %>% head()

# Create a named vector with logFC values and gene names
geneList <- res_all$logFC
names(geneList) <- rownames(res_all)
geneList <- sort(geneList, decreasing = TRUE)
sig_genes <- rownames(res_all[res_all$adj.P.Val < 0.05, ])

gene.df <- bitr(sig_genes, fromType = "SYMBOL",
        toType = c("ENTREZID"),
        OrgDb = org.Hs.eg.db)

# Filter and rebuild the geneList with ENTREZID
geneList_entrez <- geneList[gene.df$SYMBOL]
names(geneList_entrez) <- gene.df$ENTREZID
geneList_entrez <- sort(geneList_entrez, decreasing = TRUE)

gsea_res <- gseGO(geneList     = geneList_entrez,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "BP",
                  keyType      = "ENTREZID",
                  minGSSize    = 10,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)

setwd("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Figura_2")

p = dotplot(gsea_res, showCategory = 8) + 
  ggtitle("GO Biological Process")

ggsave("2G.pdf", plot = p, width = 10, height = 6)

##### END OF SCRIPT #####
