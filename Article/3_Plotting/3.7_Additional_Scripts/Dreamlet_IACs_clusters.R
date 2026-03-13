###############################################################################
###############################################################################

# Program: Dreamlet_IACs_clusters.R
# Author: Sergio Cámara Peña
# Date: 29/05/2024
# Version: FINAL

###############################################################################
###############################################################################

##### Import libraries #####
library(tidyverse)
library(cowplot)
library(patchwork)
library(dreamlet)
library(SingleCellExperiment)
library(zenith)
library(zellkonverter)
library(kableExtra)
library(scattermore)
library(EnrichmentBrowser)
library(GSEABase)
library(Cairo)
library(ggplot2)

library(reticulate)
use_python("/usr/bin/python3")
anndata <- reticulate::import("anndata")

setwd("/home/scamara/data_a/scamara/Atlas/Codigo/Rocinante_DEA/Funciones")
source("add_NTotalGenes.R")

##### Set random seed #####
set.seed(2504)

##### Read files #####
path <- "/home/scamara/data_a/scamara/Atlas/Input"
file <- paste0(path, "/Python_scVI_adata_big_V4_state4.h5ad")
adata <- anndata$read_h5ad(file)

sce <- AnnData2SCE(adata, "counts", uns = FALSE, obsm = FALSE, obsp = FALSE)
sce
assay(sce, "counts") %>% max()

print((sce %>% dim())[2])

##### Filter object #####
filtered_sce <- sce[, colData(sce)$Antigen == "Blood"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, colData(filtered_sce)$manual_celltype_annotation_high == "Monocyte-like T cells"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, colData(filtered_sce)$STATUS == "DISEASE"]
print((filtered_sce %>% dim())[2])

colData(filtered_sce)$Max_Response <- factor(colData(filtered_sce)$Max_Response, levels = c("NR", "PR", "CR"), ordered = TRUE)
colData(filtered_sce)$ICANS_Grade_Range <- factor(colData(filtered_sce)$ICANS_Grade_Range, levels = c("3-4", "1-2", "0"), ordered = TRUE)
colData(filtered_sce)$Time_Point_Ranges <- factor(colData(filtered_sce)$Time_Point_Ranges, levels = c("Infusion_Product", "<2_weeks", "2_weeks-3_months"), ordered=TRUE)
colData(filtered_sce)$Age_Range <- factor(colData(filtered_sce)$Age_Range, levels = c("<20", "20-40", "40-60", ">60"), ordered=TRUE)

# Add a new column 'IACs_differenciation' and set all values to 'Post_Infusion'
filtered_sce$IACs_differenciation <- 'Post_Infusion'

# Change values in 'IACs_differenciation' column to 'Infusion_Product' where 'Time_Point_Ranges' is 'Infusion_Product'
filtered_sce$IACs_differenciation[filtered_sce$Time_Point_Ranges == 'Infusion_Product'] <- 'Infusion_Product'

############################################################################################################################################################################################################
############################################################################################################################################################################################################

pb <- aggregateToPseudoBulk(filtered_sce,
    assay = "counts",     
    cluster_id = "IACs_differenciation",  
    sample_id = "Product_norm", 
    BPPARAM = SnowParam(6, progressbar=TRUE))

# Evaluate the specificity of each gene for each cluster
df_cts <- cellTypeSpecificity(pb)

############################################################################################################################################################################################################
############################################################################################################################################################################################################
##### Comparisons #####
############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### compare first two timepoints #####
ct.pairs_1 <- c("Post_Infusion", "Infusion_Product")

# run comparison
fit_1 <- dreamletCompareClusters(pb, ct.pairs_1, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_1 <- topTable(fit_1, coef = "compare", number = 10)
res_all <- topTable(fit_1, coef = "compare", number = Inf)

setwd("/home/scamara/data_a/scamara/Atlas/Para_Cluster_Profiler")
saveRDS(res_all, file = "dreamlet_DEG_IACs_Post_vs_Infusion.RDS")

################################
######## END OF SCRIPT #########
################################