###############################################################################
###############################################################################

# Program: Supplementary_S5.R
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
library(SMFilter)

library(reticulate)
use_python("/usr/bin/python3")
anndata <- reticulate::import("anndata")

setwd("/home/scamara/data_a/scamara/Atlas/Codigo/Rocinante_DEA/Funciones")
source("add_NTotalGenes.R")

set.seed(2504)

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### Read files #####
path <- "/home/scamara/data_a/scamara/Atlas/Input"
file <- paste0(path, "/Python_scVI_adata_big_V4_state4.h5ad")
adata <- anndata$read_h5ad(file)

sce <- AnnData2SCE(adata, "counts", uns = FALSE, obsm = FALSE, obsp = FALSE)
sce
assay(sce, "counts") %>% max()

print((sce %>% dim())[2])

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### Set PATH to save figs #####
setwd("/home/scamara/data_a/scamara/Atlas/Figuras/Suplementarias")

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### S5A - Barplot of cells per timepoint and response #####
## Filter object
filtered_sce <- sce[, colData(sce)$Antigen == "Blood"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, colData(filtered_sce)$manual_celltype_annotation_high == "CD8 cytotoxic"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, colData(filtered_sce)$STATUS != "HEALTHY"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, !((colData(filtered_sce)$Time_Point_Ranges == "Infusion_Product") & (colData(filtered_sce)$Stimulated == "YES"))]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, !((colData(filtered_sce)$Stimulated == "YES") & (colData(filtered_sce)$Stimulation_Location == "In_vitro"))]
print((filtered_sce %>% dim())[2])

# Extract metadata
meta <- as.data.frame(colData(filtered_sce))

# Reorder levels with Time_Point_Ranges variable
meta$Time_Point_Ranges <- factor(
  meta$Time_Point_Ranges,
  levels = c("Infusion_Product", "<2_weeks", "2_weeks-3_months", ">3_months"),
  ordered = TRUE
)

# Create counts table
tabla_conteos <- table(meta$Max_Response, meta$Time_Point_Ranges)
tabla_conteos

## NOTE: With this data final figure was done using Graphpad Prism 8

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### S5B - Volcanos between consecutive timepoints and S5D - Heatmap of selected genes #####
## Filter object
filtered_sce <- filtered_sce[,!is.na(filtered_sce$Max_Response)]
filtered_sce <- filtered_sce[, colData(filtered_sce)$Max_Response == "CR"]
print((filtered_sce %>% dim())[2])

colData(filtered_sce)$Time_Point_Ranges <- factor(colData(filtered_sce)$Time_Point_Ranges, levels = c("Infusion_Product", "<2_weeks", "2_weeks-3_months", ">3_months"), ordered=TRUE)
colData(filtered_sce)$Age_Range <- factor(colData(filtered_sce)$Age_Range, levels = c("<20", "20-40", "40-60", ">60"), ordered=TRUE)

colData(filtered_sce)$ScFv %>% table() # BCMA: 17863 | CD19: 59433

############################################################################################################################################################################################################

pb <- aggregateToPseudoBulk(filtered_sce,
    assay = "counts",     
    cluster_id = "Time_Point_Ranges",  
    sample_id = "Product_norm",
    BPPARAM = SnowParam(6, progressbar=TRUE))

# Evaluate the specificity of each gene for each cluster
df_cts <- cellTypeSpecificity(pb)

############################################################################################################################################################################################################
## Comparisons ##

## compare first two timepoints
ct.pairs_1 <- c("<2_weeks", "Infusion_Product")

# run comparison
fit_1 <- dreamletCompareClusters(pb, ct.pairs_1, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_1 <- topTable(fit_1, coef = "compare", number = 10)

head(res_1)

############################################################################################################################################################################################################
##### Fig S5D - compare first two timepoints summarized #####
# Rownames to filter
filas_a_filtrar <- c("DHCR7", "HMGCS1", "MSMO1")

# Filtrar el dataframe
res_1_bis <- res_1[rownames(res_1) %in% filas_a_filtrar, ]

cairo_pdf("S5D.pdf", width=10, height = 10)
dreamlet::plotHeatmap(df_cts, genes = rownames(res_1_bis))
dev.off()

############################################################################################################################################################################################################
## S5B - compare next two timepoints
ct.pairs_2 <- c("2_weeks-3_months", "<2_weeks")

# run comparison
fit_2 <- dreamletCompareClusters(pb, ct.pairs_2, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_2 <- topTable(fit_2, coef = "compare", number = 10)

head(res_2)

############################################################################################################################################################################################################
## compare last two timepoints
ct.pairs_3 <- c(">3_months", "2_weeks-3_months")

# run comparison
fit_3 <- dreamletCompareClusters(pb, ct.pairs_3, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_3 <- topTable(fit_3, coef = "compare", number = 10)

head(res_3)

############################################################################################################################################################################################################
# Make a list storing each result with a meaningful name
fitList <- list()

id <- paste0("[", ct.pairs_1[1], "]_vs_[", ct.pairs_1[2], "]")
fitList[[id]] <- fit_1

id <- paste0("[", ct.pairs_2[1], "]_vs_[", ct.pairs_2[2], "]")
fitList[[id]] <- fit_2

id <- paste0("[", ct.pairs_3[1], "]_vs_[", ct.pairs_3[2], "]")
fitList[[id]] <- fit_3

res.compare <- as.dreamletResult(fitList) # https://diseaseneurogenomics.github.io/dreamlet/reference/as.dreamletResult.html?q=dreamletCompareClusters#details
res.compare

cairo_pdf("S5B.pdf", width=10, height = 10)
plotVolcano(res.compare, coef = "compare", ncol = 4)
dev.off()

############################################################################################################################################################################################################
############################################################################################################################################################################################################
##### S5C - Volcanos between all timepoints against IP #####
pb <- aggregateToPseudoBulk(filtered_sce,
    assay = "counts",     
    cluster_id = "Time_Point_Ranges",  
    sample_id = "Product_norm",
    BPPARAM = SnowParam(6, progressbar=TRUE))

# Evaluate the specificity of each gene for each cluster
df_cts <- cellTypeSpecificity(pb)

############################################################################################################################################################################################################
## Comparisons ##
## compare first two timepoints
ct.pairs_1 <- c("<2_weeks", "Infusion_Product")

# run comparison
fit_1 <- dreamletCompareClusters(pb, ct.pairs_1, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_1 <- topTable(fit_1, coef = "compare", number = 10)

head(res_1)

############################################################################################################################################################################################################
## compare next timepoint against IP
ct.pairs_2 <- c("2_weeks-3_months", "Infusion_Product")

# run comparison
fit_2 <- dreamletCompareClusters(pb, ct.pairs_2, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_2 <- topTable(fit_2, coef = "compare", number = 10)

head(res_2)

############################################################################################################################################################################################################
## compare the last timepoint against IP
ct.pairs_3 <- c(">3_months", "Infusion_Product")

# run comparison
fit_3 <- dreamletCompareClusters(pb, ct.pairs_3, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_3 <- topTable(fit_3, coef = "compare", number = 10)

head(res_3)

############################################################################################################################################################################################################
# Make a list storing each result with a meaningful name
fitList <- list()

id <- paste0("[", ct.pairs_1[1], "]_vs_[", ct.pairs_1[2], "]")
fitList[[id]] <- fit_1

id <- paste0("[", ct.pairs_2[1], "]_vs_[", ct.pairs_2[2], "]")
fitList[[id]] <- fit_2

id <- paste0("[", ct.pairs_3[1], "]_vs_[", ct.pairs_3[2], "]")
fitList[[id]] <- fit_3

res.compare <- as.dreamletResult(fitList) # https://diseaseneurogenomics.github.io/dreamlet/reference/as.dreamletResult.html?q=dreamletCompareClusters#details
res.compare

cairo_pdf("S5C.pdf", width=10, height = 10)
plotVolcano(res.compare, coef = "compare", ncol = 4)
dev.off()

################################
######## END OF SCRIPT #########
################################