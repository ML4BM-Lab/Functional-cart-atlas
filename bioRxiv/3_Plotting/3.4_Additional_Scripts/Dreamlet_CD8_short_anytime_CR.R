###############################################################################
###############################################################################

# Program: Dreamlet_CD8_short_anytime_CR.R
# Author: Sergio Cámara Peña
# Date: 04/12/2024
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

filtered_sce <- filtered_sce[, colData(filtered_sce)$manual_celltype_annotation_high == "CD8 cytotoxic"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, colData(filtered_sce)$Time_Point_Ranges != "2_weeks-3_months"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, colData(filtered_sce)$Time_Point_Ranges != ">3_months"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, colData(filtered_sce)$STATUS != "HEALTHY"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, !((colData(filtered_sce)$Time_Point_Ranges == "Infusion_Product") & (colData(filtered_sce)$Stimulated == "YES"))]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, !((colData(filtered_sce)$Stimulated == "YES") & (colData(filtered_sce)$Stimulation_Location == "In_vitro"))]
print((filtered_sce %>% dim())[2])

colData(filtered_sce)$Max_Response %>% table(useNA = "ifany")
filtered_sce <- filtered_sce[, complete.cases(colData(filtered_sce)$Max_Response)]
print((filtered_sce %>% dim())[2])
colData(filtered_sce)$Max_Response %>% table(useNA = "ifany")

# Reorder the levels 
colData(filtered_sce)$Time_Point_Ranges <- factor(colData(filtered_sce)$Time_Point_Ranges, levels = c("Infusion_Product", "<2_weeks"), ordered = TRUE)
colData(filtered_sce)$Age_Range <- factor(colData(filtered_sce)$Age_Range, levels = c("<20", "20-40", "40-60", ">60"), ordered=TRUE)
colData(filtered_sce)$Max_Response <- factor(colData(filtered_sce)$Max_Response, levels = c("NR", "PR", "CR"), ordered=TRUE)

pb <- aggregateToPseudoBulk(filtered_sce,
    assay = "counts",     
    cluster_id = "Anytime_CR",  
    sample_id = "Product_norm", 
    BPPARAM = SnowParam(6, progressbar=TRUE))

colData(pb) %>% colnames()

form <- ~ Time_Point_Ranges + Age_Range + Sex
C = canCorPairs(form, colData(pb))

plotCorrMatrix(C)
dev.off()

plotCorrMatrix(C)
dev.off()

form = ~ Time_Point_Ranges + (1|Age_Range) + (1|Sex)

res.proc = processAssays(pb, form, 
  min.samples = 5,
  BPPARAM=SnowParam(12, progressbar=TRUE))

details(res.proc)

vp.lst = fitVarPart( res.proc, form, 
  BPPARAM=SnowParam(12, progressbar=TRUE)) 


plotVarPart(vp.lst, label.angle=60, ncol=4)   
dev.off()


res.dl = dreamlet(res.proc, form,
    BPPARAM=SnowParam(12, progressbar=TRUE))

res.dl

plotVolcano( res.dl, coef = 'Time_Point_Ranges.L', ncol=4) 
dev.off()

genes_a = c("ARAP2", "SIK3", "ZXDC", "AXIN1", "CBX4", "SMIM3", "STXBP1", "TRIP10", "LINC01943", "MSMO1")
plotGeneHeatmap(res.dl, coef = "Time_Point_Ranges.L", genes = genes_a)
dev.off()

topTable(res.dl, coef = "Time_Point_Ranges.L")

###################################################################################################################################################################################################
###################################################################################################################################################################################################

# GO_Biological_Process_2023
go.gs = getGenesets( org = "hsa", 
                  db = "enrichr", 
                  lib = "GO_Biological_Process_2023",
                  gene.id.type = "SYMBOL", 
                  return.type = "GeneSetCollection")

res.zenith = zenith_gsa( res.dl, go.gs, 
              coef = 'Time_Point_Ranges.L', 
              n_genes_min=20)

res.zenith$assay = factor(res.zenith$assay, names(res.dl))

res.zenith <- add_NTotalGenes(res.zenith, go.gs)

# View the updated res.zenith
head(res.zenith)

p1 <- plotZenithResults(res.zenith, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
p1
dev.off()

data_plot <- p1$data
data_plot <- add_NTotalGenes(data_plot, go.gs)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
ggplot(data_plot, aes(x = assay, y = Geneset, size=logFDR, color = delta)) +
  geom_point() +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Delta") +
  labs(
    x = "Comparison",
    y = "GO term",
    title = ""
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
  )
dev.off()

GO_Biological_Process_2023 <- res.zenith

###################################################################################################################################################################################################
###################################################################################################################################################################################################

# KEGG_2021_Human
go.gs = getGenesets( org = "hsa", 
                  db = "enrichr", 
                  lib = "KEGG_2021_Human",
                  gene.id.type = "SYMBOL", 
                  return.type = "GeneSetCollection")

res.zenith = zenith_gsa( res.dl, go.gs, 
              coef = 'Time_Point_Ranges.L', 
              n_genes_min=20)

res.zenith$assay = factor(res.zenith$assay, names(res.dl))

res.zenith <- add_NTotalGenes(res.zenith, go.gs)

# View the updated res.zenith
head(res.zenith)

p2 <- plotZenithResults(res.zenith, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
p2
dev.off()

data_plot <- p2$data
data_plot <- add_NTotalGenes(data_plot, go.gs)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
ggplot(data_plot, aes(x = assay, y = Geneset, size=logFDR, color = delta)) +
  geom_point() +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Delta") +
  labs(
    x = "Comparison",
    y = "GO term",
    title = ""
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
  )
dev.off()

KEGG_2021_Human <- res.zenith

###################################################################################################################################################################################################
###################################################################################################################################################################################################

# Reactome_2022
go.gs = getGenesets( org = "hsa", 
                  db = "enrichr", 
                  lib = "Reactome_2022",
                  gene.id.type = "SYMBOL", 
                  return.type = "GeneSetCollection")

res.zenith = zenith_gsa( res.dl, go.gs, 
              coef = 'Time_Point_Ranges.L', 
              n_genes_min=20)

res.zenith$assay = factor(res.zenith$assay, names(res.dl))

res.zenith <- add_NTotalGenes(res.zenith, go.gs)

# View the updated res.zenith
head(res.zenith)

p3 <- plotZenithResults(res.zenith, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
p3
dev.off()

data_plot <- p3$data
data_plot <- add_NTotalGenes(data_plot, go.gs)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
ggplot(data_plot, aes(x = assay, y = Geneset, size=logFDR, color = delta)) +
  geom_point() +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Delta") +
  labs(
    x = "Comparison",
    y = "GO term",
    title = ""
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
  )
dev.off()

Reactome_2022 <- res.zenith

###################################################################################################################################################################################################
###################################################################################################################################################################################################

# WikiPathway_2023_Human
go.gs = getGenesets( org = "hsa", 
                  db = "enrichr", 
                  lib = "WikiPathway_2023_Human",
                  gene.id.type = "SYMBOL", 
                  return.type = "GeneSetCollection")

res.zenith = zenith_gsa( res.dl, go.gs, 
              coef = 'Time_Point_Ranges.L', 
              n_genes_min=20)

res.zenith$assay = factor(res.zenith$assay, names(res.dl))

res.zenith <- add_NTotalGenes(res.zenith, go.gs)

# View the updated res.zenith
head(res.zenith)

p4 <- plotZenithResults(res.zenith, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
p4
dev.off()

data_plot <- p4$data
data_plot <- add_NTotalGenes(data_plot, go.gs)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
ggplot(data_plot, aes(x = assay, y = Geneset, size=logFDR, color = delta)) +
  geom_point() +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Delta") +
  labs(
    x = "Comparison",
    y = "GO term",
    title = ""
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
  )
dev.off()

WikiPathway_2023_Human <- res.zenith

###################################################################################################################################################################################################
###################################################################################################################################################################################################
# Custom genesets
setwd("/home/scamara/data_a/scamara/Atlas/Input")
Firmas_atlas_df <- read.csv("Firmas_Atlas.csv", header = TRUE, sep = ";")
Firmas_atlas_df %>% colnames()

Firma_Activation <- Firmas_atlas_df$Activation
Firma_Activation[Firma_Activation == ""] <- NA
Firma_Activation <- Firma_Activation[!is.na(Firma_Activation)]
Firma_Activation

Firma_Tonic <- Firmas_atlas_df$Tonic
Firma_Tonic[Firma_Tonic == ""] <- NA
Firma_Tonic <- Firma_Tonic[!is.na(Firma_Tonic)]
Firma_Tonic

Firma_Genes_Polyfun <- Firmas_atlas_df$Genes_Polyfun
Firma_Genes_Polyfun[Firma_Genes_Polyfun == ""] <- NA
Firma_Genes_Polyfun <- Firma_Genes_Polyfun[!is.na(Firma_Genes_Polyfun)]
Firma_Genes_Polyfun

Firma_Genes_Prolif <- Firmas_atlas_df$Genes_Prolif
Firma_Genes_Prolif[Firma_Genes_Prolif == ""] <- NA
Firma_Genes_Prolif <- Firma_Genes_Prolif[!is.na(Firma_Genes_Prolif)]
Firma_Genes_Prolif

Firma_Genes_Teff <- Firmas_atlas_df$Genes_Teff
Firma_Genes_Teff[Firma_Genes_Teff == ""] <- NA
Firma_Genes_Teff <- Firma_Genes_Teff[!is.na(Firma_Genes_Teff)]
Firma_Genes_Teff

Firma_GlucoliticProcess <- Firmas_atlas_df$GlucoliticProcess
Firma_GlucoliticProcess[Firma_GlucoliticProcess == ""] <- NA
Firma_GlucoliticProcess <- Firma_GlucoliticProcess[!is.na(Firma_GlucoliticProcess)]
Orig_Firma_GlucoliticProcess <- Firma_GlucoliticProcess
Firma_GlucoliticProcess <- unique(Firma_GlucoliticProcess)
Firma_GlucoliticProcess

Firma_CD8_CARHigh_Signature <- Firmas_atlas_df$CD8_CARHigh_Signature
Firma_CD8_CARHigh_Signature[Firma_CD8_CARHigh_Signature == ""] <- NA
Firma_CD8_CARHigh_Signature <- Firma_CD8_CARHigh_Signature[!is.na(Firma_CD8_CARHigh_Signature)]
Firma_CD8_CARHigh_Signature

Firma_KEGG1 <- Firmas_atlas_df$KEGG1
Firma_KEGG1[Firma_KEGG1 == ""] <- NA
Firma_KEGG1 <- Firma_KEGG1[!is.na(Firma_KEGG1)]
Firma_KEGG1

Firma_KEGG2 <- Firmas_atlas_df$KEGG2
Firma_KEGG2[Firma_KEGG2 == ""] <- NA
Firma_KEGG2 <- Firma_KEGG2[!is.na(Firma_KEGG2)]
Firma_KEGG2

Firma_Exhaustion1 <- Firmas_atlas_df$Exhaustion1
Firma_Exhaustion1[Firma_Exhaustion1 == ""] <- NA
Firma_Exhaustion1 <- Firma_Exhaustion1[!is.na(Firma_Exhaustion1)]
Orig_Firma_Exhaustion1 <- Firma_Exhaustion1
Firma_Exhaustion1 <- unique(Firma_Exhaustion1)
Firma_Exhaustion1

Firma_Exhau_2 <- Firmas_atlas_df$Exhau_2
Firma_Exhau_2[Firma_Exhau_2 == ""] <- NA
Firma_Exhau_2 <- Firma_Exhau_2[!is.na(Firma_Exhau_2)]
Firma_Exhau_2

Firma_Temeff_TnTcm <- Firmas_atlas_df$Temeff_TnTcm
Firma_Temeff_TnTcm[Firma_Temeff_TnTcm == ""] <- NA
Firma_Temeff_TnTcm <- Firma_Temeff_TnTcm[!is.na(Firma_Temeff_TnTcm)]
Firma_Temeff_TnTcm

Firma_Gluconeogenesis <- Firmas_atlas_df$Gluconeogenesis
Firma_Gluconeogenesis[Firma_Gluconeogenesis == ""] <- NA
Firma_Gluconeogenesis <- Firma_Gluconeogenesis[!is.na(Firma_Gluconeogenesis)]
Orig_Firma_Gluconeogenesis <- Firma_Gluconeogenesis
Firma_Gluconeogenesis <- unique(Firma_Gluconeogenesis)
Firma_Gluconeogenesis

Firma_TricarboxylicAcidCycle <- Firmas_atlas_df$TricarboxylicAcidCycle
Firma_TricarboxylicAcidCycle[Firma_TricarboxylicAcidCycle == ""] <- NA
Firma_TricarboxylicAcidCycle <- Firma_TricarboxylicAcidCycle[!is.na(Firma_TricarboxylicAcidCycle)]
Orig_Firma_TricarboxylicAcidCycle <- Firma_TricarboxylicAcidCycle
Firma_TricarboxylicAcidCycle <- unique(Firma_TricarboxylicAcidCycle)
Firma_TricarboxylicAcidCycle

gs1 <- GeneSet(setName="Activation", geneIds=Firma_Activation)
gs2 <- GeneSet(setName="Tonic", geneIds=Firma_Tonic)
gs3 <- GeneSet(setName="Genes Polyfun", geneIds=Firma_Genes_Polyfun)
gs4 <- GeneSet(setName="Genes Prolif", geneIds=Firma_Genes_Prolif)
gs5 <- GeneSet(setName="Genes Teff", geneIds=Firma_Genes_Teff)
gs6 <- GeneSet(setName="Glucolitic Process", geneIds=Firma_GlucoliticProcess)
gs7 <- GeneSet(setName="CD8 CARHigh Signature", geneIds=Firma_CD8_CARHigh_Signature)
gs8 <- GeneSet(setName="Oxidative Phosphorylation", geneIds=Firma_KEGG1)
gs9 <- GeneSet(setName="Citrate Cycle_TCA Cycle ", geneIds=Firma_KEGG2)
gs10 <- GeneSet(setName="Exhaustion1", geneIds=Firma_Exhaustion1)
gs11 <- GeneSet(setName="Exhaustion2", geneIds=Firma_Exhau_2)
gs12 <- GeneSet(setName="Temeff TnTcm", geneIds=Firma_Temeff_TnTcm)
gs13 <- GeneSet(setName="Gluconeogenesis", geneIds=Firma_Gluconeogenesis)
gs14 <- GeneSet(setName="TricarboxylicAcidCycle", geneIds=Firma_TricarboxylicAcidCycle)

gsc <- GeneSetCollection(gs1, gs2, gs3, gs4, gs5, gs6, gs7, gs8, gs9, gs10, gs11, gs12, gs13, gs14)

res.zenith = zenith_gsa( res.dl, gsc, 
              coef = 'Time_Point_Ranges.L', 
              n_genes_min=20)

res.zenith$assay = factor(res.zenith$assay, names(res.dl))

res.zenith <- add_NTotalGenes(res.zenith, gsc)

# View the updated res.zenith
head(res.zenith)

p5 <- plotZenithResults(res.zenith, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
p5
dev.off()

data_plot <- p5$data
data_plot <- add_NTotalGenes(data_plot, gsc)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
ggplot(data_plot, aes(x = assay, y = Geneset, size=logFDR, color = delta)) +
  geom_point() +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Delta") +
  labs(
    x = "Comparison",
    y = "GO term",
    title = ""
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
  )
dev.off()

Custom_genesets <- res.zenith

###################################################################################################################################################################################################
###################################################################################################################################################################################################
##### Filtered_Results #####

## Filter GO_Biological_Process_2023
GO_Biological_Process_2023$Geneset

# Filter by
patrones <- c("ER4809_")

# Put all in a variable with or
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
GO_Biological_Process_2023_filtered <- GO_Biological_Process_2023[grep(patron_combinado, GO_Biological_Process_2023$Geneset), ]

###################################################################################################################################################################################################

## Filter WikiPathway_2023_Human
WikiPathway_2023_Human$Geneset

# Filter by
patrones <- c("ER224_")

# Put all in a variable with or
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
WikiPathway_2023_Human_filtered <- WikiPathway_2023_Human[grep(patron_combinado, WikiPathway_2023_Human$Geneset), ]

###################################################################################################################################################################################################

## Filter Custom_genesets
Custom_genesets$Geneset

# Filter by
patrones <- c("Oxidative Phosphorylation")

# Put all in a variable with or
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
Custom_genesets_filtered <- Custom_genesets[grep(patron_combinado, Custom_genesets$Geneset), ]

###################################################################################################################################################################################################

## Final filtered
Final_result <- rbind(GO_Biological_Process_2023_filtered, WikiPathway_2023_Human_filtered, Custom_genesets_filtered)

head(Final_result)

setwd("/home/scamara/data_a/scamara/Atlas/Resultados/CD8_Short_Anytime_CR")
if(TRUE){
  write.csv(Final_result, "Final_result_Fig_2C.csv", row.names = FALSE)
} else{
  Final_result <- read.csv("Final_result_Fig_2C.csv")
}

p6 <- plotZenithResults(Final_result, Inf, Inf, sortByGeneset = FALSE)

###################################################################################################################################################################################################
###################################################################################################################################################################################################
###################################################################################################################################################################################################

##### For suplementary S7A #####
## Filter GO_Biological_Process_2023
GO_Biological_Process_2023$Geneset

# Filter by
patrones <- c("ER4809_", "ER786_", "ER781_", "ER5286_", "ER3802_", "ER2651_", "ER231_", "ER1707_")

# Put all in a variable with or
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
GO_Biological_Process_2023_filtered <- GO_Biological_Process_2023[grep(patron_combinado, GO_Biological_Process_2023$Geneset), ]

###################################################################################################################################################################################################

## Filter WikiPathway_2023_Human
WikiPathway_2023_Human$Geneset

# Filter by
patrones <- c("ER721_", "ER49_", "ER473_", "ER278_", "ER224_")

# Put all in a variable with or
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
WikiPathway_2023_Human_filtered <- WikiPathway_2023_Human[grep(patron_combinado, WikiPathway_2023_Human$Geneset), ]

###################################################################################################################################################################################################

## Filter Reactome_2022
Reactome_2022$Geneset

# Filter by
patrones <- c("ER909_", "ER74_", "ER712_", "ER1584_", "ER1531_")

# Put all in a variable with or
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
Reactome_2022_filtered <- Reactome_2022[grep(patron_combinado, Reactome_2022$Geneset), ]

###################################################################################################################################################################################################

## Filter KEGG_2021_Human
KEGG_2021_Human$Geneset

# Filter by
patrones <- c("ER198_")

# Put all in a variable with or
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
KEGG_2021_Human_filtered <- KEGG_2021_Human[grep(patron_combinado, KEGG_2021_Human$Geneset), ]

###################################################################################################################################################################################################

## Filter Custom_genesets
Custom_genesets$Geneset

# Filter by
patrones <- c("Oxidative Phosphorylation")

# Put all in a variable with or
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
Custom_genesets_filtered <- Custom_genesets[grep(patron_combinado, Custom_genesets$Geneset), ]

###################################################################################################################################################################################################

## Final filtered
Final_result <- rbind(GO_Biological_Process_2023_filtered, WikiPathway_2023_Human_filtered, Reactome_2022_filtered, KEGG_2021_Human_filtered, Custom_genesets_filtered)

head(Final_result)
setwd("/home/scamara/data_a/scamara/Atlas/Input")
saveRDS(Final_result, "Final_result_S7A.RDS")

p7 <- plotZenithResults(Final_result, Inf, Inf, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic

p7
dev.off()

################################
######## END OF SCRIPT #########
################################