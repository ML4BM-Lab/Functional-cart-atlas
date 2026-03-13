###############################################################################
###############################################################################

# Program: Dreamlet_IP_comparison.R
# Author: Sergio Cámara Peña
# Date: 22/01/2025
# Version: V4

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
library(scater)
library(coin)

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

filtered_sce <- filtered_sce[, colData(filtered_sce)$STATUS != "HEALTHY"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, colData(filtered_sce)$Time_Point_Ranges == "Infusion_Product"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, !((colData(filtered_sce)$Time_Point_Ranges == "Infusion_Product") & (colData(filtered_sce)$Stimulated == "YES"))]
print((filtered_sce %>% dim())[2])

colData(filtered_sce)$Age_Range <- factor(colData(filtered_sce)$Age_Range, levels = c("<20", "20-40", "40-60", ">60"), ordered=TRUE)

## Remove NAs to avoid any further problems
table(filtered_sce$Anytime_CR, useNA = "ifany")
filtered_sce <- filtered_sce[, !is.na(filtered_sce$Anytime_CR)]
table(filtered_sce$Anytime_CR, useNA = "ifany")

############################################################################################################################################################################################################
############################################################################################################################################################################################################

pb <- aggregateToPseudoBulk(filtered_sce,
    assay = "counts",     
    cluster_id = "Max_Response",  
    sample_id = "Product_norm",
    BPPARAM = SnowParam(6, progressbar=TRUE))

# Evaluate the specificity of each gene for each cluster
df_cts <- cellTypeSpecificity(pb)

############################################################################################################################################################################################################
############################################################################################################################################################################################################
##### Comparisons #####
############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### compare first two #####
ct.pairs_1 <- c("CR", "rest")

# run comparison
fit_1 <- dreamletCompareClusters(pb, ct.pairs_1, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_1 <- topTable(fit_1, coef = "compare", number = 10)

head(res_1)

setwd("/home/scamara/data_a/scamara/Atlas/Resultados/IP_comparison")
cairo_pdf("IP_Comparison_Heatmap_CR_rest.pdf", width=10, height = 10)
dreamlet::plotHeatmap(df_cts, genes = rownames(res_1))#, assays=colnames(df_cts)[2:3])
dev.off()

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### compare last two #####
ct.pairs_2 <- c("NR", "rest")

# run comparison
fit_2 <- dreamletCompareClusters(pb, ct.pairs_2, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_2 <- topTable(fit_2, coef = "compare", number = 10)

head(res_2)

setwd("/home/scamara/data_a/scamara/Atlas/Resultados/IP_comparison")
cairo_pdf("IP_Comparison_Heatmap_NR_rest.pdf", width=10, height = 10)
dreamlet::plotHeatmap(df_cts, genes = rownames(res_2))#, assays=colnames(df_cts)[2:3])
dev.off()

############################################################################################################################################################################################################
############################################################################################################################################################################################################

fitList <- list()

id <- paste0("[", ct.pairs_1[1], "]_vs_[", ct.pairs_1[2], "]")
fitList[[id]] <- fit_1

id <- paste0("[", ct.pairs_2[1], "]_vs_[", ct.pairs_2[2], "]")
fitList[[id]] <- fit_2

res.compare <- as.dreamletResult(fitList) # https://diseaseneurogenomics.github.io/dreamlet/reference/as.dreamletResult.html?q=dreamletCompareClusters#details
res.compare

cairo_pdf("IP_Comparison_Volcano.pdf", width=10, height = 10)
plotVolcano(res.compare, coef = "compare", ncol = 4)
dev.off()

############################################################################################################################################################################################################
############################################################################################################################################################################################################
##### Gene Set Enrichment Analysis #####
############################################################################################################################################################################################################
############################################################################################################################################################################################################

# GO_Biological_Process_2023
go.gs <- getGenesets(
  org = "hsa",
  db = "enrichr",
  lib = "GO_Biological_Process_2023",
  gene.id.type = "SYMBOL",
  return.type = "GeneSetCollection"
)

res.zenith <- zenith_gsa(res.compare, go.gs,
  coef = "compare",
  n_genes_min = 20
)

res.zenith$assay <- factor(res.zenith$assay, names(res.compare))

res.zenith <- add_NTotalGenes(res.zenith, go.gs)

# View the updated res.zenith
head(res.zenith)

p1 <- plotZenithResults(res.zenith, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf("IP_Comparison_Biological_Process_1.pdf", width=10, height = 10)
p1
dev.off()

data_plot <- p1$data
data_plot <- add_NTotalGenes(data_plot, go.gs)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
cairo_pdf("IP_Comparison_Biological_Process_2.pdf", width=10, height = 10)
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
go.gs <- getGenesets(
  org = "hsa",
  db = "enrichr",
  lib = "KEGG_2021_Human",
  gene.id.type = "SYMBOL",
  return.type = "GeneSetCollection"
)

res.zenith <- zenith_gsa(res.compare, go.gs,
  coef = "compare",
  n_genes_min = 20
)

res.zenith$assay <- factor(res.zenith$assay, names(res.compare))

res.zenith <- add_NTotalGenes(res.zenith, go.gs)

# View the updated res.zenith
head(res.zenith)

p2 <- plotZenithResults(res.zenith, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
pdf("IP_Comparison_KEGG_2021_Human_1.pdf")
p2
dev.off()

data_plot <- p2$data
data_plot <- add_NTotalGenes(data_plot, go.gs)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
pdf("IP_Comparison_KEGG_2021_Human_2.pdf")
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
go.gs <- getGenesets(
  org = "hsa",
  db = "enrichr",
  lib = "Reactome_2022",
  gene.id.type = "SYMBOL",
  return.type = "GeneSetCollection"
)

res.zenith <- zenith_gsa(res.compare, go.gs,
  coef = "compare",
  n_genes_min = 20
)

res.zenith$assay <- factor(res.zenith$assay, names(res.compare))

res.zenith <- add_NTotalGenes(res.zenith, go.gs)

# View the updated res.zenith
head(res.zenith)

p3 <- plotZenithResults(res.zenith, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf("IP_Comparison_Reactome_2022_1.pdf", width=10, height = 10)
p3
dev.off()

data_plot <- p3$data
data_plot <- add_NTotalGenes(data_plot, go.gs)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
cairo_pdf("IP_Comparison_Reactome_2022_2.pdf", width=10, height = 10)
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
go.gs <- getGenesets(
  org = "hsa",
  db = "enrichr",
  lib = "WikiPathway_2023_Human",
  gene.id.type = "SYMBOL",
  return.type = "GeneSetCollection"
)

res.zenith <- zenith_gsa(res.compare, go.gs,
  coef = "compare",
  n_genes_min = 20
)

res.zenith$assay <- factor(res.zenith$assay, names(res.compare))

res.zenith <- add_NTotalGenes(res.zenith, go.gs)

# View the updated res.zenith
head(res.zenith)

p4 <- plotZenithResults(res.zenith, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf("IP_Comparison_WikiPathway_2023_Human_1.pdf", width=10, height = 10)
p4
dev.off()

data_plot <- p4$data
data_plot <- add_NTotalGenes(data_plot, go.gs)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
cairo_pdf("IP_Comparison_WikiPathway_2023_Human_2.pdf", width=10, height = 10)
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

res.zenith = zenith_gsa( res.compare, gsc, 
              coef = 'compare', 
              n_genes_min=20)

res.zenith$assay = factor(res.zenith$assay, names(res.compare))

res.zenith <- add_NTotalGenes(res.zenith, gsc)

# View the updated res.zenith
head(res.zenith)

p5 <- plotZenithResults(res.zenith, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
setwd("/home/scamara/data_a/scamara/Atlas/Resultados/IP_comparison")
pdf("IP_Comparison_Custom_Genesets_1.pdf")
p5
dev.off()

data_plot <- p5$data
data_plot <- add_NTotalGenes(data_plot, gsc)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
pdf("IP_Comparison_Custom_Genesets_2.pdf")
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

## Filter KEGG_2021_Human
KEGG_2021_Human$Geneset

# Filter by
patrones <- c("ER69_", "ER32_")

# Put all in a variable with or
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
KEGG_2021_Human_filtered <- KEGG_2021_Human[grep(patron_combinado, KEGG_2021_Human$Geneset), ]

###################################################################################################################################################################################################

## Filter Reactome_2022
Reactome_2022$Geneset

# Filter by
patrones <- c("ER755_", "ER713_")

# Put all in a variable with or
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
Reactome_2022_filtered <- Reactome_2022[grep(patron_combinado, Reactome_2022$Geneset), ]

###################################################################################################################################################################################################

## Filter WikiPathway_2023_Human
WikiPathway_2023_Human$Geneset

# Filter by
patrones <- c("ER112_", "ER633_")

# Put all in a variable with or
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
WikiPathway_2023_Human_filtered <- WikiPathway_2023_Human[grep(patron_combinado, WikiPathway_2023_Human$Geneset), ]

###################################################################################################################################################################################################

## Filter Custom_genesets
Custom_genesets$Geneset

# Filter by
patrones <- c("Genes Prolif")

# Put all in a variable with or
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
Custom_genesets_filtered <- Custom_genesets[grep(patron_combinado, Custom_genesets$Geneset), ]

###################################################################################################################################################################################################

## Final filtered
Final_result <- rbind(KEGG_2021_Human_filtered, Reactome_2022_filtered, WikiPathway_2023_Human_filtered, Custom_genesets_filtered)

head(Final_result)

saveRDS(Final_result, "Final_filtered.RDS")

p6 <- plotZenithResults(Final_result, Inf, Inf, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
setwd("/home/scamara/data_a/scamara/Atlas/Resultados/IP_comparison")
cairo_pdf("IP_Comparison_Final_Summarized_Result_1.pdf", width=10, height = 10)
p6
dev.off()

data_plot <- p6$data
data_plot <- add_NTotalGenes(data_plot, gsc)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
cairo_pdf("IP_Comparison_Final_Summarized_Result_2.pdf", width=10, height = 10)
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

## Do these cells express more IL-10?
max(assay(filtered_sce, "counts")["IL10", ], na.rm = TRUE)
median(assay(filtered_sce, "counts")["IL10", ], na.rm = TRUE)

sum(assay(filtered_sce, "counts")["IL10", ] > 0, na.rm = TRUE)
dim(filtered_sce)[2]
sum(assay(filtered_sce, "counts")["IL10", ] > 0, na.rm = TRUE)/dim(filtered_sce)[2]

# Create the contingency table
contingency_table <- table(
  assay(filtered_sce, "counts")["IL10", ] > 0,
  filtered_sce$Max_Response
)
print(contingency_table)

contingency_table <- contingency_table[, c("CR", "NR")]
print(contingency_table)

## Perform Fisher's Exact Test - The best one
fisher.test(contingency_table)

## Perform Chi-squared test
chisq.test(contingency_table[, c("CR", "NR")])

## Permutation test
# Subset data to exclude "PR" responses from Max_Response
filtered_data <- data.frame(
  IL10_expr = assay(filtered_sce, "counts")["IL10", ] > 0,
  Max_Response = filtered_sce$Max_Response
)

# Keep only "CR" and "NR" groups
filtered_data <- filtered_data[filtered_data$Max_Response %in% c("CR", "NR"), ]

# Perform the permutation test
perm_test <- independence_test(IL10_expr ~ Max_Response, data = filtered_data, distribution = approximate(nresample = 1000))

# View the result
perm_test # Significant

# Pie chart - CR vs NR
filtered_data_2 <- filtered_data[filtered_data$IL10_expr == "TRUE", ]
tabla <- table(filtered_data_2$Max_Response)
tabla_filtrada <- tabla[tabla > 0]
pie(tabla_filtrada, main = "", col = c("forestgreen", "tomato"))
dev.off()

# Pie chart - Il-10+ vs negative
filtered_data <- data.frame(
  IL10_expr = assay(filtered_sce, "counts")["IL10", ] > 0,
  Max_Response = filtered_sce$Max_Response
)

tabla_2 <- table(filtered_data$IL10_expr)
names(tabla_2) <- c("IL-10-", "IL-10+")

# Colores para las categorías
colores <- c("tomato", "forestgreen")

# Crear el gráfico sin etiquetas (labels = NA)
pie(tabla_2, labels = NA, col = colores, main = "")

# Añadir leyenda fuera del gráfico
legend("topright", legend = names(tabla_2), fill = colores, cex = 2.8, bty = "n")

dev.off()

################################
######## END OF SCRIPT #########
################################
