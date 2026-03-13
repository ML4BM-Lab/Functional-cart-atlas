###############################################################################
###############################################################################

# Program: Dreamlet_V5_BCMA_vs_CD19_IP_All.R
# Author: Sergio Cámara Peña
# Date: 27/05/2025
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

set.seed(2504)

##### Read files #####
path <- "/home/scamara/data_a/scamara/Atlas/Input"
file <- paste0(path, "/Atlas_integ_scArches_FINAL_V5.h5ad")
adata <- anndata$read_h5ad(file)

sce <- AnnData2SCE(adata, "counts", uns = FALSE, obsm = FALSE, obsp = FALSE)
sce
assay(sce, "counts") %>% max()

print((sce %>% dim())[2])

##### Filter object #####
filtered_sce <- sce[, colData(sce)$Antigen == "Blood"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, !((colData(filtered_sce)$Time_Point_Ranges == "Infusion_Product") & (colData(filtered_sce)$Stimulated == "YES"))]
print((filtered_sce %>% dim())[2])

filtered_sce <- sce[, colData(sce)$STATUS == "DISEASE"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, !((colData(filtered_sce)$Stimulated == "YES") & (colData(filtered_sce)$Stimulation_Location == "In_vitro"))]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, !is.na(filtered_sce$Max_Response)]
filtered_sce <- filtered_sce[, colData(filtered_sce)$Max_Response == "CR"]
print((filtered_sce %>% dim())[2])

filtered_sce_2 <- filtered_sce

filtered_sce <- filtered_sce[, (colData(filtered_sce)$Time_Point_Ranges == "Infusion_Product")]
print((filtered_sce %>% dim())[2])

table(filtered_sce$ScFv)

colData(filtered_sce)$ScFv %>% table()
colData(filtered_sce)$Source %>% table()
colData(filtered_sce)$Time_Point_Ranges %>% table()

colData(filtered_sce)$Age_Range <- factor(colData(filtered_sce)$Age_Range, levels = c("<20", "20-40", "40-60", ">60"), ordered = TRUE)

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### BARPLOT WITH DATA - GENERAL #####
# Convert metadata to a data frame
meta_df <- as.data.frame(colData(filtered_sce_2))

# Ensure relevant columns exist
if (!all(c("ScFv", "Time_Point_Ranges", "Product_norm") %in% colnames(meta_df))) {
  stop("Missing required columns in colData(filtered_sce)")
}

# Remove duplicates based on ScFv, Time_Point_Ranges, and Product_norm
unique_df <- meta_df %>%
  distinct(ScFv, Time_Point_Ranges, Product_norm, .keep_all = TRUE)

# Count the unique occurrences
counts_df <- unique_df %>%
  group_by(ScFv, Time_Point_Ranges) %>%
  summarise(Sample_Count = n(), .groups = "drop")

# Convert categorical variables to factors
counts_df$ScFv <- as.factor(counts_df$ScFv)
counts_df$Time_Point_Ranges <- factor(counts_df$Time_Point_Ranges,
  levels = c(
    "Infusion_Product", "<2_weeks",
    "2_weeks-3_months", ">3_months"
  )
)

ggplot(counts_df, aes(x = Time_Point_Ranges, y = Sample_Count, fill = ScFv)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Sample_Count), vjust = -0.5, position = position_dodge(0.9)) +
  labs(
    title = "Unique Sample Counts by Time Point Ranges and ScFv",
    x = "Time Point Ranges",
    y = "Sample Count"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

## Sample distribution
counts_df_2 <- unique_df %>%
  group_by(Norm_Patient_Name, Product_norm, Time_Point_Ranges, ScFv) %>%
  summarise(.groups = "drop")

dim(counts_df_2)

counts_df_2 <- counts_df_2[counts_df_2$Time_Point_Ranges %in% c("Infusion_Product", "2_weeks-3_months"), ]
dim(counts_df_2)
print(counts_df_2, n = 100)

############################################################################################################################################################################################################
############################################################################################################################################################################################################

pb <- aggregateToPseudoBulk(filtered_sce,
  assay = "counts",
  cluster_id = "ScFv",
  sample_id = "Product_norm",
  BPPARAM = SnowParam(6, progressbar = TRUE)
)

# Evaluate the specificity of each gene for each cluster
df_cts <- cellTypeSpecificity(pb)

############################################################################################################################################################################################################
############################################################################################################################################################################################################
##### Comparisons #####
############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### compare both scFv #####
ct.pairs_1 <- c("BCMA", "CD19")

# run comparison
fit_1 <- dreamletCompareClusters(pb, ct.pairs_1, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_1 <- topTable(fit_1, coef = "compare", number = 10)

res_all <- topTable(fit_1, coef = "compare", number = Inf)

setwd("/home/scamara/data_a/scamara/Atlas/Para_Cluster_Profiler")
saveRDS(res_all, file = "Resultados_V5_BCMA_vs_CD19_IP_All.RDS")

head(res_1)

setwd("/home/scamara/data_a/scamara/Atlas/Resultados_V5/BCMA_vs_CD19_IP_All")
cairo_pdf("BCMA_vs_CD19_IP_All_Heatmap.pdf", width = 10, height = 10)
dreamlet::plotHeatmap(df_cts, genes = rownames(res_1)) # , assays=colnames(df_cts)[2:3])
dev.off()

############################################################################################################################################################################################################
############################################################################################################################################################################################################

# Make a list storing each result with a meaningful name
fitList <- list()

id <- paste0("[", ct.pairs_1[1], "]_vs_[", ct.pairs_1[2], "]")
fitList[[id]] <- fit_1

res.compare <- as.dreamletResult(fitList) # https://diseaseneurogenomics.github.io/dreamlet/reference/as.dreamletResult.html?q=dreamletCompareClusters#details
res.compare

cairo_pdf("BCMA_vs_CD19_IP_All_Volcano.pdf", width = 10, height = 10)
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

# Plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf("BCMA_vs_CD19_IP_All_Biological_Process_1.pdf", width = 10, height = 10)
plotZenithResults(res.zenith)
dev.off()

gs <- unique(res.zenith$Geneset[res.zenith$FDR < 0.05])

# keep only results of these genesets
df <- res.zenith[res.zenith$Geneset %in% gs, ]

# plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf("BCMA_vs_CD19_IP_All_Biological_Process_2.pdf", width = 10, height = 10)
plotZenithResults(df, Inf, Inf)
dev.off()

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

# plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf("BCMA_vs_CD19_IP_All_KEGG_2021_Human_1.pdf", width = 10, height = 10)
plotZenithResults(res.zenith)
dev.off()

gs <- unique(res.zenith$Geneset[res.zenith$FDR < 0.05])

# keep only results of these genesets
df <- res.zenith[res.zenith$Geneset %in% gs, ]

# plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf("BCMA_vs_CD19_IP_All_KEGG_2021_Human_2.pdf", width = 10, height = 10)
plotZenithResults(df, Inf, Inf)
dev.off()

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

# plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf("BCMA_vs_CD19_IP_All_Reactome_2022_1.pdf", width = 10, height = 10)
plotZenithResults(res.zenith)
dev.off()

gs <- unique(res.zenith$Geneset[res.zenith$FDR < 0.05])

# keep only results of these genesets
df <- res.zenith[res.zenith$Geneset %in% gs, ]

# plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf("BCMA_vs_CD19_IP_All_Reactome_2022_2.pdf", width = 10, height = 10)
plotZenithResults(df, Inf, Inf)
dev.off()

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

# plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf("BCMA_vs_CD19_IP_All_WikiPathway_2023_Human_1.pdf", width = 10, height = 10)
plotZenithResults(res.zenith)
dev.off()

gs <- unique(res.zenith$Geneset[res.zenith$FDR < 0.05])

# keep only results of these genesets
df <- res.zenith[res.zenith$Geneset %in% gs, ]

# plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf("BCMA_vs_CD19_IP_All_WikiPathway_2023_Human_2.pdf", width = 10, height = 10)
plotZenithResults(df, Inf, Inf)
dev.off()

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

Firma_Exhaustion_Summarized <- c("TIGIT", "CTLA4", "HAVCR2", "LAG3", "SPATA2") # HAVCR2 = TIM3; SPATA2 = PD1

Firma_Cytotoxicity <- c("PRF1", "IFNG", "NKG7", "GNLY", "GZMA", "GZMK", "GZMB", "GZMH", "CD7", "CCL5", "CD8A", "CCL4")

Firma_Memory <- c("TCF7", "CCR7", "SELL", "IL7R")

gs1 <- GeneSet(setName = "Firma_Activation", geneIds = Firma_Activation)
gs2 <- GeneSet(setName = "Firma_Tonic", geneIds = Firma_Tonic)
gs3 <- GeneSet(setName = "Firma_Genes_Polyfun", geneIds = Firma_Genes_Polyfun)
gs4 <- GeneSet(setName = "Firma_Genes_Prolif", geneIds = Firma_Genes_Prolif)
gs5 <- GeneSet(setName = "Firma_Genes_Teff", geneIds = Firma_Genes_Teff)
gs6 <- GeneSet(setName = "Firma_GlucoliticProcess", geneIds = Firma_GlucoliticProcess)
gs7 <- GeneSet(setName = "Firma_CD8_CARHigh_Signature", geneIds = Firma_CD8_CARHigh_Signature)
gs8 <- GeneSet(setName = "Firma_KEGG1", geneIds = Firma_KEGG1)
gs9 <- GeneSet(setName = "Firma_KEGG2", geneIds = Firma_KEGG2)
gs10 <- GeneSet(setName = "Firma_Exhaustion1", geneIds = Firma_Exhaustion1)
gs11 <- GeneSet(setName = "Firma_Exhau_2", geneIds = Firma_Exhau_2)
gs12 <- GeneSet(setName = "Firma_Temeff_TnTcm", geneIds = Firma_Temeff_TnTcm)
gs13 <- GeneSet(setName = "Firma_Gluconeogenesis", geneIds = Firma_Gluconeogenesis)
gs14 <- GeneSet(setName = "Firma_TricarboxylicAcidCycle", geneIds = Firma_TricarboxylicAcidCycle)
gs15 <- GeneSet(setName = "Exhaustion Summarized", geneIds = Firma_Exhaustion_Summarized)
gs16 <- GeneSet(setName = "Cytotoxicity", geneIds = Firma_Cytotoxicity)
gs17 <- GeneSet(setName = "Memory", geneIds = Firma_Memory)

gsc <- GeneSetCollection(gs1, gs2, gs3, gs4, gs5, gs6, gs7, gs8, gs9, gs10, gs11, gs12, gs13, gs14, gs15, gs16, gs17)

res.zenith <- zenith_gsa(res.compare, gsc,
  coef = "compare",
  n_genes_min = 2
)

res.zenith$assay <- factor(res.zenith$assay, names(res.compare))

# plot results, but with no limit based on the highest/lowest t-statistic
setwd("/home/scamara/data_a/scamara/Atlas/Resultados_V5/BCMA_vs_CD19_IP_All")
cairo_pdf("BCMA_vs_CD19_IP_All_Custom_Genesets.pdf", width = 10, height = 10)
plotZenithResults(res.zenith)
dev.off()

Custom_genesets <- res.zenith

###################################################################################################################################################################################################

## Filter Custom_genesets
Custom_genesets$Geneset

# Filter by
patrones <- c("Memory", "Cytotoxicity", "Firma_Genes_Teff", "Firma_Exhau_2")

# Put all in a variable with or
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
Custom_genesets_filtered <- Custom_genesets[grep(patron_combinado, Custom_genesets$Geneset), ]

###################################################################################################################################################################################################

## Final filtered
Final_result <- rbind(Custom_genesets_filtered)

head(Final_result)

setwd("/home/scamara/data_a/scamara/Atlas/Resultados_V5/BCMA_vs_CD19_IP_All")
if(TRUE){
  write.csv(Final_result, "Final_result_Fig_5C_IP.csv", row.names = FALSE)
} else{
  Final_result <- read.csv("Final_result_Fig_5C_IP.csv")
}

p6 <- plotZenithResults(Final_result, Inf, Inf, sortByGeneset = TRUE) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-4.5, 4.5),
    name = "t-statistic"
  )

# Plot results, but with no limit based on the highest/lowest t-statistic
setwd("/home/scamara/data_a/scamara/Atlas/Resultados_V5/BCMA_vs_CD19_IP_All")
cairo_pdf("BCMA_vs_CD19_IP_All_Final_Summarized_Result_1.pdf", width=10, height = 10)
p6
dev.off()

data_plot <- p6$data
data_plot <- add_NTotalGenes(data_plot, gsc)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
cairo_pdf("BCMA_vs_CD19_IP_All_Final_Summarized_Result_2.pdf", width=10, height = 10)
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

################################
######## END OF SCRIPT #########
################################