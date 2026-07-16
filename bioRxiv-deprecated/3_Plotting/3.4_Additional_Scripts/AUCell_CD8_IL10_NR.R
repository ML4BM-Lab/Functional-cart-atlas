###############################################################################
###############################################################################

# Program: AUCell_CD8_IL10_NR.R
# Author: Sergio Cámara Peña
# Date: 20/03/2025
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
library(AUCell)
library(reshape2)
library(ggsignif)

library(reticulate)
use_python("/usr/bin/python3")
anndata <- reticulate::import("anndata")

setwd("/home/scamara/data_a/scamara/Atlas/Codigo/Rocinante_DEA/Funciones")
source("add_NTotalGenes.R")
source("extract_gene_sets.R")

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

filtered_sce <- filtered_sce[, ((colData(filtered_sce)$Time_Point_Ranges == "Infusion_Product") & (colData(filtered_sce)$Stimulated == "NO"))]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, !is.na(filtered_sce$Max_Response)]
filtered_sce <- filtered_sce[, ((colData(filtered_sce)$Max_Response == "NR"))]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, ((colData(filtered_sce)$STATUS == "DISEASE"))]
print((filtered_sce %>% dim())[2])

colData(filtered_sce)$Age_Range <- factor(colData(filtered_sce)$Age_Range, levels = c("<20", "20-40", "40-60", ">60"), ordered=TRUE)
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, ((colData(filtered_sce)$manual_celltype_annotation_high %in% c("CD8 cytotoxic", "CD8 effector memory", "CD4 central memory", "Proliferative T cells", "Regulatory T cells", "CD4 effector memory")) )]
print((filtered_sce %>% dim())[2])

############################################################################################################################################################################################################
############################################################################################################################################################################################################
##### Load signatures #####
############################################################################################################################################################################################################
############################################################################################################################################################################################################

genes_IL10 <- c(
  "TNFRSF1B", "TNFRSF1A", "PTGS2", "CXCL2", "ICAM1", "TIMP1", "CCL22", "FCER2", 
  "TYK2", "CSF3", "CCL2", "IL10RA", "IL12B", "CD86", "IL1A", "CCL20", "IL1R2", 
  "IL1R1", "CD80", "CCR2", "IL1B", "LIF", "IL6", "IL10", "IL1RN", "IL18", "CCR5", 
  "JAK1", "CXCL1", "CCR1", "CSF2", "STAT3", "IL12A", "CXCL10", "PTAFR", "CXCL8", 
  "FPR1", "CCL19", "CSF1", "TNF", "IL10RB", "CCL5", "CCL4", "CCL3L1", "CCL3"
)

# Extract count matrix
counts_matrix <- assay(filtered_sce, "counts")

# Filter genes from the list that are in count matrix
genes_presentes <- genes_IL10[genes_IL10 %in% rownames(counts_matrix)]

# Identify genes with total expresion 0
expresion_total <- rowSums(counts_matrix[genes_presentes, , drop = FALSE])
table(expresion_total)

# Genes that have less than 50 total counts in between all cells are removed to avoid noise
genes_filtrados <- genes_presentes[expresion_total >= 50]

###################################################################################################################################################################################################
###################################################################################################################################################################################################
############## AUCell ################
###################################################################################################################################################################################################
###################################################################################################################################################################################################
# Merge all geneset objects
IL10_Signature <- GeneSet(
  setName = "IL10_Signature",
  geneIds = genes_filtrados
)

# Convert GeneSet to a named list
Gene_Sets_Merged <- list(
  IL10_Signature = geneIds(IL10_Signature)  # Extract gene names as vector
)

# Convert singlecellexperiment to sparse matrix:
exprMatrix <- assay(filtered_sce)
exprMatrix <- as(exprMatrix, "dgCMatrix")

# Calculate enrichment scores
First_Time <- TRUE
setwd("/home/scamara/data_a/scamara/Atlas/Resultados/AUCell_CD8_IL10_NR")
if (First_Time) {
  cells_AUC <- AUCell_run(exprMatrix, Gene_Sets_Merged, BPPARAM = BiocParallel::MulticoreParam(16))
  save(cells_AUC, file = "cells_AUC.RData")
} else {
  load(file = "cells_AUC.RData")
}

head(cells_AUC)

# Extract AUC scores as a matrix
auc_matrix <- as.data.frame(t(getAUC(cells_AUC)))  # Cells in rows, pathways in columns

# Rename pathway columns to start with "AUC_"
colnames(auc_matrix) <- paste0("AUC_", colnames(auc_matrix))

# Add cell name to column
auc_matrix$Cell <- rownames(auc_matrix)

# Add a dummy column to track order
auc_matrix$OrderCheck <- seq_len(nrow(auc_matrix))  # Adds 1,2,3,...

# Convert metadata to a data frame and add explicit Cell column
metadata_df <- as.data.frame(colData(filtered_sce))
metadata_df$Cell <- rownames(metadata_df)  # Create an explicit column for cell names

# Merge using the new Cell column
auc_df <- merge(auc_matrix, metadata_df, by = "Cell", sort = FALSE)

# Check if order is preserved
if (!all(auc_df$OrderCheck == seq_len(nrow(auc_df)))) {
  stop("Warning: Merging shuffled the data!")
} else {
  message("✅ Merge successful! Order is preserved.")
}

# Remove dummy column after check
auc_df$OrderCheck <- NULL

# Set rownames to the 'Cell' column
rownames(auc_df) <- auc_df$Cell

## Reshape data for ggplot (long format)
# Get the column indices for "Cell", "Max_Response", "manual_celltype_annotation_high", and AUC columns
cols_to_select <- c(
  which(colnames(auc_df) == "Cell"),
  which(colnames(auc_df) == "manual_celltype_annotation_high"),
  grep("^AUC_", colnames(auc_df))
)

# Subset the dataframe using these column indices
auc_df_clean <- auc_df[, cols_to_select]

# Melt only the AUC columns
auc_long <- melt(auc_df_clean, 
                 id.vars = c("Cell", "manual_celltype_annotation_high"), 
                 measure.vars = grep("^AUC_", colnames(auc_df_clean), value = TRUE))

# Rename columns for clarity
colnames(auc_long) <- c("Cell", "Cell_Type", "Pathway", "AUC")

# Convert Cell_Type to a factor
auc_long$Cell_Type <- factor(auc_long$Cell_Type)

# Set output directory
setwd("/home/scamara/data_a/scamara/Atlas/Resultados/AUCell_CD8_IL10_NR")

# Save plot as PDF
cairo_pdf("AUC_violin_CD8_IL10_comparison.pdf", width = 12, height = 10)
ggplot(auc_long, aes(x = Cell_Type, y = AUC, fill = Cell_Type)) +
  geom_violin() +
  scale_fill_brewer(palette = "Pastel1") +
  stat_summary(
    fun.data = "mean_sdl", fun.args = list(mult = 1),
    geom = "pointrange", color = "black"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Cell Type", y = "AUC", title = "")
dev.off()

## Save auc_long
setwd("/home/scamara/data_a/scamara/Atlas/Resultados/AUCell_CD8_IL10_NR")
saveRDS(auc_long, "auc_long_NR.RDS")

# ANOVA
anova_result <- aov(AUC ~ Cell_Type, data = auc_long)

# ANOVA summary
summary(anova_result)

print((auc_long %>% dim())[1])

###############################################################################
##### T.test mean difference CD8 vs rest #####
###############################################################################
auc_CD8 <- auc_long

auc_CD8 <- auc_CD8[!(auc_CD8$Cell_Type == "Regulatory T cells"), ]
print((auc_CD8 %>% dim())[1])

auc_CD8$Dummy <- "Rest"
auc_CD8$Dummy[(auc_CD8$Cell_Type == "CD8 cytotoxic")] <- "CD8 cytotoxic"

table(auc_CD8$Dummy)
table(auc_CD8$Cell_Type)

# Filter AUC values for each group
auc_CD8_cytotoxic <- auc_CD8$AUC[auc_CD8$Dummy == "CD8 cytotoxic"]
auc_Rest <- auc_CD8$AUC[auc_CD8$Dummy == "Rest"]

var.test(auc_CD8_cytotoxic, auc_Rest)

# 2-sample t-test
t_test_result <- t.test(auc_CD8_cytotoxic, auc_Rest, var.equal = FALSE)

# Show results
t_test_result

###############################################################################
##### T.test mean difference Regulatory vs rest #####
###############################################################################
auc_Regulatory <- auc_long

auc_Regulatory <- auc_Regulatory[!(auc_Regulatory$Cell_Type == "CD8 cytotoxic"), ]
print((auc_Regulatory %>% dim())[1])

auc_Regulatory$Dummy <- "Rest"
auc_Regulatory$Dummy[(auc_Regulatory$Cell_Type == "Regulatory T cells")] <- "Regulatory T cells"

table(auc_Regulatory$Dummy)
table(auc_Regulatory$Cell_Type)

# Filter AUC values for each group
auc_Regulatory_AUC <- auc_Regulatory$AUC[auc_Regulatory$Dummy == "Regulatory T cells"]
auc_Rest_AUC <- auc_Regulatory$AUC[auc_Regulatory$Dummy == "Rest"]

var.test(auc_Regulatory_AUC, auc_Rest_AUC)

# 2-sample t-test
t_test_result_2 <- t.test(auc_Regulatory_AUC, auc_Rest_AUC, var.equal = FALSE)

# Show results
t_test_result_2

#########################
##### END OF SCRIPT #####
#########################