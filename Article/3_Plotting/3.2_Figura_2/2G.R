###############################################################################
###############################################################################

# Program: 2G.R
# Author: Sergio Cámara Peña
# Date: 23/09/2025
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

filtered_sce <- filtered_sce[, colData(filtered_sce)$manual_celltype_annotation_high == "CD8 cytotoxic"]
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
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, colData(filtered_sce)$STATUS == "DISEASE"]
print((filtered_sce %>% dim())[2])

## Do these cells express more IL-10?
# Get unique cell types
perform_fisher_test <- function(cell_type) {
  # Subset the data for the given cell type
  cell_indices <- filtered_sce$manual_celltype_annotation_high == cell_type

  # Create contingency table
  contingency_table <- table(
    IL10_expr = assay(filtered_sce, "counts")["IL10", cell_indices] > 0,
    Response = filtered_sce$Max_Response[cell_indices]
  )
  
  # Ensure both "CR" and "NR" exist in the table
  expected_levels <- c("CR", "NR")
  contingency_table <- contingency_table[, expected_levels, drop = FALSE]
  
  # Ensure the table is 2x2 (force missing rows/columns to be zero)
  if (!all(c(TRUE, FALSE) %in% rownames(contingency_table))) {
    full_table <- matrix(0, nrow = 2, ncol = 2,
                         dimnames = list(c(FALSE, TRUE), expected_levels))
    full_table[rownames(contingency_table), ] <- contingency_table
    contingency_table <- full_table
  }
  
  # Perform Fisher's test
  test_result <- fisher.test(contingency_table)
  
  return(list(
    cell_type = cell_type,
    contingency_table = contingency_table,
    p_value = test_result$p.value,
    odds_ratio = test_result$estimate
  ))
}

cell_types <- filtered_sce$manual_celltype_annotation_high %>% unique()

# Apply function to all cell types
fisher_results <- lapply(cell_types, perform_fisher_test)

# Remove NULL results
fisher_results <- fisher_results[!sapply(fisher_results, is.null)]

# Print results
fisher_results

tab <- fisher_results[[1]]$contingency_table
totals <- colSums(tab)

perc_CR <- tab["TRUE", "CR"] / totals["CR"] * 100
perc_NR <- tab["TRUE", "NR"] / totals["NR"] * 100

print(perc_CR)
print(perc_NR)

##### NOTE: Final graph in figure has been done using the data from here in GraphPad Prism 8 #####

################################
######## END OF SCRIPT #########
################################