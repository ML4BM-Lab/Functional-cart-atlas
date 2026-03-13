###############################################################################
###############################################################################

# Program: Supplementary_Table_5.R
# Author: Sergio Cámara Peña
# Date: 04/06/2025
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
library(see)

library(reticulate)
use_python("/usr/bin/python3")
anndata <- reticulate::import("anndata")

setwd("/home/scamara/data_a/scamara/Atlas/Codigo/Rocinante_DEA/Funciones")
source("add_NTotalGenes.R")

## Set random seed
set.seed(2504)

# Go to AUCell_CD8_IL10_CR.R and AUCell_CD8_IL10_NR.R if you want to change something of the code - In there you will find all additional graphics and tests.

setwd("/home/scamara/data_a/scamara/Atlas/Resultados/AUCell_CD8_IL10_NR")
auc_long_NR <- readRDS("auc_long_NR.RDS")

setwd("/home/scamara/data_a/scamara/Atlas/Resultados/AUCell_CD8_IL10_CR")
auc_long_CR <- readRDS("auc_long_CR.RDS")

# Add condition labels
auc_long_NR$Group <- "NR"
auc_long_CR$Group <- "CR"

# Combine
auc_combined <- rbind(auc_long_NR, auc_long_CR)

## Wilcoxon additional test
wilcoxon_results <- auc_combined %>%
  group_by(Cell_Type) %>%
  summarise(
    p_value = wilcox.test(AUC[Group == "NR"], AUC[Group == "CR"])$p.value,
    NR_median = median(AUC[Group == "NR"]),
    CR_median = median(AUC[Group == "CR"]),
    N_NR = sum(Group == "NR"),
    N_CR = sum(Group == "CR")
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

wilcoxon_results <- wilcoxon_results %>%
  mutate(delta_median = NR_median - CR_median)

print(wilcoxon_results)

# Get asterisk significance
get_significance <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# Add significance column
wilcoxon_results <- wilcoxon_results %>%
  mutate(Significance = sapply(p_adj, get_significance))

# Print table
wilcoxon_results %>%
  dplyr::select(Cell_Type, N_NR, N_CR, NR_median, CR_median, delta_median, p_value, p_adj, Significance) %>%
  print(n = Inf, width = Inf)

if (FALSE) {
  write.csv(
    wilcoxon_results %>%
      dplyr::select(Cell_Type, N_NR, N_CR, NR_median, CR_median, delta_median, p_value, p_adj, Significance),
    file = "Supplementary_Table_5.csv",
    row.names = FALSE
  )
}

##################### END OF SCRIPT #####################
