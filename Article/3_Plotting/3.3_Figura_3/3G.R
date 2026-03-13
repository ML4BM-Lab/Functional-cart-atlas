###############################################################################
###############################################################################

# Program: 3G.R
# Author: Sergio Cámara Peña
# Date: 25/08/2025
# Version: FINAL

###############################################################################
###############################################################################

##### Import libraries #####
library(tidyverse)
library(cowplot)
library(patchwork)
library(SingleCellExperiment)
library(zellkonverter)
library(kableExtra)
library(scattermore)
library(EnrichmentBrowser)
library(GSEABase)
library(Cairo)
library(zenith)
library(SMFilter)

library(reticulate)
use_python("/usr/bin/python3")
anndata <- reticulate::import("anndata")

setwd("/home/scamara/data_a/scamara/Atlas/Codigo/Rocinante_DEA/Funciones")
source("add_NTotalGenes.R")

###############################################################################
###############################################################################

## Set random seed
set.seed(2504)

###############################################################################
###############################################################################

## Read files
setwd("/home/scamara/data_a/scamara/Atlas/Input")
Final_result_S7A <- readRDS("Final_result_S7A.RDS") # This object comes from "Dreamlet_CD8_short_anytime_CR.R" script

###############################################################################
###############################################################################

## Set PATH to save figs
setwd("/home/scamara/data_a/scamara/Atlas/Figuras/Figura_3")

###############################################################################
###############################################################################

## 3G - Summary of S7A short time point comparison heatmap
head(Final_result_S7A)

Final_result_3G <- Final_result_S7A
selected_genesets <- c(
  "ER74_Activation_Of_Gene_Expression_By_SREBF_(SREBP)_R-HSA-2426168",
  "ER712_Immunoregulatory_Interactions_Between_A_Lymphoid_And_A_non-Lymphoid_Cell_R-HSA-198933",
  "ER49_IL_1_Signaling_Pathway_WP195",
  "ER278_Angiopoietin_Like_Protein_8_Regulatory_Pathway_WP3915",
  "ER2651_Peptidyl-Threonine_Phosphorylation_(GO0018107)",
  "ER4809_Regulation_Of_Transforming_Growth_Factor_Beta_Receptor_Signaling_Pathway_(GO0017015)",
  "ER224_Wnt_Signaling_Pathway_WP363",
  "Oxidative Phosphorylation"
)

Final_result_3G_sel <- Final_result_3G %>%
  filter(Geneset %in% selected_genesets) %>%
  mutate(Geneset = factor(Geneset, levels = selected_genesets))  # respeta tu orden

p3G <- plotZenithResults(Final_result_3G_sel, Inf, Inf, sortByGeneset = FALSE)

cairo_pdf("3G_selected_genesets.pdf", width=10, height=10)
p3G
dev.off()

################################
######## END OF SCRIPT #########
################################