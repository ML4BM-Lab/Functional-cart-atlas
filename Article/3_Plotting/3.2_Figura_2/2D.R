###############################################################################
###############################################################################

# Program: 2D.R
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
setwd("/diska/data/scamara/Atlas/Resultados/CD8_IP_comparison")
Final_result_2D <- readRDS("Final_filtered.RDS") # This object comes from "Dreamlet_IP_comparison.R" script

###############################################################################
###############################################################################

## Set PATH to save figs
setwd("/home/scamara/data_a/scamara/Atlas/Figuras/Figura_2")

###############################################################################
###############################################################################

## 2D - Summary 
head(Final_result_2D)

Final_result_2D_noPR <- Final_result_2D %>%
  dplyr::mutate(assay = as.character(assay)) %>%
  dplyr::filter(assay != "[PR]_vs_[rest]")

head(Final_result_2D_noPR)

table(Final_result_2D$assay)
table(Final_result_2D_noPR$assay)

p6 <- plotZenithResults(Final_result_2D_noPR, Inf, Inf, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf("2D.pdf", width=10, height = 10)
p6
dev.off()

################################
######## END OF SCRIPT #########
################################