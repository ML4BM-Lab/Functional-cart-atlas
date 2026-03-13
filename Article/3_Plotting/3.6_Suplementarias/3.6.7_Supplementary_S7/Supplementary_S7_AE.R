###############################################################################
###############################################################################

# Program: Supplementary_S7_AE.R
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
Final_result_S7E <- readRDS("Final_result_S7E.RDS") # This object comes from "Dreamlet_IACs_vs_Rest.R" script

###############################################################################
###############################################################################

## Set PATH to save figs
setwd("/home/scamara/data_a/scamara/Atlas/Figuras/Suplementarias")

###############################################################################
###############################################################################

## S7A - Complete short time point comparison heatmap
head(Final_result_S7A)

pS7A <- plotZenithResults(Final_result_S7A, Inf, Inf, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf("S7A_CD8_short_Anytime_CR_Final_Summarized_Result_Long.pdf", width=10, height = 10)
pS7A
dev.off()

###############################################################################
###############################################################################

## S7E - IACs vs Rest summarized comparison heatmap
head(Final_result_S7E)

pS7E <- plotZenithResults(Final_result_S7E, Inf, Inf, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf("S7E_IACs_vs_Rest_Summarized_Result.pdf", width=10, height = 10)
pS7E
dev.off()

################################
######## END OF SCRIPT #########
################################