###############################################################################
###############################################################################

# Program: 5C.R
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
setwd("/home/scamara/data_a/scamara/Atlas/Resultados_V5/BCMA_vs_CD19_IP_All")
Final_result_IP <- read.csv("Final_result_Fig_5C_IP.csv") # This object is generated in Dreamlet_V5_BCMA_vs_CD19_IP_All.R script
Final_result_IP

setwd("/home/scamara/data_a/scamara/Atlas/Resultados_V5/BCMA_vs_CD19_MID")
Final_result_MID <- read.csv("Final_result_Fig_5C_MID.csv") # This object is generated in Dreamlet_V5_BCMA_vs_CD19_MID.R script
Final_result_MID

# The values obtained here have been used to generate the final figure in GraphPad Prism 8.0.1

################################
######## END OF SCRIPT #########
################################