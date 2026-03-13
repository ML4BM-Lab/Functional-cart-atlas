###############################################################################
###############################################################################

# Program: 3H.R
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
setwd("/home/scamara/data_a/scamara/Atlas/Resultados/CD8_Short_Anytime_CR")
Final_result <- read.csv("Final_result_Fig_3H.csv") # This object is generated in Dreamlet_CD8_short_anytime_CR.R script
Final_result

# The values obtained here have been used to generate the final figure in GraphPad Prism 8.0.1

################################
######## END OF SCRIPT #########
################################