###############################################################################
###############################################################################

# Program: Supplementary_Table_8.R
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
Final_result <- read.csv("Final_result_Fig_2C.csv") # This object is generated in Dreamlet_CD8_short_anytime_CR.R script
Final_result

library(dplyr)

df <- Final_result

results <- df %>%
  group_by(Geneset) %>%
  summarise(
    delta1 = dplyr::first(delta),
    se1 = dplyr::first(se),
    delta2 = dplyr::last(delta),
    se2 = dplyr::last(se),
    .groups = "drop"
  ) %>%
  mutate(
    Z = (delta1 - delta2) / sqrt(se1^2 + se2^2),
    pval = 1 - pnorm(abs(Z)), # One-tail test
    sig = case_when(
      pval <= 0.001  ~ "***",
      pval <= 0.01   ~ "**",
      pval <= 0.05   ~ "*",
      TRUE           ~ "ns"
    )
  )

results
write.csv(results, "Supplementary_Table_8.csv")

################################
######## END OF SCRIPT #########
################################