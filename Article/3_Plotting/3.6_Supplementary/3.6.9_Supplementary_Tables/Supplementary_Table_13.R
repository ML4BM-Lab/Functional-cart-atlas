###############################################################################
###############################################################################

# Program: Supplementary_Table_13.R
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
Final_result_IP <- read.csv("Final_result_Fig_2J_IP.csv") # This object is generated in Dreamlet_V5_BCMA_vs_CD19_IP_All.R script
Final_result_IP

setwd("/home/scamara/data_a/scamara/Atlas/Resultados_V5/BCMA_vs_CD19_MID")
Final_result_MID <- read.csv("Final_result_Fig_2J_MID.csv") # This object is generated in Dreamlet_V5_BCMA_vs_CD19_MID.R script
Final_result_MID

# Add ID column
Final_result_IP  <- Final_result_IP  %>% mutate(assay = "IP")
Final_result_MID <- Final_result_MID %>% mutate(assay = "MID")

# Merge all in one df
Final_result_all <- bind_rows(Final_result_IP, Final_result_MID)

# Compare between IP and MID
comparison <- Final_result_all %>%
  group_by(Geneset) %>%
  summarise(
    delta_IP = delta[assay == "IP"],
    se_IP    = se[assay == "IP"],
    delta_MID = delta[assay == "MID"],
    se_MID    = se[assay == "MID"],
    .groups = "drop"
  ) %>%
  mutate(
    Z = (delta_IP - delta_MID) / sqrt(se_IP^2 + se_MID^2),
    pval = 1 - pnorm(abs(Z)), # One-tail test
    sig = case_when(
      pval <= 0.001  ~ "***",
      pval <= 0.01   ~ "**",
      pval <= 0.05   ~ "*",
      TRUE           ~ "ns"
    )
  )

comparison
setwd("/home/scamara/data_a/scamara/Atlas/Resultados_V5/Supplementary_Table_13")
write.csv(comparison, "Supplementary_Table_13.csv")

################################
######## END OF SCRIPT #########
################################