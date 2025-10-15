###############################################################################
###############################################################################

# Program: 1G.R
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

# Define which categories do I want to keep
tipos_keep <- c("CD4 central memory", "CD8 cytotoxic", "Regulatory T cells")

auc_long_NR <- auc_long_NR %>%
  dplyr::filter(Cell_Type %in% tipos_keep)

auc_long_CR <- auc_long_CR %>%
  dplyr::filter(Cell_Type %in% tipos_keep)

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

# Save plot as PDF
setwd("/home/scamara/data_a/scamara/Atlas/Figuras/Figura_1")
cairo_pdf("Figura_1G.pdf", width = 12, height = 10)

ggplot(auc_combined, aes(x = Cell_Type, y = AUC, fill = Group)) +
  geom_violinhalf(
    data = subset(auc_combined, Group == "NR"),
    alpha = 0.6, color = NA
  ) +
  geom_violinhalf(
    data = subset(auc_combined, Group == "CR"),
    alpha = 0.6, color = NA, flip = TRUE
  ) +
  stat_summary(
    fun.data = "mean_sdl", fun.args = list(mult = 1),
    geom = "pointrange", position = position_dodge(width = 0.6),
    color = "black"
  ) +
  scale_fill_manual(
    values = c("CR" = "#006400", "NR" = "#fc8d62")
  ) +
  coord_cartesian(ylim = c(0, 0.27)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Cell Type", y = "AUC", title = "", fill = "Response")

dev.off()

################################
######## END OF SCRIPT #########
################################