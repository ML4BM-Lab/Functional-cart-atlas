###############################################################################
###############################################################################

# Program: 2E.R
# Author: Sergio Cámara Peña
# Date: 04/06/2025
# Version: FINAL
# Machine: Rocinante

###############################################################################
###############################################################################

## Command-line and environment path configuration

.path_script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.path_script_dir <- if (length(.path_script_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", .path_script_arg[[1]]), mustWork = FALSE))
} else {
    getwd()
}
.find_article_dir <- function(path) {
    current <- normalizePath(path, mustWork = FALSE)
    repeat {
        if (basename(current) == "Article") {
            return(current)
        }
        article_child <- file.path(current, "Article")
        if (dir.exists(article_child)) {
            return(normalizePath(article_child, mustWork = FALSE))
        }
        parent <- dirname(current)
        if (identical(parent, current)) {
            return(normalizePath(path, mustWork = FALSE))
        }
        current <- parent
    }
}
.article_dir <- .find_article_dir(.path_script_dir)
.path_args <- if (interactive()) character() else commandArgs(trailingOnly = TRUE)
.get_path_arg <- function(option, environment, fallback) {
    equals_prefix <- paste0(option, "=")
    equals_match <- grep(paste0("^", equals_prefix), .path_args, value = TRUE)
    if (length(equals_match) > 0) {
        return(sub(paste0("^", equals_prefix), "", equals_match[[1]]))
    }
    option_index <- match(option, .path_args)
    if (!is.na(option_index)) {
        if (option_index == length(.path_args)) {
            stop(paste("Missing value for", option), call. = FALSE)
        }
        return(.path_args[[option_index + 1]])
    }
    environment_value <- Sys.getenv(environment, unset = "")
    if (nzchar(environment_value)) {
        return(environment_value)
    }
    fallback
}
if (!interactive() && any(.path_args %in% c("-h", "--help"))) {
    cat("Path options:\n")
    cat("  --project-dir DIR   Project root (env: CART_ATLAS_PROJECT_DIR; default: Article directory)\n")
    cat("  --python-path FILE  Python executable for reticulate (env: CART_ATLAS_PYTHON; default: python3 on PATH)\n")
    quit(status = 0)
}
project_dir <- normalizePath(
    .get_path_arg("--project-dir", "CART_ATLAS_PROJECT_DIR", .article_dir),
    mustWork = FALSE
)
if (!dir.exists(project_dir)) {
    stop(
        paste(
            "Project directory does not exist:", project_dir,
            "- set --project-dir or CART_ATLAS_PROJECT_DIR"
        ),
        call. = FALSE
    )
}
python_path <- .get_path_arg("--python-path", "CART_ATLAS_PYTHON", Sys.which("python3"))
if (!nzchar(python_path) || !file.exists(python_path)) {
    stop(
        paste(
            "Python executable does not exist:", python_path,
            "- set --python-path or CART_ATLAS_PYTHON"
        ),
        call. = FALSE
    )
}

.input_path <- function(directory, ...) {
    path <- file.path(directory, ...)
    if (!file.exists(path) && !dir.exists(path)) {
        stop(paste("Required input path does not exist:", path), call. = FALSE)
    }
    path
}
.output_path <- function(directory, ...) {
    path <- file.path(directory, ...)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    path
}

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
use_python(python_path)
anndata <- reticulate::import("anndata")

.current_dir <- file.path(.article_dir, "3_Plotting", "Functions")
source(.input_path(.current_dir, "add_NTotalGenes.R"))

## Set random seed
set.seed(2504)

# Go to AUCell_CD8_IL10_CR.R and AUCell_CD8_IL10_NR.R if you want to change something of the code - In there you will find all additional graphics and tests.

.current_dir <- file.path(project_dir, "Resultados", "AUCell_CD8_IL10_NR")
auc_long_NR <- readRDS(.input_path(.current_dir, "auc_long_NR.RDS"))

.current_dir <- file.path(project_dir, "Resultados", "AUCell_CD8_IL10_CR")
auc_long_CR <- readRDS(.input_path(.current_dir, "auc_long_CR.RDS"))

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
.current_dir <- file.path(project_dir, "Figuras", "Figura_2")
cairo_pdf(.output_path(.current_dir, "Figura_2E.pdf"), width = 12, height = 10)

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