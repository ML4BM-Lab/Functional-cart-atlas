###############################################################################
###############################################################################

# Program: Supplementary_Table_14.R
# Author: Sergio Cámara Peña
# Date: 04/12/2024
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

library(reticulate)
use_python(python_path)
anndata <- reticulate::import("anndata")

.current_dir <- file.path(project_dir, "Codigo", "Rocinante_DEA", "Funciones")
source(.input_path(.current_dir, "add_NTotalGenes.R"))

set.seed(2504)

##### Read files #####
.current_dir <- file.path(project_dir, "Resultados_V5", "BCMA_vs_CD19_IP_All")
Final_result_IP <- read.csv(.input_path(.current_dir, "Final_result_Fig_2J_IP.csv")) # This object is generated in Dreamlet_V5_BCMA_vs_CD19_IP_All.R script
Final_result_IP

.current_dir <- file.path(project_dir, "Resultados_V5", "BCMA_vs_CD19_MID")
Final_result_MID <- read.csv(.input_path(.current_dir, "Final_result_Fig_2J_MID.csv")) # This object is generated in Dreamlet_V5_BCMA_vs_CD19_MID.R script
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
.current_dir <- file.path(project_dir, "Resultados_V5", "Supplementary_Table_14")
write.csv(comparison, .output_path(.current_dir, "Supplementary_Table_14.csv"))

################################
######## END OF SCRIPT #########
################################