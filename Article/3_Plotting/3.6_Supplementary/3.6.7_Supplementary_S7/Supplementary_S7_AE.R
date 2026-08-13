###############################################################################
###############################################################################

# Program: Supplementary_S7_AE.R
# Author: Sergio Cámara Peña
# Date: 25/08/2025
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
use_python(python_path)
anndata <- reticulate::import("anndata")

.current_dir <- file.path(.article_dir, "3_Plotting", "Functions")
source(.input_path(.current_dir, "add_NTotalGenes.R"))

###############################################################################
###############################################################################

## Set random seed
set.seed(2504)

###############################################################################
###############################################################################

## Read files
.current_dir <- file.path(project_dir, "Input")
Final_result_S7A <- readRDS(.input_path(.current_dir, "Final_result_S7A.RDS")) # This object comes from "Dreamlet_CD8_short_anytime_CR.R" script
Final_result_S7E <- readRDS(.input_path(.current_dir, "Final_result_S7E.RDS")) # This object comes from "Dreamlet_IACs_vs_Rest.R" script

###############################################################################
###############################################################################

## Set PATH to save figs
.current_dir <- file.path(project_dir, "Figuras", "Suplementarias")
###############################################################################
###############################################################################

## S7A - Complete short time point comparison heatmap
head(Final_result_S7A)

pS7A <- plotZenithResults(Final_result_S7A, Inf, Inf, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf(.output_path(.current_dir, "S7A_CD8_short_Anytime_CR_Final_Summarized_Result_Long.pdf"), width=10, height = 10)
pS7A
dev.off()

###############################################################################
###############################################################################

## S7E - IACs vs Rest summarized comparison heatmap
head(Final_result_S7E)

pS7E <- plotZenithResults(Final_result_S7E, Inf, Inf, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf(.output_path(.current_dir, "S7E_IACs_vs_Rest_Summarized_Result.pdf"), width=10, height = 10)
pS7E
dev.off()

################################
######## END OF SCRIPT #########
################################