###############################################################################
###############################################################################

# Program: Supplementary_Table_6.R
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
    if (file.exists(path) || dir.exists(path)) {
        return(path)
    }

    atlas_filenames <- c(
        "Python_scVI_adata_big_V4_state4.h5ad",
        "Python_scVI_adata_big_V4_state4_Normalized.h5ad",
        "Atlas_integ_scArches_FINAL_V5.h5ad"
    )
    atlas_directories <- c(
        file.path(project_dir, "Input"),
        file.path(project_dir, "Resultados", "Joined_datasets", "Raw_Atlas")
    )
    if (basename(path) %in% atlas_filenames && directory %in% atlas_directories) {
        alternate_paths <- file.path(
            atlas_directories[atlas_directories != directory],
            basename(path)
        )
        for (alternate_path in alternate_paths) {
            if (file.exists(alternate_path)) {
                return(alternate_path)
            }
        }
        checked_paths <- c(path, alternate_paths)
        stop(
            paste(
                "Required atlas input file does not exist. Checked:",
                paste(checked_paths, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    stop(paste("Required input path does not exist:", path), call. = FALSE)
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

.current_dir <- file.path(project_dir, "Codigo", "Rocinante_DEA", "Funciones")
source(.input_path(.current_dir, "add_NTotalGenes.R"))

## Set random seed
set.seed(2504)

##### Supplementary Table 6 #####
## Read files
path <- file.path(project_dir, "Input")
file <- .input_path(path, "Python_scVI_adata_big_V4_state4.h5ad")
adata <- anndata$read_h5ad(file)

sce <- AnnData2SCE(adata, "counts", uns = FALSE, obsm = FALSE, obsp = FALSE)
sce
assay(sce, "counts") %>% max()

print((sce %>% dim())[2])

## Filter object
filtered_sce <- sce[, colData(sce)$Antigen == "Blood"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, colData(filtered_sce)$Time_Point_Ranges == "Infusion_Product"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, !((colData(filtered_sce)$Time_Point_Ranges == "Infusion_Product") & (colData(filtered_sce)$Stimulated == "YES"))]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, colData(filtered_sce)$STATUS == "DISEASE"]
print((filtered_sce %>% dim())[2])

colData(filtered_sce)$Age_Range <- factor(colData(filtered_sce)$Age_Range, levels = c("<20", "20-40", "40-60", ">60"), ordered=TRUE)

## Remove NAs to avoid any further problems
table(filtered_sce$Anytime_CR, useNA = "ifany")
filtered_sce <- filtered_sce[, !is.na(filtered_sce$Anytime_CR)]
table(filtered_sce$Anytime_CR, useNA = "ifany")

## Do these cells express more IL-10?
# Get unique cell types
perform_fisher_test <- function(cell_type) {
  # Subset the data for the given cell type
  cell_indices <- filtered_sce$manual_celltype_annotation_high == cell_type

  # Create contingency table
  contingency_table <- table(
    IL10_expr = assay(filtered_sce, "counts")["IL10", cell_indices] > 0,
    Response = filtered_sce$Max_Response[cell_indices]
  )
  
  # Ensure both "CR" and "NR" exist in the table
  expected_levels <- c("CR", "NR")
  contingency_table <- contingency_table[, expected_levels, drop = FALSE]
  
  # Ensure the table is 2x2 (force missing rows/columns to be zero)
  if (!all(c(TRUE, FALSE) %in% rownames(contingency_table))) {
    full_table <- matrix(0, nrow = 2, ncol = 2,
                         dimnames = list(c(FALSE, TRUE), expected_levels))
    full_table[rownames(contingency_table), ] <- contingency_table
    contingency_table <- full_table
  }
  
  # Perform Fisher's test
  test_result <- fisher.test(contingency_table)
  
  return(list(
    cell_type = cell_type,
    contingency_table = contingency_table,
    p_value = test_result$p.value,
    odds_ratio = test_result$estimate
  ))
}

cell_types <- filtered_sce$manual_celltype_annotation_high %>% unique()

# Apply function to all cell types
fisher_results <- lapply(cell_types, perform_fisher_test)

# Remove NULL results
fisher_results <- fisher_results[!sapply(fisher_results, is.null)]

# Print results
fisher_results

################################
######## END OF SCRIPT #########
################################