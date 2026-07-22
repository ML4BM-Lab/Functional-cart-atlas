###############################################################################
###############################################################################

# Program: Dreamlet_IACs_clusters.R
# Author: Sergio Cámara Peña
# Date: 29/05/2024
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
library(ggplot2)

library(reticulate)
use_python(python_path)
anndata <- reticulate::import("anndata")

.current_dir <- file.path(project_dir, "Codigo", "Rocinante_DEA", "Funciones")
source(.input_path(.current_dir, "add_NTotalGenes.R"))

##### Set random seed #####
set.seed(2504)

##### Read files #####
path <- file.path(project_dir, "Input")
file <- .input_path(path, "Python_scVI_adata_big_V4_state4.h5ad")
adata <- anndata$read_h5ad(file)

sce <- AnnData2SCE(adata, "counts", uns = FALSE, obsm = FALSE, obsp = FALSE)
sce
assay(sce, "counts") %>% max()

print((sce %>% dim())[2])

##### Filter object #####
filtered_sce <- sce[, colData(sce)$Antigen == "Blood"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, colData(filtered_sce)$manual_celltype_annotation_high == "Monocyte-like T cells"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, colData(filtered_sce)$STATUS == "DISEASE"]
print((filtered_sce %>% dim())[2])

colData(filtered_sce)$Max_Response <- factor(colData(filtered_sce)$Max_Response, levels = c("NR", "PR", "CR"), ordered = TRUE)
colData(filtered_sce)$ICANS_Grade_Range <- factor(colData(filtered_sce)$ICANS_Grade_Range, levels = c("3-4", "1-2", "0"), ordered = TRUE)
colData(filtered_sce)$Time_Point_Ranges <- factor(colData(filtered_sce)$Time_Point_Ranges, levels = c("Infusion_Product", "<2_weeks", "2_weeks-3_months"), ordered=TRUE)
colData(filtered_sce)$Age_Range <- factor(colData(filtered_sce)$Age_Range, levels = c("<20", "20-40", "40-60", ">60"), ordered=TRUE)

# Add a new column 'IACs_differenciation' and set all values to 'Post_Infusion'
filtered_sce$IACs_differenciation <- 'Post_Infusion'

# Change values in 'IACs_differenciation' column to 'Infusion_Product' where 'Time_Point_Ranges' is 'Infusion_Product'
filtered_sce$IACs_differenciation[filtered_sce$Time_Point_Ranges == 'Infusion_Product'] <- 'Infusion_Product'

############################################################################################################################################################################################################
############################################################################################################################################################################################################

pb <- aggregateToPseudoBulk(filtered_sce,
    assay = "counts",     
    cluster_id = "IACs_differenciation",  
    sample_id = "Product_norm", 
    BPPARAM = SnowParam(6, progressbar=TRUE))

# Evaluate the specificity of each gene for each cluster
df_cts <- cellTypeSpecificity(pb)

############################################################################################################################################################################################################
############################################################################################################################################################################################################
##### Comparisons #####
############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### compare first two timepoints #####
ct.pairs_1 <- c("Post_Infusion", "Infusion_Product")

# run comparison
fit_1 <- dreamletCompareClusters(pb, ct.pairs_1, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_1 <- topTable(fit_1, coef = "compare", number = 10)
res_all <- topTable(fit_1, coef = "compare", number = Inf)

res_1

# Save
output_dir <- file.path(project_dir, "Resultados", "Dreamlet_Vs_PyDESeq2")
output_file <- .output_path(output_dir, "Dreamlet_de_genes_IACs_clusters.csv")
write.csv(res_1, output_file, row.names = TRUE)

# Significant Dreamlet genes
genes_dreamlet <- rownames(res_all[res_all$adj.P.Val < 0.05,])

.current_dir <- file.path(project_dir, "Para_Cluster_Profiler")

# Save gene list as txt, one gene per line
writeLines(genes_dreamlet, .output_path(.current_dir, "genes_dreamlet.txt"))

# Save RDS
saveRDS(res_all, file = .output_path(.current_dir, "dreamlet_DEG_IACs_Post_vs_Infusion.RDS"))

################################
######## END OF SCRIPT #########
################################
