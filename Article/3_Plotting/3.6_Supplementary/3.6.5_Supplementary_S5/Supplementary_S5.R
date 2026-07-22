###############################################################################
###############################################################################

# Program: Supplementary_S5.R
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
library(SMFilter)

library(reticulate)
use_python(python_path)
anndata <- reticulate::import("anndata")

.current_dir <- file.path(project_dir, "Codigo", "Rocinante_DEA", "Funciones")
source(.input_path(.current_dir, "add_NTotalGenes.R"))

set.seed(2504)

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### Read files #####
path <- file.path(project_dir, "Input")
file <- .input_path(path, "Python_scVI_adata_big_V4_state4.h5ad")
adata <- anndata$read_h5ad(file)

sce <- AnnData2SCE(adata, "counts", uns = FALSE, obsm = FALSE, obsp = FALSE)
sce
assay(sce, "counts") %>% max()

print((sce %>% dim())[2])

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### Set PATH to save figs #####
.current_dir <- file.path(project_dir, "Figuras", "Suplementarias")
############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### S5A - Barplot of cells per timepoint and response #####
## Filter object
filtered_sce <- sce[, colData(sce)$Antigen == "Blood"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, colData(filtered_sce)$manual_celltype_annotation_high == "CD8 cytotoxic"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, colData(filtered_sce)$STATUS != "HEALTHY"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, !((colData(filtered_sce)$Time_Point_Ranges == "Infusion_Product") & (colData(filtered_sce)$Stimulated == "YES"))]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, !((colData(filtered_sce)$Stimulated == "YES") & (colData(filtered_sce)$Stimulation_Location == "In_vitro"))]
print((filtered_sce %>% dim())[2])

# Extract metadata
meta <- as.data.frame(colData(filtered_sce))

# Reorder levels with Time_Point_Ranges variable
meta$Time_Point_Ranges <- factor(
  meta$Time_Point_Ranges,
  levels = c("Infusion_Product", "<2_weeks", "2_weeks-3_months", ">3_months"),
  ordered = TRUE
)

# Create counts table
tabla_conteos <- table(meta$Max_Response, meta$Time_Point_Ranges)
tabla_conteos

## NOTE: With this data final figure was done using Graphpad Prism 8

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### S5B - Volcanos between consecutive timepoints and S5D - Heatmap of selected genes #####
## Filter object
filtered_sce <- filtered_sce[,!is.na(filtered_sce$Max_Response)]
filtered_sce <- filtered_sce[, colData(filtered_sce)$Max_Response == "CR"]
print((filtered_sce %>% dim())[2])

colData(filtered_sce)$Time_Point_Ranges <- factor(colData(filtered_sce)$Time_Point_Ranges, levels = c("Infusion_Product", "<2_weeks", "2_weeks-3_months", ">3_months"), ordered=TRUE)
colData(filtered_sce)$Age_Range <- factor(colData(filtered_sce)$Age_Range, levels = c("<20", "20-40", "40-60", ">60"), ordered=TRUE)

colData(filtered_sce)$ScFv %>% table() # BCMA: 17863 | CD19: 59433

############################################################################################################################################################################################################

pb <- aggregateToPseudoBulk(filtered_sce,
    assay = "counts",     
    cluster_id = "Time_Point_Ranges",  
    sample_id = "Product_norm",
    BPPARAM = SnowParam(6, progressbar=TRUE))

# Evaluate the specificity of each gene for each cluster
df_cts <- cellTypeSpecificity(pb)

############################################################################################################################################################################################################
## Comparisons ##

## compare first two timepoints
ct.pairs_1 <- c("<2_weeks", "Infusion_Product")

# run comparison
fit_1 <- dreamletCompareClusters(pb, ct.pairs_1, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_1 <- topTable(fit_1, coef = "compare", number = 10)

head(res_1)

############################################################################################################################################################################################################
##### Fig S5D - compare first two timepoints summarized #####
# Rownames to filter
filas_a_filtrar <- c("DHCR7", "HMGCS1", "MSMO1")

# Filter the dataframe
res_1_bis <- res_1[rownames(res_1) %in% filas_a_filtrar, ]

cairo_pdf(.output_path(.current_dir, "S5D.pdf"), width=10, height = 10)
dreamlet::plotHeatmap(df_cts, genes = rownames(res_1_bis))
dev.off()

############################################################################################################################################################################################################
## S5B - compare next two timepoints
ct.pairs_2 <- c("2_weeks-3_months", "<2_weeks")

# run comparison
fit_2 <- dreamletCompareClusters(pb, ct.pairs_2, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_2 <- topTable(fit_2, coef = "compare", number = 10)

head(res_2)

############################################################################################################################################################################################################
## compare last two timepoints
ct.pairs_3 <- c(">3_months", "2_weeks-3_months")

# run comparison
fit_3 <- dreamletCompareClusters(pb, ct.pairs_3, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_3 <- topTable(fit_3, coef = "compare", number = 10)

head(res_3)

############################################################################################################################################################################################################
# Make a list storing each result with a meaningful name
fitList <- list()

id <- paste0("[", ct.pairs_1[1], "]_vs_[", ct.pairs_1[2], "]")
fitList[[id]] <- fit_1

id <- paste0("[", ct.pairs_2[1], "]_vs_[", ct.pairs_2[2], "]")
fitList[[id]] <- fit_2

id <- paste0("[", ct.pairs_3[1], "]_vs_[", ct.pairs_3[2], "]")
fitList[[id]] <- fit_3

res.compare <- as.dreamletResult(fitList) # https://diseaseneurogenomics.github.io/dreamlet/reference/as.dreamletResult.html?q=dreamletCompareClusters#details
res.compare

cairo_pdf(.output_path(.current_dir, "S5B.pdf"), width=10, height = 10)
plotVolcano(res.compare, coef = "compare", ncol = 4)
dev.off()

############################################################################################################################################################################################################
############################################################################################################################################################################################################
##### S5C - Volcanos between all timepoints against IP #####
pb <- aggregateToPseudoBulk(filtered_sce,
    assay = "counts",     
    cluster_id = "Time_Point_Ranges",  
    sample_id = "Product_norm",
    BPPARAM = SnowParam(6, progressbar=TRUE))

# Evaluate the specificity of each gene for each cluster
df_cts <- cellTypeSpecificity(pb)

############################################################################################################################################################################################################
## Comparisons ##
## compare first two timepoints
ct.pairs_1 <- c("<2_weeks", "Infusion_Product")

# run comparison
fit_1 <- dreamletCompareClusters(pb, ct.pairs_1, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_1 <- topTable(fit_1, coef = "compare", number = 10)

head(res_1)

############################################################################################################################################################################################################
## compare next timepoint against IP
ct.pairs_2 <- c("2_weeks-3_months", "Infusion_Product")

# run comparison
fit_2 <- dreamletCompareClusters(pb, ct.pairs_2, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_2 <- topTable(fit_2, coef = "compare", number = 10)

head(res_2)

############################################################################################################################################################################################################
## compare the last timepoint against IP
ct.pairs_3 <- c(">3_months", "Infusion_Product")

# run comparison
fit_3 <- dreamletCompareClusters(pb, ct.pairs_3, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_3 <- topTable(fit_3, coef = "compare", number = 10)

head(res_3)

############################################################################################################################################################################################################
# Make a list storing each result with a meaningful name
fitList <- list()

id <- paste0("[", ct.pairs_1[1], "]_vs_[", ct.pairs_1[2], "]")
fitList[[id]] <- fit_1

id <- paste0("[", ct.pairs_2[1], "]_vs_[", ct.pairs_2[2], "]")
fitList[[id]] <- fit_2

id <- paste0("[", ct.pairs_3[1], "]_vs_[", ct.pairs_3[2], "]")
fitList[[id]] <- fit_3

res.compare <- as.dreamletResult(fitList) # https://diseaseneurogenomics.github.io/dreamlet/reference/as.dreamletResult.html?q=dreamletCompareClusters#details
res.compare

cairo_pdf(.output_path(.current_dir, "S5C.pdf"), width=10, height = 10)
plotVolcano(res.compare, coef = "compare", ncol = 4)
dev.off()

################################
######## END OF SCRIPT #########
################################