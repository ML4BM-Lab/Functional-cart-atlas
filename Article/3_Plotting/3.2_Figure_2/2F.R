###############################################################################
###############################################################################

# Program: 2F.R
# Author: Sergio Cámara Peña
# Date: 23/07/2025
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

library(reticulate)
use_python(python_path)
anndata <- reticulate::import("anndata")

.current_dir <- file.path(project_dir, "Codigo", "Rocinante_DEA", "Funciones")
source(.input_path(.current_dir, "add_NTotalGenes.R"))

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

# filtered_sce <- filtered_sce[, colData(filtered_sce)$manual_celltype_annotation_high == "CD8 cytotoxic"]
# print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, ((colData(filtered_sce)$Time_Point_Ranges == "Infusion_Product") & (colData(filtered_sce)$Stimulated == "NO"))]
print((filtered_sce %>% dim())[2])

colData(filtered_sce)$Age_Range <- factor(colData(filtered_sce)$Age_Range, levels = c("<20", "20-40", "40-60", ">60"), ordered=TRUE)

## Remove NAs to avoid any further problems
table(filtered_sce$Anytime_CR, useNA = "ifany")
filtered_sce <- filtered_sce[, !is.na(filtered_sce$Anytime_CR)]
table(filtered_sce$Anytime_CR, useNA = "ifany")
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, colData(filtered_sce)$STATUS == "DISEASE"]
print((filtered_sce %>% dim())[2])

############################################################################################################################################################################################################
############################################################################################################################################################################################################

pb <- aggregateToPseudoBulk(filtered_sce,
    assay = "counts",     
    cluster_id = "Anytime_CR",  
    sample_id = "Product_norm",
    BPPARAM = SnowParam(6, progressbar=TRUE))

# Evaluate the specificity of each gene for each cluster
df_cts <- cellTypeSpecificity(pb)

############################################################################################################################################################################################################
############################################################################################################################################################################################################
##### Comparisons #####
############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### compare first two #####
ct.pairs_1 <- c("Yes", "No")

# run comparison
fit_1 <- dreamletCompareClusters(pb, ct.pairs_1, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_1 <- topTable(fit_1, coef = "compare", number = 10)
res <- topTable(fit_1, coef = "compare", number = Inf)

head(res_1)

############################################################################################################################################################################################################
############################################################################################################################################################################################################

# Filter to keep only IL10 pathway genes
genes_IL10 <- c(
  "TNFRSF1B", "TNFRSF1A", "PTGS2", "CXCL2", "ICAM1", "TIMP1", "CCL22", "FCER2", 
  "TYK2", "CSF3", "CCL2", "IL10RA", "IL12B", "CD86", "IL1A", "CCL20", "IL1R2", 
  "IL1R1", "CD80", "CCR2", "IL1B", "LIF", "IL6", "IL10", "IL1RN", "IL18", "CCR5", 
  "JAK1", "CXCL1", "CCR1", "CSF2", "STAT3", "IL12A", "CXCL10", "PTAFR", "CXCL8", 
  "FPR1", "CCL19", "CSF1", "TNF", "IL10RB", "CCL5", "CCL4", "CCL3L1", "CCL3"
)

res %>% head()

# Make sure that res has a Gene field with the names
res$Gene <- rownames(res)

# Mark the IL10 genes
res <- res %>%
  mutate(
    IL10_gene = ifelse(Gene %in% genes_IL10, "IL10-related", "Other"),
    color = ifelse(IL10_gene == "IL10-related", "red", "grey")
  )

# Calculate the -log10(pvalue)
res <- res %>%
  mutate(logP = -log10(P.Value))

# Select the most highly expressed IL10-related genes to label
res_to_label <- res %>%
  filter(IL10_gene == "IL10-related") %>%
  arrange(P.Value) %>%
  slice_head(n = 10)  # The 10 most significant genes

# Volcano plot with ggplot2
.current_dir <- file.path(project_dir, "Figuras", "Figura_2")
CairoPDF(.output_path(.current_dir, "Figura_2F.pdf"), width = 10, height = 10)
ggplot(res, aes(x = logFC, y = -log10(P.Value))) +
  # Plot the gray points (other genes) first
  geom_point(data = res %>% filter(IL10_gene == "Other"),
             color = "grey", alpha = 0.6, size = 2) +
  # Then the red points (IL10 genes)
  geom_point(data = res %>% filter(IL10_gene == "IL10-related"),
             color = "darkred", alpha = 0.8, size = 2) +
  # Add labels only to the most significant IL10-related genes
  geom_text(data = res_to_label,
            aes(label = Gene),
            vjust = -1,
            size = 3) +
  theme_minimal() +
  labs(
    x = "log2 Fold Change",
    y = "-log10(P-value)",
    title = "Volcano plot - IL10 genes highlighted"
  )
dev.off()

##### END OF SCRIPT #####
