###############################################################################
###############################################################################

# Program: Jordana_et_al_Workflow.R
# Author: Sergio Cámara Peña
# Date: 2023
# Version: FINAL
# Machine: Margaret

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

##### Load libraries #####
library(Seurat)
library(tidyverse)
library(cowplot)
library(scales)
library(patchwork)


##### Set seed #####
set.seed(2504)


##### Load data #####
.current_dir <- file.path(project_dir, "New_Datasets", "Jordana_et_al", "Datos_lorea")
Jordana_Seurat <- readRDS(.input_path(.current_dir, "seurat_new.rds"))


##### Metadata columns to leave #####
# [1] "orig.ident"       "nCount_RNA"       "nFeature_RNA"     "Product"
# [5] "log10GenesPerUMI" "mitoRatio"        "S.Score"          "G2M.Score"
# [9] "Phase"            "Product_norm"


##### Clean metadata #####
# Remove sample with 1 CAR
Jordana_Seurat
Jordana_Seurat <- subset(Jordana_Seurat, subset = orig.ident != "HCB08_PB_M3")
Jordana_Seurat

# Remove columns from the metadata of the Seurat object
Jordana_Seurat@meta.data <- Jordana_Seurat@meta.data %>% select(-nUMI, -nGene, -RPSRatio, -RNA_snn_res.0.2, -RNA_snn_res.0.4, -RNA_snn_res.0.6, -RNA_snn_res.0.8, -seurat_clusters, -subcluster, -subcluster2, -cell_id, -Location, -type, -Time, -Sample, -Patient)

# Rename columns from the metadata of the Seurat object
Jordana_Seurat@meta.data <- Jordana_Seurat@meta.data %>%
    rename(
        Product = orig.ident,
        Lorea_Annotation = Idents
    )

Jordana_Seurat@meta.data$orig.ident <- "Jordana_et_al"
Jordana_Seurat@meta.data %>% head()

## Change forbiden names
# Step 1: Replace "HCB08" with "P2"
Jordana_Seurat@meta.data$Product <- str_replace(Jordana_Seurat@meta.data$Product, "HCB08", "P2")

# Step 2: Replace "HVA02" with "P1"
Jordana_Seurat@meta.data$Product <- str_replace(Jordana_Seurat@meta.data$Product, "HVA02", "P1")

# Step 3: Replace "IBS08" with "P3"
Jordana_Seurat@meta.data$Product <- str_replace(Jordana_Seurat@meta.data$Product, "IBS08", "P3")

table(paste0(Jordana_Seurat@meta.data$Product, "_", Jordana_Seurat@meta.data$ID))

# Remove columns from the metadata of the Seurat object
Jordana_Seurat@meta.data <- Jordana_Seurat@meta.data %>% select(-ID)

Jordana_Seurat@meta.data$Product_norm <- paste0("Jor_", Jordana_Seurat@meta.data$Product)

Jordana_Seurat@meta.data %>% head()

# Re-order columns from the metadata of the Seurat object
Jordana_Seurat@meta.data <- Jordana_Seurat@meta.data %>%
    select(
        orig.ident, nCount_RNA, nFeature_RNA, Product, log10GenesPerUMI, mitoRatio, S.Score, G2M.Score, Phase, Product_norm,
        everything()
    ) # 'everything()' keeps remaining columns

## Change cell names to avoid forbbiden names
# Get current cell names
old_names <- Cells(Jordana_Seurat)

# Step 1, 2, and 3: Create new cell names by replacing parts of the old names
new_names <- old_names %>%
  str_replace("HCB08", "Jor_P2") %>%
  str_replace("HVA02", "Jor_P1") %>%
  str_replace("IBS08", "Jor_P3")

# Create a named vector with old names as names and new names as values
name_mapping <- setNames(new_names, old_names)

Jordana_Seurat@meta.data %>% head()

# Rename cells
Jordana_Seurat <- RenameCells(Jordana_Seurat, new.names = name_mapping)

Jordana_Seurat@meta.data %>% head()
Jordana_Seurat@meta.data %>% tail()

Jordana_Seurat@meta.data["Product_norm"] %>% unique()

# Check everything is fine
Check_1 <- sub("_([^_]+)$", "", rownames(Jordana_Seurat@meta.data))
Check_2 <- Jordana_Seurat@meta.data$Product_norm

all(Check_1 == Check_2)

##### Save RDS object state 1 #####
.current_dir <- file.path(project_dir, "Resultados", "Jordana_et_al", "RDS")
saveRDS(Jordana_Seurat, .output_path(.current_dir, "PostQC_CellRanger_Jordana_RDS.rds"))


##### Normalization #####
.current_dir <- file.path(project_dir, "Resultados", "Jordana_et_al", "RDS")
Jordana_Seurat <- readRDS(.input_path(.current_dir, "PostQC_CellRanger_Jordana_RDS.rds"))

Jordana_Seurat <- NormalizeData(Jordana_Seurat)
Jordana_Seurat <- FindVariableFeatures(Jordana_Seurat,
    selection.method = "vst",
    nfeatures = 2000)
Jordana_Seurat <- ScaleData(Jordana_Seurat)

saveRDS(Jordana_Seurat, file = .output_path(.current_dir, "Normalized_CellRanger_Jordana_RDS.rds"))


##### Merge ---- Worst case scenario #####
.current_dir <- file.path(project_dir, "Resultados", "Jordana_et_al", "RDS")
Jordana_Seurat <- readRDS(.input_path(.current_dir, "Normalized_CellRanger_Jordana_RDS.rds"))

# Run PCA
Jordana_Seurat <- RunPCA(object = Jordana_Seurat, reduction.name = "pca_wo_integ")

# Run UMAP
Jordana_Seurat <- RunUMAP(Jordana_Seurat,
    dims = 1:30,
    reduction = "pca_wo_integ",
    reduction.name = "umap_wo_integ"
)

# Plot UMAP
pdf(.output_path(project_dir, "Resultados", "Jordana_et_al", "Plots", "Worst_case_scenario", "WO_integ_Seurat.pdf"))
DimPlot(Jordana_Seurat,
    group.by = "Product",
    reduction = "umap_wo_integ"
)
dev.off()

# Quality control metrics
p <- FeaturePlot(Jordana_Seurat,
    reduction = "umap_wo_integ",
    features = c("nCount_RNA", "nFeature_RNA", "log10GenesPerUMI", "mitoRatio"),
    order = TRUE
)

pdf(.output_path(project_dir, "Resultados", "Jordana_et_al", "Plots", "Worst_case_scenario", "QC_metrics_WO_integ_Seurat.pdf"))
p + plot_annotation(title = paste0(Jordana_Seurat@meta.data$orig.ident %>% unique()), theme = theme(plot.title = element_text(size = 16)))
dev.off()

rm(p)