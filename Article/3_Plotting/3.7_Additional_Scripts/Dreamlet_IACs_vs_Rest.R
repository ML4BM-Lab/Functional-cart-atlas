###############################################################################
###############################################################################

# Program: Dreamlet_IACs_vs_Rest.R
# Author: Sergio Cámara Peña
# Date: 25/09/2024
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
use_python(python_path)
anndata <- reticulate::import("anndata")

.current_dir <- file.path(project_dir, "Codigo", "Rocinante_DEA", "Funciones")
source(.input_path(.current_dir, "add_NTotalGenes.R"))

set.seed(2504)

##### Read files #####
path <- file.path(project_dir, "Input")
file <- .input_path(path, "Python_scVI_adata_big_V4_state4.h5ad")
adata <- anndata$read_h5ad(file)

sce <- AnnData2SCE(adata, "counts", uns = FALSE, obsm = TRUE, obsp = TRUE)
sce
assay(sce, "counts") %>% max()

print((sce %>% dim())[2])

##### Filter object #####
filtered_sce <- sce[, colData(sce)$Antigen == "Blood"]
print((filtered_sce %>% dim())[2])

colData(filtered_sce)$Age_Range <- factor(colData(filtered_sce)$Age_Range, levels = c("<20", "20-40", "40-60", ">60"), ordered=TRUE)

colData(filtered_sce)$manual_celltype_annotation_high %>% table()

# Add a new column in colData
colData(filtered_sce)$IACs_dummy <- ifelse(
  colData(filtered_sce)$manual_celltype_annotation_high == "Monocyte-like T cells",
  "IACs",
  "Rest_of_Cells"
)

colData(filtered_sce)$IACs_dummy %>% table()

############################################################################################################################################################################################################
############################################################################################################################################################################################################

pb <- aggregateToPseudoBulk(filtered_sce,
    assay = "counts",     
    cluster_id = "IACs_dummy",  
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
ct.pairs_1 <- c("IACs", "Rest_of_Cells")

# run comparison
fit_1 <- dreamletCompareClusters(pb, ct.pairs_1, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_1 <- topTable(fit_1, coef = "compare", number = 10)

head(res_1)

.current_dir <- file.path(project_dir, "Resultados", "IACs_vs_Rest")
cairo_pdf(.output_path(.current_dir, "IACs_vs_Rest_Heatmap.pdf"), width=10, height = 10)
dreamlet::plotHeatmap(df_cts, genes = rownames(res_1))#, assays=colnames(df_cts)[2:3])
dev.off()

fitList <- list()

id <- paste0("[", ct.pairs_1[1], "]_vs_[", ct.pairs_1[2], "]")
fitList[[id]] <- fit_1

res.compare <- as.dreamletResult(fitList) # https://diseaseneurogenomics.github.io/dreamlet/reference/as.dreamletResult.html?q=dreamletCompareClusters#details
res.compare

############################################################################################################################################################################################################
############################################################################################################################################################################################################
##### Gene Set Enrichment Analysis #####
############################################################################################################################################################################################################
############################################################################################################################################################################################################

# GO_Biological_Process_2023
go.gs <- getGenesets(
  org = "hsa",
  db = "enrichr",
  lib = "GO_Biological_Process_2023",
  gene.id.type = "SYMBOL",
  return.type = "GeneSetCollection"
)

res.zenith <- zenith_gsa(res.compare, go.gs,
  coef = "compare",
  n_genes_min = 20
)

res.zenith$assay <- factor(res.zenith$assay, names(res.compare))

res.zenith <- add_NTotalGenes(res.zenith, go.gs)

# View the updated res.zenith
head(res.zenith)

p1 <- plotZenithResults(res.zenith, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf(.output_path(.current_dir, "IACs_vs_Rest_Biological_Process_1.pdf"), width=10, height = 10)
p1
dev.off()

data_plot <- p1$data
data_plot <- add_NTotalGenes(data_plot, go.gs)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
cairo_pdf(.output_path(.current_dir, "IACs_vs_Rest_Biological_Process_2.pdf"), width=10, height = 10)
ggplot(data_plot, aes(x = assay, y = Geneset, size=logFDR, color = delta)) +
  geom_point() +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Delta") +
  labs(
    x = "Comparison",
    y = "GO term",
    title = ""
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
  )
dev.off()

GO_Biological_Process_2023 <- res.zenith

###################################################################################################################################################################################################
###################################################################################################################################################################################################

# KEGG_2021_Human
go.gs <- getGenesets(
  org = "hsa",
  db = "enrichr",
  lib = "KEGG_2021_Human",
  gene.id.type = "SYMBOL",
  return.type = "GeneSetCollection"
)

res.zenith <- zenith_gsa(res.compare, go.gs,
  coef = "compare",
  n_genes_min = 20
)

res.zenith$assay <- factor(res.zenith$assay, names(res.compare))

res.zenith <- add_NTotalGenes(res.zenith, go.gs)

# View the updated res.zenith
head(res.zenith)

p2 <- plotZenithResults(res.zenith, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
pdf(.output_path(.current_dir, "IACs_vs_Rest_KEGG_2021_Human_1.pdf"))
p2
dev.off()

data_plot <- p2$data
data_plot <- add_NTotalGenes(data_plot, go.gs)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
pdf(.output_path(.current_dir, "IACs_vs_Rest_KEGG_2021_Human_2.pdf"))
ggplot(data_plot, aes(x = assay, y = Geneset, size=logFDR, color = delta)) +
  geom_point() +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Delta") +
  labs(
    x = "Comparison",
    y = "GO term",
    title = ""
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
  )
dev.off()

KEGG_2021_Human <- res.zenith

###################################################################################################################################################################################################
###################################################################################################################################################################################################

# Reactome_2022
go.gs <- getGenesets(
  org = "hsa",
  db = "enrichr",
  lib = "Reactome_2022",
  gene.id.type = "SYMBOL",
  return.type = "GeneSetCollection"
)

res.zenith <- zenith_gsa(res.compare, go.gs,
  coef = "compare",
  n_genes_min = 20
)

res.zenith$assay <- factor(res.zenith$assay, names(res.compare))

res.zenith <- add_NTotalGenes(res.zenith, go.gs)

# View the updated res.zenith
head(res.zenith)

p3 <- plotZenithResults(res.zenith, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf(.output_path(.current_dir, "IACs_vs_Rest_Reactome_2022_1.pdf"), width=10, height = 10)
p3
dev.off()

data_plot <- p3$data
data_plot <- add_NTotalGenes(data_plot, go.gs)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
cairo_pdf(.output_path(.current_dir, "IACs_vs_Rest_Reactome_2022_2.pdf"), width=10, height = 10)
ggplot(data_plot, aes(x = assay, y = Geneset, size=logFDR, color = delta)) +
  geom_point() +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Delta") +
  labs(
    x = "Comparison",
    y = "GO term",
    title = ""
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
  )
dev.off()

Reactome_2022 <- res.zenith

###################################################################################################################################################################################################
###################################################################################################################################################################################################

# WikiPathway_2023_Human
go.gs <- getGenesets(
  org = "hsa",
  db = "enrichr",
  lib = "WikiPathway_2023_Human",
  gene.id.type = "SYMBOL",
  return.type = "GeneSetCollection"
)

res.zenith <- zenith_gsa(res.compare, go.gs,
  coef = "compare",
  n_genes_min = 20
)

res.zenith$assay <- factor(res.zenith$assay, names(res.compare))

res.zenith <- add_NTotalGenes(res.zenith, go.gs)

# View the updated res.zenith
head(res.zenith)

p4 <- plotZenithResults(res.zenith, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf(.output_path(.current_dir, "IACs_vs_Rest_WikiPathway_2023_Human_1.pdf"), width=10, height = 10)
p4
dev.off()

data_plot <- p4$data
data_plot <- add_NTotalGenes(data_plot, go.gs)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
cairo_pdf(.output_path(.current_dir, "IACs_vs_Rest_WikiPathway_2023_Human_2.pdf"), width=10, height = 10)
ggplot(data_plot, aes(x = assay, y = Geneset, size=logFDR, color = delta)) +
  geom_point() +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Delta") +
  labs(
    x = "Comparison",
    y = "GO term",
    title = ""
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
  )
dev.off()

WikiPathway_2023_Human <- res.zenith

###################################################################################################################################################################################################
###################################################################################################################################################################################################
# Custom genesets
.current_dir <- file.path(project_dir, "Input")
Firmas_atlas_df <- read.csv(.input_path(.current_dir, "Firmas_Atlas.csv"), header = TRUE, sep = ";")
Firmas_atlas_df %>% colnames()

Firma_Activation <- Firmas_atlas_df$Activation
Firma_Activation[Firma_Activation == ""] <- NA
Firma_Activation <- Firma_Activation[!is.na(Firma_Activation)]
Firma_Activation

Firma_Tonic <- Firmas_atlas_df$Tonic
Firma_Tonic[Firma_Tonic == ""] <- NA
Firma_Tonic <- Firma_Tonic[!is.na(Firma_Tonic)]
Firma_Tonic

Firma_Genes_Polyfun <- Firmas_atlas_df$Genes_Polyfun
Firma_Genes_Polyfun[Firma_Genes_Polyfun == ""] <- NA
Firma_Genes_Polyfun <- Firma_Genes_Polyfun[!is.na(Firma_Genes_Polyfun)]
Firma_Genes_Polyfun

Firma_Genes_Prolif <- Firmas_atlas_df$Genes_Prolif
Firma_Genes_Prolif[Firma_Genes_Prolif == ""] <- NA
Firma_Genes_Prolif <- Firma_Genes_Prolif[!is.na(Firma_Genes_Prolif)]
Firma_Genes_Prolif

Firma_Genes_Teff <- Firmas_atlas_df$Genes_Teff
Firma_Genes_Teff[Firma_Genes_Teff == ""] <- NA
Firma_Genes_Teff <- Firma_Genes_Teff[!is.na(Firma_Genes_Teff)]
Firma_Genes_Teff

Firma_GlucoliticProcess <- Firmas_atlas_df$GlucoliticProcess
Firma_GlucoliticProcess[Firma_GlucoliticProcess == ""] <- NA
Firma_GlucoliticProcess <- Firma_GlucoliticProcess[!is.na(Firma_GlucoliticProcess)]
Orig_Firma_GlucoliticProcess <- Firma_GlucoliticProcess
Firma_GlucoliticProcess <- unique(Firma_GlucoliticProcess)
Firma_GlucoliticProcess

Firma_CD8_CARHigh_Signature <- Firmas_atlas_df$CD8_CARHigh_Signature
Firma_CD8_CARHigh_Signature[Firma_CD8_CARHigh_Signature == ""] <- NA
Firma_CD8_CARHigh_Signature <- Firma_CD8_CARHigh_Signature[!is.na(Firma_CD8_CARHigh_Signature)]
Firma_CD8_CARHigh_Signature

Firma_KEGG1 <- Firmas_atlas_df$KEGG1
Firma_KEGG1[Firma_KEGG1 == ""] <- NA
Firma_KEGG1 <- Firma_KEGG1[!is.na(Firma_KEGG1)]
Firma_KEGG1

Firma_KEGG2 <- Firmas_atlas_df$KEGG2
Firma_KEGG2[Firma_KEGG2 == ""] <- NA
Firma_KEGG2 <- Firma_KEGG2[!is.na(Firma_KEGG2)]
Firma_KEGG2

Firma_Exhaustion1 <- Firmas_atlas_df$Exhaustion1
Firma_Exhaustion1[Firma_Exhaustion1 == ""] <- NA
Firma_Exhaustion1 <- Firma_Exhaustion1[!is.na(Firma_Exhaustion1)]
Orig_Firma_Exhaustion1 <- Firma_Exhaustion1
Firma_Exhaustion1 <- unique(Firma_Exhaustion1)
Firma_Exhaustion1

Firma_Exhau_2 <- Firmas_atlas_df$Exhau_2
Firma_Exhau_2[Firma_Exhau_2 == ""] <- NA
Firma_Exhau_2 <- Firma_Exhau_2[!is.na(Firma_Exhau_2)]
Firma_Exhau_2

Firma_Temeff_TnTcm <- Firmas_atlas_df$Temeff_TnTcm
Firma_Temeff_TnTcm[Firma_Temeff_TnTcm == ""] <- NA
Firma_Temeff_TnTcm <- Firma_Temeff_TnTcm[!is.na(Firma_Temeff_TnTcm)]
Firma_Temeff_TnTcm

Firma_Gluconeogenesis <- Firmas_atlas_df$Gluconeogenesis
Firma_Gluconeogenesis[Firma_Gluconeogenesis == ""] <- NA
Firma_Gluconeogenesis <- Firma_Gluconeogenesis[!is.na(Firma_Gluconeogenesis)]
Orig_Firma_Gluconeogenesis <- Firma_Gluconeogenesis
Firma_Gluconeogenesis <- unique(Firma_Gluconeogenesis)
Firma_Gluconeogenesis

Firma_TricarboxylicAcidCycle <- Firmas_atlas_df$TricarboxylicAcidCycle
Firma_TricarboxylicAcidCycle[Firma_TricarboxylicAcidCycle == ""] <- NA
Firma_TricarboxylicAcidCycle <- Firma_TricarboxylicAcidCycle[!is.na(Firma_TricarboxylicAcidCycle)]
Orig_Firma_TricarboxylicAcidCycle <- Firma_TricarboxylicAcidCycle
Firma_TricarboxylicAcidCycle <- unique(Firma_TricarboxylicAcidCycle)
Firma_TricarboxylicAcidCycle

gs1 <- GeneSet(setName="Activation", geneIds=Firma_Activation)
gs2 <- GeneSet(setName="Tonic", geneIds=Firma_Tonic)
gs3 <- GeneSet(setName="Genes Polyfun", geneIds=Firma_Genes_Polyfun)
gs4 <- GeneSet(setName="Genes Prolif", geneIds=Firma_Genes_Prolif)
gs5 <- GeneSet(setName="Genes Teff", geneIds=Firma_Genes_Teff)
gs6 <- GeneSet(setName="Glucolitic Process", geneIds=Firma_GlucoliticProcess)
gs7 <- GeneSet(setName="CD8 CARHigh Signature", geneIds=Firma_CD8_CARHigh_Signature)
gs8 <- GeneSet(setName="Oxidative Phosphorylation", geneIds=Firma_KEGG1)
gs9 <- GeneSet(setName="Citrate Cycle_TCA Cycle ", geneIds=Firma_KEGG2)
gs10 <- GeneSet(setName="Exhaustion1", geneIds=Firma_Exhaustion1)
gs11 <- GeneSet(setName="Exhaustion2", geneIds=Firma_Exhau_2)
gs12 <- GeneSet(setName="Temeff TnTcm", geneIds=Firma_Temeff_TnTcm)
gs13 <- GeneSet(setName="Gluconeogenesis", geneIds=Firma_Gluconeogenesis)
gs14 <- GeneSet(setName="TricarboxylicAcidCycle", geneIds=Firma_TricarboxylicAcidCycle)

gsc <- GeneSetCollection(gs1, gs2, gs3, gs4, gs5, gs6, gs7, gs8, gs9, gs10, gs11, gs12, gs13, gs14)

res.zenith = zenith_gsa( res.compare, gsc, 
              coef = 'compare', 
              n_genes_min=20)

res.zenith$assay = factor(res.zenith$assay, names(res.compare))

res.zenith <- add_NTotalGenes(res.zenith, gsc)

# View the updated res.zenith
head(res.zenith)

p5 <- plotZenithResults(res.zenith, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
.current_dir <- file.path(project_dir, "Resultados", "IACs_vs_Rest")
pdf(.output_path(.current_dir, "IACs_vs_Rest_Custom_Genesets_1.pdf"))
p5
dev.off()

data_plot <- p5$data
data_plot <- add_NTotalGenes(data_plot, gsc)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
pdf(.output_path(.current_dir, "IACs_vs_Rest_Custom_Genesets_2.pdf"))
ggplot(data_plot, aes(x = assay, y = Geneset, size=logFDR, color = delta)) +
  geom_point() +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Delta") +
  labs(
    x = "Comparison",
    y = "GO term",
    title = ""
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
  )
dev.off()

Custom_genesets <- res.zenith

###################################################################################################################################################################################################
###################################################################################################################################################################################################
##### Filtered_Results #####

## Filter Reactome_2022
Reactome_2022$Geneset

# Filter by
patrones <- c("ER750_")

# Put all in a variable with or
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
Reactome_2022_filtered <- Reactome_2022[grep(patron_combinado, Reactome_2022$Geneset), ]

###################################################################################################################################################################################################

## Filter WikiPathway_2023_Human
WikiPathway_2023_Human$Geneset

# Filter by
patrones <- c("ER296_", "ER290_")

# Put all in a variable with or
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
WikiPathway_2023_Human_filtered <- WikiPathway_2023_Human[grep(patron_combinado, WikiPathway_2023_Human$Geneset), ]

###################################################################################################################################################################################################

## Final filtered
Final_result <- rbind(Reactome_2022_filtered, WikiPathway_2023_Human_filtered)

saveRDS(Final_result, .output_path(project_dir, "Input", "Final_result_S7E.RDS"))

head(Final_result)

p6 <- plotZenithResults(Final_result, Inf, Inf, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
cairo_pdf(.output_path(.current_dir, "IACs_vs_Rest_Summarized_Result_1.pdf"), width=10, height = 10)
p6
dev.off()

data_plot <- p6$data
data_plot <- add_NTotalGenes(data_plot, gsc)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
cairo_pdf(.output_path(.current_dir, "IACs_vs_Rest_Summarized_Result_2.pdf"), width=10, height = 10)
ggplot(data_plot, aes(x = assay, y = Geneset, size=logFDR, color = delta)) +
  geom_point() +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Delta") +
  labs(
    x = "Comparison",
    y = "GO term",
    title = ""
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
  )
dev.off()

################################
######## END OF SCRIPT #########
################################