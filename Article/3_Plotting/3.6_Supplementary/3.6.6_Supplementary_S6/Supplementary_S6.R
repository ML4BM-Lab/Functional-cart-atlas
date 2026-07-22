###############################################################################
###############################################################################

# Program: Supplementary_S6.R
# Author: Sergio Cámara Peña
# Date: 29/05/2024
# Version: FINAL
# Additional info: To see where the file being read comes from, check the Dreamlet_CD8_clusters.R script
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
.current_dir <- file.path(project_dir, "Input")
Final_result_orig <- readRDS(.input_path(.current_dir, "Final_result_orig.RDS")) # Comes from the Dreamlet_CD8_clusters.R script
res.compare <- readRDS(.input_path(.current_dir, "res_compare.RDS")) # Comes from the Dreamlet_CD8_clusters.R script

##### Load custom genesets #####
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

Firma_Exhaustion_Summarized <- c("TIGIT", "CTLA4", "HAVCR2", "LAG3", "SPATA2") # HAVCR2 = TIM3; SPATA2 = PD1

Firma_Cytotoxicity <- c("PRF1", "IFNG", "NKG7", "GNLY", "GZMA", "GZMK", "GZMB", "GZMH", "CD7", "CCL5", "CD8A", "CCL4")

Firma_Memory <- c("TCF7", "CCR7", "SELL", "IL7R")

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
gs15 <- GeneSet(setName="Exhaustion Summarized", geneIds=Firma_Exhaustion_Summarized)
gs16 <- GeneSet(setName="Cytotoxicity", geneIds=Firma_Cytotoxicity)
gs17 <- GeneSet(setName="Memory", geneIds=Firma_Memory)

gsc <- GeneSetCollection(gs1, gs2, gs3, gs4, gs5, gs6, gs7, gs8, gs9, gs10, gs11, gs12, gs13, gs14, gs15, gs16, gs17)

###### Trend chart ######
lista_genesets <- Final_result_orig$Geneset %>% unique()

## GO_Biological_Process_2023
go.gs <- getGenesets(
  org = "hsa",
  db = "enrichr",
  lib = "GO_Biological_Process_2023",
  gene.id.type = "SYMBOL",
  return.type = "GeneSetCollection"
)

res.zenith_Niveles_GOBP <- zenith_gsa(res.compare, go.gs,
  coef = "cellClusterbaseline",
  n_genes_min = 20
)

res.zenith_filtrado_GOBP <- res.zenith_Niveles_GOBP %>%
  filter(Geneset %in% lista_genesets)

## KEGG_2021_Human
go.gs <- getGenesets(
  org = "hsa",
  db = "enrichr",
  lib = "KEGG_2021_Human",
  gene.id.type = "SYMBOL",
  return.type = "GeneSetCollection"
)

res.zenith_Niveles_KEGG <- zenith_gsa(res.compare, go.gs,
  coef = "cellClusterbaseline",
  n_genes_min = 20
)

res.zenith_filtrado_KEGG <- res.zenith_Niveles_KEGG %>%
  filter(Geneset %in% lista_genesets)

## Reactome_2022
go.gs <- getGenesets(
  org = "hsa",
  db = "enrichr",
  lib = "Reactome_2022",
  gene.id.type = "SYMBOL",
  return.type = "GeneSetCollection"
)

res.zenith_Niveles_Reactome <- zenith_gsa(res.compare, go.gs,
  coef = "cellClusterbaseline",
  n_genes_min = 20
)

res.zenith_filtrado_Reactome <- res.zenith_Niveles_Reactome %>%
  filter(Geneset %in% lista_genesets)

## WikiPathway_2023_Human
go.gs <- getGenesets(
  org = "hsa",
  db = "enrichr",
  lib = "WikiPathway_2023_Human",
  gene.id.type = "SYMBOL",
  return.type = "GeneSetCollection"
)

res.zenith_Niveles_WikiPathway <- zenith_gsa(res.compare, go.gs,
  coef = "cellClusterbaseline",
  n_genes_min = 20
)

res.zenith_filtrado_WikiPathway <- res.zenith_Niveles_WikiPathway %>%
  filter(Geneset %in% lista_genesets)

## Custom genesets
res.zenith_Niveles_Custom <- zenith_gsa(res.compare, gsc,
  coef = "cellClusterbaseline",
  n_genes_min = 2
)

res.zenith_filtrado_Custom <- res.zenith_Niveles_Custom %>%
  filter(Geneset %in% lista_genesets)

## Merge
Baseline <- rbind(res.zenith_filtrado_GOBP, res.zenith_filtrado_KEGG, res.zenith_filtrado_Reactome, res.zenith_filtrado_WikiPathway, res.zenith_filtrado_Custom)

Baseline_Summarized <- Baseline %>% dplyr::select(c("assay", ,"coef", "Geneset", "delta"))

## Final results and comparisons
Baseline_Summarized_FINAL <- Baseline_Summarized[Baseline_Summarized$assay=="[<2_weeks]_vs_[Infusion_Product]",]

Baseline_Summarized_FINAL

## Deltas
Final_result_orig[Final_result_orig$Geneset=="ER856_Cytoplasmic_Translation_(GO0002181)",]
Final_result_orig[Final_result_orig$Geneset=="ER1707_Mitochondrial_Gene_Expression_(GO0140053)",]
Final_result_orig[Final_result_orig$Geneset=="ER4954_Response_To_Interferon-Beta_(GO0035456)",]
Final_result_orig[Final_result_orig$Geneset=="ER3628_Protein_Dephosphorylation_(GO0006470)",]
Final_result_orig[Final_result_orig$Geneset=="ER2701_Phosphatidylinositol-Mediated_Signaling_(GO0048015)",]
Final_result_orig[Final_result_orig$Geneset=="ER229_Proteasome",]
Final_result_orig[Final_result_orig$Geneset=="ER111_Glycolysis__Gluconeogenesis",]
Final_result_orig[Final_result_orig$Geneset=="ER68_Cysteine_and_methionine_metabolism",]
Final_result_orig[Final_result_orig$Geneset=="ER1412_SRP-dependent_Cotranslational_Protein_Targeting_To_Membrane_R-HSA-1799339",]
Final_result_orig[Final_result_orig$Geneset=="ER1188_RAC1_GTPase_Cycle_R-HSA-9013149",]
Final_result_orig[Final_result_orig$Geneset=="ER712_Immunoregulatory_Interactions_Between_A_Lymphoid_And_A_non-Lymphoid_Cell_R-HSA-198933",]
Final_result_orig[Final_result_orig$Geneset=="ER493_Cytoplasmic_Ribosomal_Proteins_WP477",]
Final_result_orig[Final_result_orig$Geneset=="ER450_Cancer_Immunotherapy_By_PD_1_Blockade_WP4585",]
Final_result_orig[Final_result_orig$Geneset=="ER65_TNF_Related_Weak_Inducer_Of_Apoptosis_TWEAK_Signaling_Pathway_WP2036",]
Final_result_orig[Final_result_orig$Geneset=="Genes Prolif",]

##### NOTE: Final graphs in figure have been done using the data from here in GraphPad Prism 8 #####

################################
######## END OF SCRIPT #########
################################