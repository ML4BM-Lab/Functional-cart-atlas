###############################################################################
###############################################################################

# Program: Supplementary_S7F.R
# Author: Sergio Cámara Peña
# Date: 26/08/2025
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

##### Load needed libraries #####
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

###############################################################################
###############################################################################

##### Set seed #####
set.seed(2504)

###############################################################################
###############################################################################

##### Load data #####
.current_dir <- file.path(project_dir, "Resultados", "Joined_datasets", "Dreamlet_to_ClusterProfiler")
res_all <- readRDS(.input_path(.current_dir, "dreamlet_DEG_IACs_Post_vs_Infusion.RDS")) # This object comes from "Dreamlet_IACs_clusters.R" script
res_all %>% head()

.current_dir <- file.path(project_dir, "Resultados", "Joined_datasets", "Raw_Atlas")
###############################################################################
###############################################################################

##### Set PATH to save figs
.current_dir <- file.path(project_dir, "Resultados_Figuras", "Suplementarias")
###############################################################################
###############################################################################

### S7F - GSEA GO dotplot of IACs IP vs Post
## Create a vector named with logFC and gene symbols
geneList <- res_all$logFC
names(geneList) <- rownames(res_all)
geneList <- sort(geneList, decreasing = TRUE)
sig_genes <- rownames(res_all[res_all$adj.P.Val < 0.05, ])

gene.df <- bitr(sig_genes, fromType = "SYMBOL",
        toType = c("ENTREZID"),
        OrgDb = org.Hs.eg.db)

## Filter and re-do geneList with ENTREZID
geneList_entrez <- geneList[gene.df$SYMBOL]
names(geneList_entrez) <- gene.df$ENTREZID
geneList_entrez <- sort(geneList_entrez, decreasing = TRUE)

gsea_res <- gseGO(geneList     = geneList_entrez,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "BP",
                  keyType      = "ENTREZID",
                  minGSSize    = 10,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)

p = dotplot(gsea_res, showCategory = 20) + 
  ggtitle("GSEA - GO Biological Process")

ggsave(.output_path(.current_dir, "S7F_GSEA_GO_dotplot_IACs_Pre_vs_Post_20.pdf"), plot = p, width = 10, height = 10)

################################
######## END OF SCRIPT #########
################################