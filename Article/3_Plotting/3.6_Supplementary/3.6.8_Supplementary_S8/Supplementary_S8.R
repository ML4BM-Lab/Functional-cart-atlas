###############################################################################
###############################################################################

# Program: Supplementary_S8.R
# Author: Sergio Cámara Peña
# Date: 27/05/2025
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

.generated_input_paths <- c(
    file.path(project_dir, "Resultados_V5", "BCMA_vs_CD19_MID", "Resultados_V5_BCMA_vs_CD19_MID.RDS")
)

.input_path <- function(directory, ...) {
    path <- file.path(directory, ...)
    if (file.exists(path) || dir.exists(path)) {
        return(path)
    }

    if (path %in% .generated_input_paths) {
        input_path <- file.path(project_dir, "Input", basename(path))
        if (file.exists(input_path)) {
            return(input_path)
        }
        stop(
            paste("Required generated input file does not exist. Checked:", paste(c(path, input_path), collapse = ", ")),
            call. = FALSE
        )
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
library(SingleCellExperiment)
library(zellkonverter)
library(kableExtra)
library(scater)
library(scattermore)
library(EnrichmentBrowser)
library(GSEABase)
library(Cairo)
library(SMFilter)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

library(reticulate)
use_python(python_path)
anndata <- reticulate::import("anndata")

############################################################################################################################################################################################################
############################################################################################################################################################################################################

# Set random seed
set.seed(2504)

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### Read files #####
path <- file.path(project_dir, "Resultados", "Joined_datasets", "Raw_Atlas")
file <- .input_path(path, "Atlas_integ_scArches_FINAL_V5.h5ad")
adata <- anndata$read_h5ad(file)

sce <- AnnData2SCE(adata, "counts", uns = FALSE, obsm = TRUE, obsp = FALSE)
sce
assay(sce, "counts") %>% max()

print((sce %>% dim())[2])

.current_dir <- file.path(project_dir, "Resultados", "Joined_datasets", "Dreamlet_to_ClusterProfiler")
res_all <- readRDS(.input_path(.current_dir, "Resultados_V5_BCMA_vs_CD19_MID.RDS")) # Generated in Dreamlet_V5_BCMA_vs_CD19_MID.R

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### Set PATH to save figs #####
.current_dir <- file.path(project_dir, "Resultados_Figuras", "Suplementarias")
############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### Filter object #####
filtered_sce <- sce[, colData(sce)$Antigen == "Blood"]
print((filtered_sce %>% dim())[2])

filtered_sce <- sce[, colData(sce)$STATUS == "DISEASE"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, !((colData(filtered_sce)$Time_Point_Ranges == "Infusion_Product") & (colData(filtered_sce)$Stimulated == "YES"))]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, !((colData(filtered_sce)$Stimulated == "YES") & (colData(filtered_sce)$Stimulation_Location == "In_vitro"))]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, !is.na(filtered_sce$Max_Response)]
filtered_sce <- filtered_sce[, colData(filtered_sce)$Max_Response == "CR"]
print((filtered_sce %>% dim())[2])

colData(filtered_sce)$ScFv %>% table()
colData(filtered_sce)$Source %>% table()
colData(filtered_sce)$Time_Point_Ranges %>% table()

colData(filtered_sce)$Age_Range <- factor(colData(filtered_sce)$Age_Range, levels = c("<20", "20-40", "40-60", ">60"), ordered = TRUE)

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### S8A - BARPLOT WITH DATA #####
# Convert metadata to a data frame
meta_df <- as.data.frame(colData(filtered_sce))

# Ensure relevant columns exist
if (!all(c("ScFv", "Time_Point_Ranges", "Product_norm") %in% colnames(meta_df))) {
  stop("Missing required columns in colData(filtered_sce)")
}

# Remove duplicates based on ScFv, Time_Point_Ranges, and Product_norm
unique_df <- meta_df %>%
  distinct(ScFv, Time_Point_Ranges, Product_norm, .keep_all = TRUE)

# Count the unique occurrences
counts_df <- unique_df %>%
  group_by(ScFv, Time_Point_Ranges) %>%
  summarise(Sample_Count = n(), .groups = "drop")

# Convert categorical variables to factors
counts_df$ScFv <- as.factor(counts_df$ScFv)
counts_df$Time_Point_Ranges <- factor(counts_df$Time_Point_Ranges,
  levels = c(
    "Infusion_Product", "<2_weeks",
    "2_weeks-3_months", ">3_months"
  )
)

pdf(.output_path(.current_dir, "S8A.pdf"))
ggplot(counts_df, aes(x = Time_Point_Ranges, y = Sample_Count, fill = ScFv)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Sample_Count), 
            vjust = -0.5, 
            position = position_dodge(0.9)) +
  labs(
    title = "Unique Sample Counts by Time Point Ranges and ScFv",
    x = "Time Point Ranges",
    y = "Sample Count"
  ) +
  scale_fill_manual(values = c(
    "CD19" = "#66c2a5",
    "BCMA" = "#fc8d62"
  )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## Sample distribution
counts_df_2 <- unique_df %>%
  group_by(Norm_Patient_Name, Product_norm, Time_Point_Ranges, ScFv) %>%
  summarise(.groups = "drop")

dim(counts_df_2)

counts_df_2 <- counts_df_2[counts_df_2$Time_Point_Ranges %in% c("Infusion_Product", "2_weeks-3_months"), ]
dim(counts_df_2)
print(counts_df_2, n = 100)

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### S8B - UMAP original atlas + Jordana's dataset #####
sce2 <- sce

colData(sce2)$dataset_group <- ifelse(
  colData(sce2)$orig_ident == "Jordana_et_al",
  "Jordana",
  "Original_atlas"
)

reducedDimNames(sce2)[reducedDimNames(sce2) == "X_umap"] <- "UMAP"

pdf(.output_path(.current_dir, "S8B.pdf"), width = 8, height = 4)
plotUMAP(
  sce2,
  colour_by = "dataset_group",
  point_size = 0.1,
  other_fields = "dataset_group"
) +
  facet_wrap(~dataset_group, nrow = 1) +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  )
dev.off()

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### S8C, S8D, and S8E in its own python script "Supplementary_S8C,D,E.py" #####

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### S8F - Stacked barplot with cell composition by group #####
filtered_sce2 <- filtered_sce[, 
  colData(filtered_sce)$Time_Point_Ranges %in% c("Infusion_Product", "2_weeks-3_months")
]

colData(filtered_sce2)$manual_celltype_annotation_high %>% unique()
colData(filtered_sce2)$dummy_CD4_CD8 <- ifelse(
  grepl("^CD4", colData(filtered_sce2)$manual_celltype_annotation_high) |
    colData(filtered_sce2)$manual_celltype_annotation_high == "Regulatory T cells",
  "CD4",
  ifelse(
    grepl("^CD8", colData(filtered_sce2)$manual_celltype_annotation_high),
    "CD8",
    NA
  )
)

percentage_table <- as.data.frame(colData(filtered_sce2)) %>%
  filter(!is.na(dummy_CD4_CD8)) %>%
  group_by(Time_Point_Ranges, ScFv, dummy_CD4_CD8) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(freq = n / sum(n) * 100) %>%
  ungroup()

percentage_table

table_wide <- percentage_table %>%
  tidyr::pivot_wider(
    id_cols = c(Time_Point_Ranges, ScFv),
    names_from = dummy_CD4_CD8,
    values_from = freq,
    values_fill = 0
  )

table_wide

## This data was used to generate the final figure using Graphpad Prism 8

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### S8G - GSEA analysis with clusterProfiler #####

## Check data
res_all %>% head()

# Create a vector named with logFC and gene symbols
geneList <- res_all$logFC
names(geneList) <- rownames(res_all)
geneList <- sort(geneList, decreasing = TRUE)
sig_genes <- rownames(res_all[res_all$adj.P.Val < 0.05, ])

gene.df <- bitr(sig_genes, fromType = "SYMBOL",
        toType = c("ENTREZID"),
        OrgDb = org.Hs.eg.db)

# Filter and re-do geneList with ENTREZID
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

p = dotplot(gsea_res, showCategory = 5) + 
  ggtitle("GSEA - GO Biological Process")

ggsave(.output_path(.current_dir, "S8G.pdf"), plot = p, width = 10, height = 6)

################################
######## END OF SCRIPT #########
################################