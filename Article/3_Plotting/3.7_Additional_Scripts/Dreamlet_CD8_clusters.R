###############################################################################
###############################################################################

# Program: Dreamlet_CD8_clusters.R
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

filtered_sce <- filtered_sce[, colData(filtered_sce)$manual_celltype_annotation_high == "CD8 cytotoxic"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[,!is.na(filtered_sce$Max_Response)]
filtered_sce <- filtered_sce[, colData(filtered_sce)$Max_Response == "CR"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, colData(filtered_sce)$STATUS == "DISEASE"]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, !((colData(filtered_sce)$Time_Point_Ranges == "Infusion_Product") & (colData(filtered_sce)$Stimulated == "YES"))]
print((filtered_sce %>% dim())[2])

filtered_sce <- filtered_sce[, !((colData(filtered_sce)$Stimulated == "YES") & (colData(filtered_sce)$Stimulation_Location == "In_vitro"))]
print((filtered_sce %>% dim())[2])

colData(filtered_sce)$Time_Point_Ranges <- factor(colData(filtered_sce)$Time_Point_Ranges, levels = c("Infusion_Product", "<2_weeks", "2_weeks-3_months", ">3_months"), ordered=TRUE)
colData(filtered_sce)$Age_Range <- factor(colData(filtered_sce)$Age_Range, levels = c("<20", "20-40", "40-60", ">60"), ordered=TRUE)

colData(filtered_sce)$ScFv %>% table() # BCMA: 17863 | CD19: 59433

Metadata_ScFv <- colData(filtered_sce)

# Subset the dataset to retain only the specified columns
Metadata_ScFv_subset <- Metadata_ScFv[, c("Product_norm", "ScFv")]

# Remove duplicate rows
Metadata_ScFv_unique <- unique(Metadata_ScFv_subset)

# View the result
head(Metadata_ScFv_unique)

# Count occurrences of each ScFv category
counts_scfv <- as.data.frame(table(Metadata_ScFv_unique$ScFv))
colnames(counts_scfv) <- c("ScFv", "Count")

# Remove rows where the count is 0
counts_scfv <- counts_scfv %>% filter(Count > 0)

# Create a bar plot with a color palette
ggplot(data = counts_scfv, aes(x = ScFv, y = Count, fill = ScFv)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, size = 4) +
  scale_fill_brewer(palette = "Set1") +  # Use a color palette
  labs(title = "Frequency of ScFv Categories",
       x = "ScFv",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, max(counts_scfv$Count) + 10)
dev.off()

############################################################################################################################################################################################################
############################################################################################################################################################################################################

pb <- aggregateToPseudoBulk(filtered_sce,
    assay = "counts",     
    cluster_id = "Time_Point_Ranges",  
    sample_id = "Product_norm",
    BPPARAM = SnowParam(6, progressbar=TRUE))

# Evaluate the specificity of each gene for each cluster
df_cts <- cellTypeSpecificity(pb)

##### Calculate frobenius distances to know distance #####
Infusion_Product_matrix <- as.matrix(assays(pb)[["Infusion_Product"]])
Two_weeks_matrix <- as.matrix(assays(pb)[["<2_weeks"]])
Two_weeks_3_months_matrix <- as.matrix(assays(pb)[["2_weeks-3_months"]])
Three_months_matrix <- as.matrix(assays(pb)[[">3_months"]])

FDist2(Infusion_Product_matrix,Two_weeks_matrix)
FDist2(Two_weeks_matrix,Two_weeks_3_months_matrix)
FDist2(Two_weeks_3_months_matrix,Three_months_matrix)

FDist2(Infusion_Product_matrix,Two_weeks_3_months_matrix)
FDist2(Infusion_Product_matrix,Three_months_matrix)

############################################################################################################################################################################################################
############################################################################################################################################################################################################
##### Comparisons #####
############################################################################################################################################################################################################
############################################################################################################################################################################################################

.current_dir <- file.path(project_dir, "Resultados", "CD8_clusters")

##### compare first two timepoints #####
ct.pairs_1 <- c("<2_weeks", "Infusion_Product")

# run comparison
fit_1 <- dreamletCompareClusters(pb, ct.pairs_1, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_1 <- topTable(fit_1, coef = "compare", number = 10)
res_all_1 <- topTable(fit_1, coef = "compare", number = Inf)

# Significant Dreamlet genes
genes_dreamlet_1 <- rownames(res_all_1[res_all_1$adj.P.Val < 0.05, ])

# Save gene list as txt, one gene per line
writeLines(genes_dreamlet_1, .output_path(.current_dir, "genes_dreamlet_1.txt"))

res_1

# Save
output_dir <- file.path(project_dir, "Resultados", "Dreamlet_Vs_PyDESeq2")
output_file <- .output_path(output_dir, "Dreamlet_de_genes_CD8_Comparison_1.csv")
write.csv(res_1, output_file, row.names = TRUE)

dreamlet::plotHeatmap(df_cts, genes = rownames(res_1))#, assays=colnames(df_cts)[2:3])
dev.off()

##### compare first two timepoints summarized #####
# Rownames to filter
filas_a_filtrar <- c("DHCR7", "HMGCS1", "MSMO1")

# Filter the df
res_1_bis <- res_1[rownames(res_1) %in% filas_a_filtrar, ]

dreamlet::plotHeatmap(df_cts, genes = rownames(res_1_bis))#, assays=colnames(df_cts)[2:3])
dev.off()

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### compare next two timepoints #####
ct.pairs_2 <- c("2_weeks-3_months", "<2_weeks")

# run comparison
fit_2 <- dreamletCompareClusters(pb, ct.pairs_2, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_2 <- topTable(fit_2, coef = "compare", number = 10)
res_all_2 <- topTable(fit_2, coef = "compare", number = Inf)

# Significant Dreamlet genes
genes_dreamlet_2 <- rownames(res_all_2[res_all_2$adj.P.Val < 0.05, ])

# Save gene list as txt, one gene per line
writeLines(genes_dreamlet_2, .output_path(.current_dir, "genes_dreamlet_2.txt"))

res_2

# Save
output_dir <- file.path(project_dir, "Resultados", "Dreamlet_Vs_PyDESeq2")
output_file <- .output_path(output_dir, "Dreamlet_de_genes_CD8_Comparison_2.csv")
write.csv(res_2, output_file, row.names = TRUE)

dreamlet::plotHeatmap(df_cts, genes = rownames(res_2))
dev.off()

############################################################################################################################################################################################################
############################################################################################################################################################################################################

##### compare last two timepoints #####
ct.pairs_3 <- c(">3_months", "2_weeks-3_months")

# run comparison
fit_3 <- dreamletCompareClusters(pb, ct.pairs_3, method = "none")

# Extract top 10 differentially expressed genes
# The coefficient 'compare' is the value logFC between test and baseline:
# compare = cellClustertest - cellClusterbaseline
res_3 <- topTable(fit_3, coef = "compare", number = 10)
res_all_3 <- topTable(fit_3, coef = "compare", number = Inf)

# Significant Dreamlet genes
genes_dreamlet_3 <- rownames(res_all_3[res_all_3$adj.P.Val < 0.05, ])

# Save gene list as txt, one gene per line
writeLines(genes_dreamlet_3, .output_path(.current_dir, "genes_dreamlet_3.txt"))

res_3

# Save
output_dir <- file.path(project_dir, "Resultados", "Dreamlet_Vs_PyDESeq2")
output_file <- .output_path(output_dir, "Dreamlet_de_genes_CD8_Comparison_3.csv")
write.csv(res_3, output_file, row.names = TRUE)

dreamlet::plotHeatmap(df_cts, genes = rownames(res_3))
dev.off()

############################################################################################################################################################################################################
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

.current_dir <- file.path(project_dir, "Input")
saveRDS(res.compare, .output_path(.current_dir, "res_compare.RDS"))
res.compare

plotVolcano(res.compare, coef = "compare", ncol = 4)
dev.off()

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
p1
dev.off()

data_plot <- p1$data
data_plot <- add_NTotalGenes(data_plot, go.gs)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
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
p2
dev.off()

data_plot <- p2$data
data_plot <- add_NTotalGenes(data_plot, go.gs)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
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
p3
dev.off()

data_plot <- p3$data
data_plot <- add_NTotalGenes(data_plot, go.gs)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
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
p4
dev.off()

data_plot <- p4$data
data_plot <- add_NTotalGenes(data_plot, go.gs)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
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

res.zenith = zenith_gsa( res.compare, gsc, 
              coef = 'compare', 
              n_genes_min=2)

res.zenith$assay = factor(res.zenith$assay, names(res.compare))

res.zenith <- add_NTotalGenes(res.zenith, gsc)

# View the updated res.zenith
head(res.zenith)

p5 <- plotZenithResults(res.zenith, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
p5
dev.off()

data_plot <- p5$data
data_plot <- add_NTotalGenes(data_plot, gsc)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
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
    plot.title = element_text(hjust = 0.5)
    ,
  )
dev.off()

Custom_genesets <- res.zenith


###################################################################################################################################################################################################
###################################################################################################################################################################################################
###################################################################################################################################################################################################
###################################################################################################################################################################################################
###################################################################################################################################################################################################


##### Filtered_Results for main Figure #####

## Filter GO_Biological_Process_2023
GO_Biological_Process_2023$Geneset

# Filter by
patrones <- c("ER231_", "ER781_")

# Put all in a variable with OR
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
GO_Biological_Process_2023_filtered <- GO_Biological_Process_2023[grep(patron_combinado, GO_Biological_Process_2023$Geneset), ]

###################################################################################################################################################################################################

## Filter Reactome_2022
Reactome_2022$Geneset

# Filter by
patrones <- c("ER1439_")

# Put all in a variable with OR
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
Reactome_2022_filtered <- Reactome_2022[grep(patron_combinado, Reactome_2022$Geneset), ]

###################################################################################################################################################################################################

## Filter Custom_genesets
Custom_genesets$Geneset

# Filter by
patrones <- c("Cytotoxicity")

# Put all in a variable with OR
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
Custom_genesets_filtered <- Custom_genesets[grep(patron_combinado, Custom_genesets$Geneset), ]

###################################################################################################################################################################################################

## Final filtered
Final_result <- rbind(GO_Biological_Process_2023_filtered, Reactome_2022_filtered, Custom_genesets_filtered)

Final_result_orig_main <- Final_result
head(Final_result)

.current_dir <- file.path(project_dir, "Input")
saveRDS(Final_result_orig_main, .output_path(.current_dir, "Final_result_orig_main.RDS"))

head(Final_result)

p6 <- plotZenithResults(Final_result, Inf, Inf, sortByGeneset = FALSE)

# Plot results, but with no limit based on the highest/lowest t-statistic
p6
dev.off()

data_plot <- p6$data
data_plot <- add_NTotalGenes(data_plot, gsc)
data_plot$GeneRatio <- data_plot$NGenes/data_plot$NTotalGenes
data_plot$logFDR <- -log10(data_plot$FDR)

# Horizontal dotplot
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


###################################################################################################################################################################################################
###################################################################################################################################################################################################
###################################################################################################################################################################################################
###################################################################################################################################################################################################
###################################################################################################################################################################################################


##### Filtered_Results for supplementary figure #####

## Filter GO_Biological_Process_2023
GO_Biological_Process_2023$Geneset

# Filter by
patrones <- c("ER856_", "ER4954_", "ER3628_", "ER2701_", "ER1707_")

# Put all in a variable with OR
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
GO_Biological_Process_2023_filtered <- GO_Biological_Process_2023[grep(patron_combinado, GO_Biological_Process_2023$Geneset), ]

###################################################################################################################################################################################################

## Filter KEGG_2021_Human
KEGG_2021_Human$Geneset

# Filter by
patrones <- c("ER68_", "ER229_", "ER111_")

# Put all in a variable with OR
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
KEGG_2021_Human_filtered <- KEGG_2021_Human[grep(patron_combinado, KEGG_2021_Human$Geneset), ]

###################################################################################################################################################################################################

## Filter Reactome_2022
Reactome_2022$Geneset

# Filter by
patrones <- c("ER1188_", "ER1412_", "ER712_")

# Put all in a variable with OR
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
Reactome_2022_filtered <- Reactome_2022[grep(patron_combinado, Reactome_2022$Geneset), ]

###################################################################################################################################################################################################

## Filter WikiPathway_2023_Human
WikiPathway_2023_Human$Geneset

# Filter by
patrones <- c("ER65_", "ER493_", "ER450_")

# Put all in a variable with OR
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
WikiPathway_2023_Human_filtered <- WikiPathway_2023_Human[grep(patron_combinado, WikiPathway_2023_Human$Geneset), ]

###################################################################################################################################################################################################

## Filter Custom_genesets
Custom_genesets$Geneset

# Filter by
patrones <- c("Genes Prolif")

# Put all in a variable with OR
patron_combinado <- paste(patrones, collapse = "|")

# Filter dataframe
Custom_genesets_filtered <- Custom_genesets[grep(patron_combinado, Custom_genesets$Geneset), ]

###################################################################################################################################################################################################

## Final filtered
Final_result <- rbind(GO_Biological_Process_2023_filtered, Reactome_2022_filtered, KEGG_2021_Human_filtered, WikiPathway_2023_Human_filtered, Custom_genesets_filtered)

Final_result_orig <- Final_result
head(Final_result)

.current_dir <- file.path(project_dir, "Input")
saveRDS(Final_result_orig, .output_path(.current_dir, "Final_result_orig.RDS"))

################################
######## END OF SCRIPT #########
################################
