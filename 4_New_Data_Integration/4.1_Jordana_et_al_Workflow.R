###############################################################################
###############################################################################

# Program: 4.1_Jordana_et_al_Workflow.R
# Author: Sergio Cámara Peña
# Date: 2024
# Version: FINAL

###############################################################################
###############################################################################

##### Load required libraries #####
library(Seurat)
library(tidyverse)
library(cowplot)
library(scales)
library(patchwork)


##### Set seed #####
set.seed(2504)


##### Load data #####
setwd("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/New_Datasets/Jordana_et_al/Datos_lorea")
Jordana_Seurat <- readRDS("seurat_new.rds")


##### Metadata columns to leave #####
# [1] "orig.ident"       "nCount_RNA"       "nFeature_RNA"     "Product"
# [5] "log10GenesPerUMI" "mitoRatio"        "S.Score"          "G2M.Score"
# [9] "Phase"            "Product_norm"


##### Clean metadata #####
# Remove sample with 1 CAR
Jordana_Seurat
Jordana_Seurat <- subset(Jordana_Seurat, subset = orig.ident != "P2_PB_M3")
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

# Remove additional columns from the metadata of the Seurat object
Jordana_Seurat@meta.data <- Jordana_Seurat@meta.data %>% select(-ID)

Jordana_Seurat@meta.data$Product_norm <- paste0("Jor_", Jordana_Seurat@meta.data$Product)

Jordana_Seurat@meta.data %>% head()

# Re-order columns from the metadata of the Seurat object
Jordana_Seurat@meta.data <- Jordana_Seurat@meta.data %>%
    select(
        orig.ident, nCount_RNA, nFeature_RNA, Product, log10GenesPerUMI, mitoRatio, S.Score, G2M.Score, Phase, Product_norm,
        everything()
    ) # 'everything()' keeps remaining columns

Jordana_Seurat@meta.data %>% head()
Jordana_Seurat@meta.data %>% tail()

Jordana_Seurat@meta.data["Product_norm"] %>% unique()

##### Save RDS object state 1 #####
setwd("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Jordana_et_al/RDS")
saveRDS(Jordana_Seurat, "PostQC_CellRanger_Jordana_RDS.rds")


##### Normalization #####
setwd("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Jordana_et_al/RDS")
Jordana_Seurat <- readRDS("PostQC_CellRanger_Jordana_RDS.rds")

Jordana_Seurat <- NormalizeData(Jordana_Seurat)
Jordana_Seurat <- FindVariableFeatures(Jordana_Seurat,
    selection.method = "vst",
    nfeatures = 2000)
Jordana_Seurat <- ScaleData(Jordana_Seurat)

saveRDS(Jordana_Seurat, file = "Normalized_CellRanger_Jordana_RDS.rds")


##### Merge ---- Worst case scenario #####
setwd("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Jordana_et_al/RDS")
Jordana_Seurat <- readRDS("Normalized_CellRanger_Jordana_RDS.rds")

# Run PCA
Jordana_Seurat <- RunPCA(object = Jordana_Seurat, reduction.name = "pca_wo_integ")

# Run UMAP
Jordana_Seurat <- RunUMAP(Jordana_Seurat,
    dims = 1:30,
    reduction = "pca_wo_integ",
    reduction.name = "umap_wo_integ"
)

# Plot UMAP
pdf("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Jordana_et_al/Plots/Worst_case_scenario/WO_integ_Seurat.pdf")
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

pdf("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Jordana_et_al/Plots/Worst_case_scenario/QC_metrics_WO_integ_Seurat.pdf")
p + plot_annotation(title = paste0(Jordana_Seurat@meta.data$orig.ident %>% unique()), theme = theme(plot.title = element_text(size = 16)))
dev.off()

rm(p)

##### END OF SCRIPT #####
