##### Carga de librerias #####
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