###############################################################################
###############################################################################

# Program: 4.2_Integrated_datasets_Workflow_plus_Jordanas.R
# Author: Sergio Cámara Peña
# Date: 29/10/2024
# Version: FINAL

###############################################################################
###############################################################################

##### Load required libraries #####
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(cowplot)
library(patchwork)


##### Set seed #####
set.seed(2504)


##### Switches #####
Create_h5ad <- FALSE


##### Loading the different datasets #####
## Bai et al dataset
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Bai_et_al/RDS")
Seurat_list_Bai <- readRDS("Normalized_CellRanger_Bai_RDS.rds")
colnames(Seurat_list_Bai[[1]])
names(Seurat_list_Bai) <- c("Bai_CR_basal", "Bai_CR_CD19_cocult", "Bai_CR_MESO_cocult", "Bai_HD_CD19_cocult", "Bai_HD_MESO_cocult", "Bai_HD_basal", "Bai_NR_basal", "Bai_NR_CD19_cocult", "Bai_NR_MESO_cocult")

for (i in seq_along(Seurat_list_Bai)) {
    Seurat_list_Bai[[i]] <- RenameCells(Seurat_list_Bai[[i]], new.names = paste0(names(Seurat_list_Bai)[i], "_", colnames(Seurat_list_Bai[[i]])))
    Seurat_list_Bai[[i]]$Product_norm <- names(Seurat_list_Bai)[i]
}

colnames(Seurat_list_Bai[[1]])

# V5 - Remove CAR custom genes #
counts <- NULL
Seurat_list_Bai_V5 <- list()
for (i in seq_along(Seurat_list_Bai)){
    counts <- GetAssayData(Seurat_list_Bai[[i]], assay = "RNA")
    counts <- counts[-(which(rownames(counts) == "CAR")),]
    Seurat_list_Bai_V5[[i]] <- subset(Seurat_list_Bai[[i]], features = rownames(counts))
    counts <- NULL

    names(Seurat_list_Bai_V5)[i] <- names(Seurat_list_Bai)[i]

    print(paste0("DONE for ", names(Seurat_list_Bai_V5)[i]))
}

rm(Seurat_list_Bai)
Seurat_list_Bai_V5[[1]]@meta.data %>% head()


## Boroughs et al dataset
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Boroughs_et_al/RDS")
Seurat_list_Boroughs <- readRDS("Normalized_CellRanger_Boroughs_RDS.rds")

names(Seurat_list_Boroughs) <- c(
    "Bor_D1_28Z_CD19_Stim", "Bor_D1_28Z_No_Stim", "Bor_D1_BBZ_CD19_Stim", "Bor_D1_BBZ_No_Stim", "Bor_D1_Z_CD19_Stim", "Bor_D1_Z_No_Stim",
    "Bor_D2_28Z_CD19_Stim", "Bor_D2_28Z_No_Stim", "Bor_D2_BBZ_CD19_Stim", "Bor_D2_BBZ_No_Stim", "Bor_D2_Z_CD19_Stim", "Bor_D2_Z_No_Stim"
)

for (i in seq_along(Seurat_list_Boroughs)) {
    Seurat_list_Boroughs[[i]] <- RenameCells(Seurat_list_Boroughs[[i]], new.names = paste0(names(Seurat_list_Boroughs)[i], "_", colnames(Seurat_list_Boroughs[[i]])))
    Seurat_list_Boroughs[[i]]$Product_norm <- names(Seurat_list_Boroughs)[i]
}

colnames(Seurat_list_Boroughs[[1]])
Seurat_list_Boroughs[[1]]@meta.data %>% colnames()


## Deng et al dataset
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/RDS")
Seurat_list_Deng <- readRDS("Normalized_CellRanger_Deng_RDS.rds")

names(Seurat_list_Deng) <- c(
    "Den_Pt_14", "Den_Pt_15", "Den_Pt_16", "Den_Pt_18", "Den_Pt_20", "Den_Pt_21", "Den_Pt_26", "Den_Pt_27",
    "Den_Pt_28", "Den_Pt_33", "Den_Pt_34", "Den_Pt_37", "Den_Pt_38", "Den_Pt_40", "Den_Pt_41", "Den_Pt_42",
    "Den_Pt_43", "Den_Pt_49", "Den_Pt_50", "Den_Pt_54", "Den_Pt_55", "Den_Pt_56", "Den_Pt_59", "Den_Pt_64"
)

for (i in seq_along(Seurat_list_Deng)) {
    Seurat_list_Deng[[i]] <- RenameCells(Seurat_list_Deng[[i]], new.names = paste0(names(Seurat_list_Deng)[i], "_", colnames(Seurat_list_Deng[[i]])))
    Seurat_list_Deng[[i]]$Product_norm <- names(Seurat_list_Deng)[i]
}

colnames(Seurat_list_Deng[[1]])

# V5 - Remove CAR custom genes #
counts <- NULL
Seurat_list_Deng_V5 <- list()
for (i in seq_along(Seurat_list_Deng)){
    counts <- GetAssayData(Seurat_list_Deng[[i]], assay = "RNA")
    counts <- counts[-(which(rownames(counts) == "FMC63-28Z")),]
    Seurat_list_Deng_V5[[i]] <- subset(Seurat_list_Deng[[i]], features = rownames(counts))
    counts <- NULL

    names(Seurat_list_Deng_V5)[i] <- names(Seurat_list_Deng)[i]

    print(paste0("DONE for ", names(Seurat_list_Deng_V5)[i]))
}

rm(Seurat_list_Deng)


## Good et al dataset
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Good_et_al/RDS")
Seurat_list_Good <- readRDS("Normalized_CellRanger_Good_RDS.rds")

names(Seurat_list_Good) <- c("Goo_Pt_110", "Goo_Pt_116", "Goo_Pt_125", "Goo_Pt_129", "Goo_Pt_245", "Goo_Pt_253", "Goo_Pt_263", "Goo_Pt_276", "Goo_Pt_282")

for (i in seq_along(Seurat_list_Good)) {
    Seurat_list_Good[[i]] <- RenameCells(Seurat_list_Good[[i]], new.names = paste0(names(Seurat_list_Good)[i], "_", colnames(Seurat_list_Good[[i]])))
    Seurat_list_Good[[i]]$Product_norm <- names(Seurat_list_Good)[i]
}

colnames(Seurat_list_Good[[1]])

# V5 - Remove CAR custom genes #
counts <- NULL
Seurat_list_Good_V5 <- list()
for (i in seq_along(Seurat_list_Good)){
    counts <- GetAssayData(Seurat_list_Good[[i]], assay = "RNA")
    counts <- counts[-(which(rownames(counts) %in% c("axicel", "CAR19-28Z-Mackall"))),]
    Seurat_list_Good_V5[[i]] <- subset(Seurat_list_Good[[i]], features = rownames(counts))
    counts <- NULL

    names(Seurat_list_Good_V5)[i] <- names(Seurat_list_Good)[i]

    print(paste0("DONE for ", names(Seurat_list_Good_V5)[i]))
}

rm(Seurat_list_Good)


## Lynn et al dataset
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Lynn_et_al/RDS")
Seurat_list_Lynn <- readRDS("Normalized_CellRanger_Lynn_RDS.rds")

names(Seurat_list_Lynn) <- c("Lyn_Exp1_CD19", "Lyn_Exp1_GD2", "Lyn_Exp2_Cont", "Lyn_Exp2_JUN")

for (i in seq_along(Seurat_list_Lynn)) {
    Seurat_list_Lynn[[i]] <- RenameCells(Seurat_list_Lynn[[i]], new.names = paste0(names(Seurat_list_Lynn)[i], "_", colnames(Seurat_list_Lynn[[i]])))
    Seurat_list_Lynn[[i]]$Product_norm <- names(Seurat_list_Lynn)[i]
}

colnames(Seurat_list_Lynn[[1]])


## Melenhorst et al dataset
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Melenhorst_et_al/RDS")
Seurat_list_Melenhorst <- readRDS("Normalized_CellRanger_Melenhorst_RDS.rds")

names(Seurat_list_Melenhorst) <- c("Mel_PT1_M12", "Mel_PT1_M15", "Mel_PT1_Y9", "Mel_PT2_M3", "Mel_PT2_Y3", "Mel_PT2_Y6_5")

for (i in seq_along(Seurat_list_Melenhorst)) {
    Seurat_list_Melenhorst[[i]] <- RenameCells(Seurat_list_Melenhorst[[i]], new.names = paste0(names(Seurat_list_Melenhorst)[i], "_", colnames(Seurat_list_Melenhorst[[i]])))
    Seurat_list_Melenhorst[[i]]$Product_norm <- names(Seurat_list_Melenhorst)[i]
}

colnames(Seurat_list_Melenhorst[[1]])

# V5 - Remove CAR custom genes #
counts <- NULL
Seurat_list_Melenhorst_V5 <- list()
for (i in seq_along(Seurat_list_Melenhorst)){
    counts <- GetAssayData(Seurat_list_Melenhorst[[i]], assay = "RNA")
    counts <- counts[-(which(rownames(counts) == "muCD19BBz")),]
    Seurat_list_Melenhorst_V5[[i]] <- subset(Seurat_list_Melenhorst[[i]], features = rownames(counts))
    counts <- NULL

    names(Seurat_list_Melenhorst_V5)[i] <- names(Seurat_list_Melenhorst)[i]

    print(paste0("DONE for ", names(Seurat_list_Melenhorst_V5)[i]))
}

rm(Seurat_list_Melenhorst)


## Sheih et al dataset
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Sheih_et_al/RDS")
Seurat_list_Sheih <- readRDS("Normalized_CellRanger_Sheih_RDS.rds")

names(Seurat_list_Sheih) <- c(
    "She_CLL_1_d21", "She_CLL_1_d38", "She_CLL_1_d112", "She_CLL_1_IP", "She_CLL_2_d12", "She_CLL_2_d29", "She_CLL_2_d83", "She_CLL_2_IP",
    "She_NHL_6_d12", "She_NHL_6_d29", "She_NHL_6_d102", "She_NHL_7_d12", "She_NHL_7_d28", "She_NHL_7_d89", "She_NHL_7_IP"
)

for (i in seq_along(Seurat_list_Sheih)) {
    Seurat_list_Sheih[[i]] <- RenameCells(Seurat_list_Sheih[[i]], new.names = paste0(names(Seurat_list_Sheih)[i], "_", colnames(Seurat_list_Sheih[[i]])))
    Seurat_list_Sheih[[i]]$Product_norm <- names(Seurat_list_Sheih)[i]
}

colnames(Seurat_list_Sheih[[1]])


## Wang et al dataset
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/RDS")
Seurat_list_Wang <- readRDS("Normalized_CellRanger_Wang_RDS.rds")

names(Seurat_list_Wang) <- c("Wan_PD1", "Wan_PD2", "Wan_PD3", "Wan_SPD1", "Wan_SPD2", "Wan_SPD3")

for (i in seq_along(Seurat_list_Wang)) {
    Seurat_list_Wang[[i]] <- RenameCells(Seurat_list_Wang[[i]], new.names = paste0(names(Seurat_list_Wang)[i], "_", colnames(Seurat_list_Wang[[i]])))
    Seurat_list_Wang[[i]]$Product_norm <- names(Seurat_list_Wang)[i]
}

colnames(Seurat_list_Wang[[1]])

# V5 - Remove CAR custom genes #
counts <- NULL
Seurat_list_Wang_V5 <- list()
for (i in seq_along(Seurat_list_Wang)){
    counts <- GetAssayData(Seurat_list_Wang[[i]], assay = "RNA")
    counts <- counts[-(which(rownames(counts) == "IFNB1-5LTR")),]
    Seurat_list_Wang_V5[[i]] <- subset(Seurat_list_Wang[[i]], features = rownames(counts))
    counts <- NULL

    names(Seurat_list_Wang_V5)[i] <- names(Seurat_list_Wang)[i]

    print(paste0("DONE for ", names(Seurat_list_Wang_V5)[i]))
}

rm(Seurat_list_Wang)


## Xhangolli et al dataset
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/RDS")
Seurat_list_Xhangolli <- readRDS("Normalized_CellRanger_Xhangolli_RDS.rds")

for (i in seq_along(Seurat_list_Xhangolli)) {
    Seurat_list_Xhangolli[[i]] <- RenameCells(Seurat_list_Xhangolli[[i]], new.names = paste0("Xha_", colnames(Seurat_list_Xhangolli[[i]])))
    Seurat_list_Xhangolli[[i]]$Product_norm <- paste0("Xha_", names(Seurat_list_Xhangolli)[i])
}

colnames(Seurat_list_Xhangolli[[1]])


## Li X. et al dataset
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Li_X_et_al/RDS")
Seurat_list_Li_X <- readRDS("Normalized_CellRanger_Li_RDS.rds")

names(Seurat_list_Li_X) <- c("LiX_IP", "LiX_PP", "LiX_RP")

for (i in seq_along(Seurat_list_Li_X)) {
    Seurat_list_Li_X[[i]] <- RenameCells(Seurat_list_Li_X[[i]], new.names = paste0(names(Seurat_list_Li_X)[i], "_", colnames(Seurat_list_Li_X[[i]])))
    Seurat_list_Li_X[[i]]$Product_norm <- names(Seurat_list_Li_X)[i]
}

colnames(Seurat_list_Li_X[[1]])


## Rodriguez-Marquez et al dataset
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Rodriguez-Marquez_et_al/RDS")
Seurat_list_Rodriguez_Marquez <- readRDS("Normalized_CellRanger_Rodrguez_Marquez_RDS.rds")

a <- names(Seurat_list_Rodriguez_Marquez)

names(Seurat_list_Rodriguez_Marquez) <- c("Rod_D10", "Rod_D10_d7", "Rod_D14", "Rod_D14_d7", "Rod_D18", "Rod_D18_d7")

for (i in seq_along(Seurat_list_Rodriguez_Marquez)) {
    Seurat_list_Rodriguez_Marquez[[i]] <- RenameCells(Seurat_list_Rodriguez_Marquez[[i]], new.names = paste0(names(Seurat_list_Rodriguez_Marquez)[i], "_", colnames(Seurat_list_Rodriguez_Marquez[[i]])))
    Seurat_list_Rodriguez_Marquez[[i]]$Product_norm <- names(Seurat_list_Rodriguez_Marquez)[i]
}

colnames(Seurat_list_Rodriguez_Marquez[[1]])


## Haradvala et al dataset
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/RDS")
Seurat_list_Haradvala <- readRDS("Normalized_CellRanger_Haradvala_RDS.rds")

names(Seurat_list_Haradvala)

for (i in seq_along(Seurat_list_Haradvala)) {
    Seurat_list_Haradvala[[i]] <- RenameCells(Seurat_list_Haradvala[[i]], new.names = paste0(names(Seurat_list_Haradvala)[i], "_", colnames(Seurat_list_Haradvala[[i]])))
    Seurat_list_Haradvala[[i]]$Product_norm <- names(Seurat_list_Haradvala)[i]
}

colnames(Seurat_list_Haradvala[[1]])

# V5 - Remove CAR custom genes #
counts <- NULL
Seurat_list_Haradvala_V5 <- list()
for (i in seq_along(Seurat_list_Haradvala)){
    counts <- GetAssayData(Seurat_list_Haradvala[[i]], assay = "RNA")
    counts <- counts[-(which(rownames(counts) %in% c("Yescarta", "Kymriah"))),]
    Seurat_list_Haradvala_V5[[i]] <- subset(Seurat_list_Haradvala[[i]], features = rownames(counts))
    counts <- NULL

    names(Seurat_list_Haradvala_V5)[i] <- names(Seurat_list_Haradvala)[i]

    print(paste0("DONE for ", names(Seurat_list_Haradvala_V5)[i]))
}

rm(Seurat_list_Haradvala)


## Li X. Cancer cell letter et al dataset
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Li_X_Cancer_Cell_letter_et_al/RDS")
Seurat_list_Li_X_letter <- readRDS("Normalized_CellRanger_Li_X_letter_RDS.rds")

names(Seurat_list_Li_X_letter)

for (i in seq_along(Seurat_list_Li_X_letter)) {
    Seurat_list_Li_X_letter[[i]] <- RenameCells(Seurat_list_Li_X_letter[[i]], new.names = paste0(names(Seurat_list_Li_X_letter)[i], "_", colnames(Seurat_list_Li_X_letter[[i]])))
    Seurat_list_Li_X_letter[[i]]$Product_norm <- names(Seurat_list_Li_X_letter)[i]
}

colnames(Seurat_list_Li_X_letter[[1]])

# V5 - Remove CAR custom genes #
counts <- NULL
Seurat_list_Li_X_letter_V5 <- list()
for (i in seq_along(Seurat_list_Li_X_letter)){
    counts <- GetAssayData(Seurat_list_Li_X_letter[[i]], assay = "RNA")
    counts <- counts[-(which(rownames(counts) %in% c("Yescarta", "Kymriah"))),]
    Seurat_list_Li_X_letter_V5[[i]] <- subset(Seurat_list_Li_X_letter[[i]], features = rownames(counts))
    counts <- NULL

    names(Seurat_list_Li_X_letter_V5)[i] <- names(Seurat_list_Li_X_letter)[i]

    print(paste0("DONE for ", names(Seurat_list_Li_X_letter_V5)[i]))
}

rm(Seurat_list_Li_X_letter)


## Jordana et al dataset
setwd("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Jordana_et_al/RDS")
Seurat_Jordana <- readRDS("Normalized_CellRanger_Jordana_RDS.rds")

Seurat_Jordana@meta.data


####################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################

### Merge ---- Worst case scenario
Bai_merged <- merge(
    x = Seurat_list_Bai_V5[[1]],
    y = c(Seurat_list_Bai_V5[[2]], Seurat_list_Bai_V5[[3]], Seurat_list_Bai_V5[[4]], Seurat_list_Bai_V5[[5]], Seurat_list_Bai_V5[[6]], Seurat_list_Bai_V5[[7]], Seurat_list_Bai_V5[[8]], Seurat_list_Bai_V5[[9]])
)

Seurat_merged <- merge(x = Bai_merged, y = c(
    Seurat_list_Boroughs, Seurat_list_Deng_V5, Seurat_list_Good_V5, Seurat_list_Lynn, Seurat_list_Melenhorst_V5,
    Seurat_list_Sheih, Seurat_list_Wang_V5, Seurat_list_Xhangolli, Seurat_list_Li_X, Seurat_list_Rodriguez_Marquez,
    Seurat_list_Haradvala_V5, Seurat_list_Li_X_letter_V5, Seurat_Jordana
))

# Seurat_Jordana # +47855 cells
Seurat_merged # 491997 cells

## Remove NON-CD3+ population in all datasets
Seurat_merged <- subset(x = Seurat_merged, subset = CD3D > 0)
Seurat_merged # 455628 cells

## Keep doing the normal workflow
Seurat_merged <- FindVariableFeatures(Seurat_merged,
    selection.method = "vst",
    nfeatures = 2000
)

Seurat_merged <- ScaleData(Seurat_merged)

# Run PCA
Seurat_merged <- RunPCA(object = Seurat_merged)

# Elbow plot to pick up only relevant dimensions
ElbowPlot(Seurat_merged, ndims = 50, reduction = "pca")
dev.off()

# Run UMAP
Seurat_merged <- RunUMAP(Seurat_merged,
    dims = 1:30,
    reduction = "pca",
    reduction.name = "umap_wo_integ"
)

colors_palette <- c("#CDFD02", "#6C35FF", "#FBA3B4", "#CBC7FE", "#2B097E", "#FE4F11", "#8E0D0D", "#A9FCFA", "#D35E82", "#CECECE", "#FFC49C", "#383838", "#00A57E", "#7D8570")

# Plot UMAP
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V5/WO_integ/WO_integ_Seurat.pdf")
DimPlot(Seurat_merged,
    group.by = "orig.ident",
    reduction = "umap_wo_integ",
    raster = FALSE,
    cols = colors_palette
)
dev.off()

# Quality control metrics
p <- FeaturePlot(Seurat_merged,
    reduction = "umap_wo_integ",
    features = c("nCount_RNA", "nFeature_RNA", "log10GenesPerUMI", "mitoRatio"),
    order = TRUE,
    raster = FALSE
)

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V5/WO_integ/QC_metrics_WO_integ_Seurat.pdf")
p + plot_annotation(title = "All datasets merged", theme = theme(plot.title = element_text(size = 16)))
dev.off()

rm(p)

# Cell cycle plot
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V5/WO_integ/WO_integ_Seurat_cellcycle_scoring.pdf")
DimPlot(Seurat_merged,
    group.by = "Phase",
    reduction = "umap_wo_integ",
    raster = FALSE,
)

DimPlot(Seurat_merged,
    group.by = "Phase",
    reduction = "pca",
    raster = FALSE,
)
dev.off()

saveRDS(Seurat_merged, "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/RDS/WO_integ/V5/Seurat_merged_WO_integ.RDS")

## Add cell type information with Celltypist
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Python-Celltypist/V5")

## Create or load h5ad file if needed
if (Create_h5ad) {
    # Here I remove unnecesary slots that produce some errors in Python
    Seurat_merged2 <- CreateSeuratObject(Seurat_merged$RNA@counts)
    Seurat_merged2@meta.data <- Seurat_merged@meta.data
    SaveH5Seurat(Seurat_merged2, filename = "Seurat_merged.h5Seurat")
    Convert("Seurat_merged.h5Seurat", dest = "h5ad")
    print("Correctamente guardado")
} else {
    load_Seurat_merged <- LoadH5Seurat("Seurat_merged.h5Seurat")
    Seurat_merged <- load_Seurat_merged
    rm(load_Seurat_merged)
    print("Correctamente cargado")
}

############## PAUSE ##############
# Here, you should switch to Python ("4.3_Celltypist_annotation_plus_Jordanas.py") to obtain the CellTypist results
###################################


## Load and add metadata obtained with CellTypist
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Python-Celltypist/V5")
celltypist_metadata <- read.csv(file = "celltypist_metadata_table.csv", row.names = 1)
celltypist_metadata <- celltypist_metadata %>%
    dplyr::select("predicted_labels", "over_clustering", "majority_voting", "conf_score") %>%
    rename_with(~ paste0("celltypist_", .), everything())

Seurat_merged <- readRDS("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/RDS/WO_integ/V5/Seurat_merged_WO_integ.RDS")
Seurat_merged <- AddMetaData(Seurat_merged, metadata = celltypist_metadata)

## Celltypist annotations plot
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Plots/V5/WO_integ/Celltypist_ann_WO_integ_Seurat.pdf")
DimPlot(Seurat_merged, reduction = "umap_wo_integ", group.by = "celltypist_majority_voting", pt.size = 0.3, raster = FALSE)
dev.off()

## Save again h5ad object to include celltypist data
Seurat_merged3 <- CreateSeuratObject(Seurat_merged$RNA@counts)
Seurat_merged3@meta.data <- Seurat_merged@meta.data
SaveH5Seurat(Seurat_merged3, filename = "Seurat_merged_With_Celltypist.h5Seurat")
Convert("Seurat_merged_With_Celltypist.h5Seurat", dest = "h5ad")

## Save state 2 in RDS
saveRDS(Seurat_merged, "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/RDS/WO_integ/V5/Seurat_merged_WO_integ_state2.RDS")

##### END OF SCRIPT #####