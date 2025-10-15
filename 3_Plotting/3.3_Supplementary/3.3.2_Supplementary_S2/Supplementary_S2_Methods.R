###############################################################################
###############################################################################

# Program: Supplementary_S2_Methods.R
# Author: Sergio Cámara Peña
# Date: 09/06/2025
# Version: FINAL

###############################################################################
###############################################################################

##### Loading needed libraries #####
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(cowplot)
library(patchwork)
library(harmony)
library(rliger)
library(SeuratWrappers)
library(future)
library(STACAS)
library(dittoSeq)
library(reticulate)
library(sceasy)

##### Set seed #####
set.seed(2504)

##### Colour palette #####
palette <- c(
    "#F8766D",
    "#00BA38",
    "#FF9900",
    "#619CFF",
    "#F564E3",
    "#B79F00"
)

##### Switches #####
Primera_vez <- FALSE

##### Loading the different datasets
if (Primera_vez) {
    ## Xhangolli et al dataset
    setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/RDS")
    Seurat_list_Xhangolli <- readRDS("Normalized_CellRanger_Xhangolli_RDS.rds")

    for (i in seq_along(Seurat_list_Xhangolli)) {
        Seurat_list_Xhangolli[[i]] <- RenameCells(Seurat_list_Xhangolli[[i]], new.names = paste0("Xha_", colnames(Seurat_list_Xhangolli[[i]])))
        Seurat_list_Xhangolli[[i]]$Product_norm <- paste0("Xha_", names(Seurat_list_Xhangolli)[i])
    }

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

    ## Wang et al dataset
    setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/RDS")
    Seurat_list_Wang <- readRDS("Normalized_CellRanger_Wang_RDS.rds")

    names(Seurat_list_Wang) <- c("Wan_PD1", "Wan_PD2", "Wan_PD3", "Wan_SPD1", "Wan_SPD2", "Wan_SPD3")

    for (i in seq_along(Seurat_list_Wang)) {
        Seurat_list_Wang[[i]] <- RenameCells(Seurat_list_Wang[[i]], new.names = paste0(names(Seurat_list_Wang)[i], "_", colnames(Seurat_list_Wang[[i]])))
        Seurat_list_Wang[[i]]$Product_norm <- names(Seurat_list_Wang)[i]
    }

    colnames(Seurat_list_Wang[[1]])

    ## Lynn et al dataset
    setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Lynn_et_al/RDS")
    Seurat_list_Lynn <- readRDS("Normalized_CellRanger_Lynn_RDS.rds")

    names(Seurat_list_Lynn) <- c("Lyn_Exp1_CD19", "Lyn_Exp1_GD2", "Lyn_Exp2_Cont", "Lyn_Exp2_JUN")

    for (i in seq_along(Seurat_list_Lynn)) {
        Seurat_list_Lynn[[i]] <- RenameCells(Seurat_list_Lynn[[i]], new.names = paste0(names(Seurat_list_Lynn)[i], "_", colnames(Seurat_list_Lynn[[i]])))
        Seurat_list_Lynn[[i]]$Product_norm <- names(Seurat_list_Lynn)[i]
    }

    colnames(Seurat_list_Lynn[[1]])

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

    ## Bai et al dataset --- Keep only healthy
    setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Bai_et_al/RDS")
    Seurat_list_Bai <- readRDS("Normalized_CellRanger_Bai_RDS.rds")
    Seurat_list_Bai <- c(Seurat_list_Bai$Healthy_donor_3T3_CD19_coculture, Seurat_list_Bai$Healthy_donor_3T3_MESO_coculture, Seurat_list_Bai$Healthy_donor_Basal_CAR_T)
    # names(Seurat_list_Bai) <- c("Healthy_donor_3T3_CD19_coculture", "Healthy_donor_3T3_MESO_coculture", "Healthy_donor_Basal_CAR_T")

    names(Seurat_list_Bai) <- c("Bai_HD_CD19_cocult", "Bai_HD_MESO_cocult", "Bai_HD_basal")

    for (i in seq_along(Seurat_list_Bai)) {
        Seurat_list_Bai[[i]] <- RenameCells(Seurat_list_Bai[[i]], new.names = paste0(names(Seurat_list_Bai)[i], "_", colnames(Seurat_list_Bai[[i]])))
        Seurat_list_Bai[[i]]$Product_norm <- names(Seurat_list_Bai)[i]
    }

    colnames(Seurat_list_Bai[[1]])

    ##### Trying different integration methods #####
    ### Merge ---- Worst case scenario
    Rodriguez_marquez_merged <- merge(
        x = Seurat_list_Rodriguez_Marquez[[1]],
        y = c(Seurat_list_Rodriguez_Marquez[[2]], Seurat_list_Rodriguez_Marquez[[3]], Seurat_list_Rodriguez_Marquez[[4]], Seurat_list_Rodriguez_Marquez[[5]], Seurat_list_Rodriguez_Marquez[[6]])
    )

    Seurat_merged <- merge(x = Rodriguez_marquez_merged, y = c(Seurat_list_Xhangolli, Seurat_list_Wang, Seurat_list_Lynn, Seurat_list_Boroughs, Seurat_list_Bai))

    if (FALSE) {
        # Save version with raw counts
        SaveH5Seurat(
            Seurat_merged,
            filename = "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/scVI/Sin_GT_With_Python/Seurat_merged_RAW_for_Py.h5Seurat"
        )
        Convert(
            "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/scVI/Sin_GT_With_Python/Seurat_merged_RAW_for_Py.h5Seurat",
            dest = "h5ad",
            overwrite = TRUE
        )
    }

    # Read CSV with annotations
    annotations <- read.csv("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/cell_annotation_manual.csv")

    # Check cell_name is rowname
    rownames(annotations) <- annotations$cell_name

    # Filter Seurat object
    Seurat_filtered <- subset(Seurat_merged, cells = annotations$cell_name)

    annotations$cell_name <- NULL

    Seurat_filtered <- AddMetaData(Seurat_filtered, metadata = annotations)
    Seurat_merged <- Seurat_filtered

    is.na(Seurat_merged$manual_celltype_annotation_high) %>% sum()

    Seurat_merged <- FindVariableFeatures(Seurat_merged,
        selection.method = "vst",
        nfeatures = 2000
    )

    Seurat_merged <- ScaleData(Seurat_merged)

    # Run PCA
    Seurat_merged <- RunPCA(object = Seurat_merged, reduction.name = "pca_wo_integ")

    # Run UMAP
    Seurat_merged <- RunUMAP(Seurat_merged,
        dims = 1:30,
        reduction = "pca_wo_integ",
        reduction.name = "umap_wo_integ"
    )

    # Plot UMAP
    pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/Merged_WO_integration/Sin_GT/WO_integ_Seurat.pdf")
    DimPlot(Seurat_merged,
        group.by = "orig.ident",
        reduction = "umap_wo_integ",
        raster = FALSE
    )
    dev.off()

    # Quality control metrics
    p <- FeaturePlot(Seurat_merged,
        reduction = "umap_wo_integ",
        features = c("nCount_RNA", "nFeature_RNA", "log10GenesPerUMI", "mitoRatio"),
        order = TRUE,
        raster = FALSE
    )

    pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/Merged_WO_integration/Sin_GT/QC_metrics_WO_integ_Seurat.pdf")
    p + plot_annotation(title = "Healthy data merged", theme = theme(plot.title = element_text(size = 16)))
    dev.off()

    rm(p)


    # Set path
    path_base <- "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/Merged_WO_integration/Sin_GT/Seurat_merged"

    # Save .h5Seurat
    SaveH5Seurat(Seurat_merged, filename = paste0(path_base, ".h5Seurat"))

    # Transform into .h5ad
    Convert(paste0(path_base, ".h5Seurat"), dest = "h5ad", overwrite = TRUE)

    saveRDS(Seurat_merged, "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/Merged_WO_integration/Sin_GT/Seurat_merged.RDS")
} else {
    setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/Merged_WO_integration/Sin_GT")
    Seurat_merged <- readRDS("Seurat_merged.RDS")
}

pdf("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Suplementarias/Suplementaria_S2_1_WO_integ_Seurat.pdf")
DimPlot(Seurat_merged, reduction = "umap_wo_integ", group.by = "orig.ident", pt.size = 0.3, raster = FALSE, cols = palette) +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
    )
dev.off()

#######################################################################################################################################################################################
#######################################################################################################################################################################################

### scVI is done in python script ("Supplementary_S2_scVI.py")

#######################################################################################################################################################################################
#######################################################################################################################################################################################

### Harmony
if (Primera_vez) {
    Seurat_harmony <- Seurat_merged
    Seurat_harmony <- RunHarmony(Seurat_harmony, "Product_norm",
        reduction = "pca_wo_integ",
        dims.use = 1:30, reduction.save = "harmony"
    )

    # Run UMAP
    Seurat_harmony <- RunUMAP(Seurat_harmony,
        dims = 1:30,
        reduction = "harmony",
        reduction.name = "harmony_umap"
    )

    Seurat_harmony <- FindNeighbors(Seurat_harmony, reduction = "harmony", dims = 1:30)

    Seurat_harmony <- FindClusters(Seurat_harmony)

    saveRDS(Seurat_harmony, "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/Harmony/Sin_GT/Seurat_harmony.RDS")

    # Set path
    path_base <- "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/Harmony/Sin_GT/Seurat_harmony"

    # Save .h5Seurat
    SaveH5Seurat(Seurat_harmony, filename = paste0(path_base, ".h5Seurat"))

    # Transform into .h5ad
    Convert(paste0(path_base, ".h5Seurat"), dest = "h5ad", overwrite = TRUE)
} else {
    Seurat_harmony <- readRDS("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/Harmony/Sin_GT/Seurat_harmony.RDS")
}

pdf("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Suplementarias/Suplementaria_S2_3_Harmony.pdf")
DimPlot(Seurat_harmony, reduction = "harmony_umap", group.by = "orig.ident", pt.size = 0.3, raster = FALSE, cols = palette) +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
    )
dev.off()

#######################################################################################################################################################################################
#######################################################################################################################################################################################

### LIGER
if (Primera_vez) {
    Seurat_liger <- Seurat_merged

    Seurat_liger <- ScaleData(Seurat_liger, split.by = "Product_norm", do.center = FALSE)
    Seurat_liger@meta.data$Marked_cells <- NULL

    Seurat_liger <- RunOptimizeALS(Seurat_liger, k = 20, lambda = 5, split.by = "Product_norm")
    Seurat_liger <- RunQuantileNorm(Seurat_liger, split.by = "Product_norm")

    Seurat_liger <- FindNeighbors(Seurat_liger, reduction = "iNMF", dims = 1:20)
    Seurat_liger <- FindClusters(Seurat_liger)

    Seurat_liger <- RunUMAP(Seurat_liger, dims = 1:ncol(Seurat_liger[["iNMF"]]), reduction = "iNMF", reduction.name = "liger_umap")

    saveRDS(Seurat_liger, "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/LIGER/Sin_GT/Seurat_liger.RDS")

    # Set path
    path_base <- "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/LIGER/Sin_GT/Seurat_liger"

    # Save .h5Seurat
    SaveH5Seurat(Seurat_liger, filename = paste0(path_base, ".h5Seurat"))

    # Transform into .h5ad
    Convert(paste0(path_base, ".h5Seurat"), dest = "h5ad", overwrite = TRUE)
} else {
    Seurat_liger <- readRDS("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/LIGER/Sin_GT/Seurat_liger.RDS")
}

pdf("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Suplementarias/Suplementaria_S2_4_LIGER.pdf")
DimPlot(Seurat_liger, reduction = "liger_umap", group.by = "orig.ident", pt.size = 0.3, raster = FALSE, cols = palette) +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
    )
dev.off()

#######################################################################################################################################################################################
#######################################################################################################################################################################################

### STACAS
if (Primera_vez) {
    Seurat_STACAS <- SplitObject(Seurat_merged, split.by = "Product_norm")

    Seurat_STACAS_integ <- Seurat_STACAS %>%
        Run.STACAS(dims = 1:30, anchor.features = 2000) %>%
        RunUMAP(dims = 1:30)

    Seurat_STACAS_integ2 <- Seurat_STACAS_integ

    Seurat_STACAS_integ2 <- FindNeighbors(Seurat_STACAS_integ2, dims = 1:30)
    Seurat_STACAS_integ2 <- FindClusters(Seurat_STACAS_integ2)

    saveRDS(Seurat_STACAS_integ2, "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/STACAS/Sin_GT/Seurat_STACAS_integ2.RDS")

    # Set path
    path_base <- "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/STACAS/Sin_GT/Seurat_STACAS_integ2"

    # Save .h5Seurat
    SaveH5Seurat(Seurat_STACAS_integ2, filename = paste0(path_base, ".h5Seurat"))

    # Transform into .h5ad
    Convert(paste0(path_base, ".h5Seurat"), dest = "h5ad", overwrite = TRUE)
} else {
    Seurat_STACAS_integ2 <- readRDS("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/STACAS/Sin_GT/Seurat_STACAS_integ2.RDS")
}

pdf("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Suplementarias/Suplementaria_S2_5_STACAS.pdf")
DimPlot(Seurat_STACAS_integ2, group.by = "orig.ident", pt.size = 0.3, raster = FALSE, cols = palette) +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
    )
dev.off()

#######################################################################################################################################################################################
#######################################################################################################################################################################################

### Seurat (RPCA)
if (Primera_vez) {
    Seurat_RPCA <- SplitObject(Seurat_merged, split.by = "Product_norm")

    features <- SelectIntegrationFeatures(object.list = Seurat_RPCA)
    Seurat_RPCA <- lapply(X = Seurat_RPCA, FUN = function(x) {
        x <- ScaleData(x, features = features)
        x <- RunPCA(x, features = features)
    })

    RPCA_integ_anch <- FindIntegrationAnchors(
        object.list = Seurat_RPCA,
        anchor.features = features, reduction = "rpca"
    ) # You can increase the strength of alignment by increasing the k.anchor parameter, which is set to 5 by default.

    Seurat_RPCA_integ <- IntegrateData(anchorset = RPCA_integ_anch)

    DefaultAssay(Seurat_RPCA_integ) <- "integrated"

    Seurat_RPCA_integ <- ScaleData(Seurat_RPCA_integ)
    Seurat_RPCA_integ <- RunPCA(Seurat_RPCA_integ, reduction.name = "RPCA_pca")

    Seurat_RPCA_integ <- RunUMAP(Seurat_RPCA_integ, reduction = "RPCA_pca", dims = 1:30, reduction.name = "RPCA_umap")
    Seurat_RPCA_integ <- FindNeighbors(Seurat_RPCA_integ, reduction = "RPCA_pca", dims = 1:30)
    Seurat_RPCA_integ <- FindClusters(Seurat_RPCA_integ)
    saveRDS(Seurat_RPCA_integ, "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/Seurat_RPCA/Sin_GT/Seurat_RPCA_integ.RDS")

    # Set path
    path_base <- "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/Seurat_RPCA/Sin_GT/Seurat_RPCA_integ"

    # Save .h5Seurat
    SaveH5Seurat(Seurat_RPCA_integ, filename = paste0(path_base, ".h5Seurat"))

    # Transform into .h5ad
    Convert(paste0(path_base, ".h5Seurat"), dest = "h5ad", overwrite = TRUE)
} else {
    Seurat_RPCA_integ <- readRDS("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/Seurat_RPCA/Sin_GT/Seurat_RPCA_integ.RDS")
}

pdf("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Suplementarias/Suplementaria_S2_6_RPCA_Seurat.pdf")
DimPlot(Seurat_RPCA_integ, reduction = "RPCA_umap", group.by = "orig.ident", pt.size = 0.3, raster = FALSE, cols = palette) +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
    )
dev.off()

#######################################################################################################################################################################################
#######################################################################################################################################################################################

### fastMNN
if (Primera_vez) {
    Seurat_fastMNN <- Seurat_merged

    Seurat_fastMNN <- FindVariableFeatures(Seurat_fastMNN,
        selection.method = "vst",
        nfeatures = 2000
    )

    Seurat_fastMNN <- ScaleData(Seurat_fastMNN)
    Seurat_fastMNN$Marked_cells <- NULL
    Seurat_fastMNN2 <- DietSeurat(Seurat_fastMNN)
    Seurat_fastMNN2 <- RunFastMNN(SplitObject(Seurat_fastMNN2, split.by = "Product_norm"))
    Seurat_fastMNN2 <- RunUMAP(Seurat_fastMNN2, reduction = "mnn", dims = 1:30)

    Seurat_fastMNN2 <- FindNeighbors(Seurat_fastMNN2, reduction = "mnn", dims = 1:30)
    Seurat_fastMNN2 <- FindClusters(Seurat_fastMNN2)

    saveRDS(Seurat_fastMNN2, "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/fastMNN/Sin_GT/Seurat_fastMNN2.RDS")

    # Set path
    path_base <- "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/fastMNN/Sin_GT/Seurat_fastMNN2"

    # Save .h5Seurat
    SaveH5Seurat(Seurat_fastMNN2, filename = paste0(path_base, ".h5Seurat"))

    # Transform into .h5ad
    Convert(paste0(path_base, ".h5Seurat"), dest = "h5ad", overwrite = TRUE)
} else {
    Seurat_fastMNN2 <- readRDS("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration_methods_lab/fastMNN/Sin_GT/Seurat_fastMNN2.RDS")
}

pdf("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Suplementarias/Suplementaria_S2_7_fastMNN2.pdf")
DimPlot(Seurat_fastMNN2, reduction = "umap", group.by = "orig.ident", pt.size = 0.3, raster = FALSE, cols = palette) +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
    )
dev.off()

#######################################################################################################################################################################################
#######################################################################################################################################################################################

##### End of script #####
