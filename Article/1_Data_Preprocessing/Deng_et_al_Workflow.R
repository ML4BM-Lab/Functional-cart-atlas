###############################################################################
###############################################################################

# Program: Deng_et_al_Workflow.R
# Author: Sergio Cámara Peña
# Date: 2023
# Version: FINAL

###############################################################################
###############################################################################


##### Load required libraries #####
library(DropletUtils)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(DoubletFinder)
library(cowplot)
library(scales)
library(patchwork)


##### Set seed #####
set.seed(2504)


###### Data loading from CellRanger ######
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Count_Matrices/Cell_Ranger/")
Patient_list <- gsub("_cellranger_count.out", "", dir()[dir() %>% grep(pattern = ".out")])

Seurat_list_Deng <- list()
Sum_CAR <- list()
Table_CAR <- list()
for (contador in Patient_list) {
    cell_ranger_dir <- paste0("./run_", contador, "_Deng_et_al/outs/filtered_feature_bc_matrix")
    cell_ranger_data <- Read10X(data.dir = cell_ranger_dir)
    Seurat_list_Deng[[contador]] <- CreateSeuratObject(counts = cell_ranger_data, project = "Deng_et_al", min.features = 200, min.cells = 3)
    Seurat_list_Deng[[contador]]$Product <- contador
    Sum_CAR[[contador]] <- sum(Seurat_list_Deng[[contador]]$RNA@counts[which(rownames(Seurat_list_Deng[[contador]]) == "FMC63-28Z"), ] > 0)
    Table_CAR[[contador]] <- table(Seurat_list_Deng[[contador]]$RNA@counts[which(rownames(Seurat_list_Deng[[contador]]) == "FMC63-28Z"), Seurat_list_Deng[[contador]]$RNA@counts[which(rownames(Seurat_list_Deng[[contador]]) == "FMC63-28Z"), ] > 0])

    Seurat_list_Deng[[contador]] <- RenameCells(Seurat_list_Deng[[contador]], new.names = gsub("-1", "", colnames(Seurat_list_Deng[[contador]])))

    print(Seurat_list_Deng[[contador]])
    print(Sum_CAR[[contador]])
    print(paste0(contador, " DONE"))
}

rm(contador, cell_ranger_dir, cell_ranger_data)


##### Quality control #####

### Add more metadata: novelty score and mito ratio
for (contador in seq_along(Seurat_list_Deng)) {
    Seurat_list_Deng[[contador]]$log10GenesPerUMI <- log10(Seurat_list_Deng[[contador]]$nFeature_RNA) / log10(Seurat_list_Deng[[contador]]$nCount_RNA)

    Seurat_list_Deng[[contador]]$mitoRatio <- (PercentageFeatureSet(object = Seurat_list_Deng[[contador]], pattern = "^MT-")) / 100

    Seurat_list_Deng[[contador]]$log10GenesPerUMI <- log10(Seurat_list_Deng[[contador]]$nFeature_RNA) / log10(Seurat_list_Deng[[contador]]$nCount_RNA)
}
rm(contador)

### Pre-filter
minCov <- 1000 # if a sample has a good coverage, then it doesn't set a lower thresold for nCount.
countLOW <- list()
countHIGH <- list()
featureLOW <- list()
featureHIGH <- list()

if (exists("Pre_filter")) {
    print("Prefiltering DONE")
} else {
    for (contador in seq_along(Seurat_list_Deng)) {
        if (min(Seurat_list_Deng[[contador]]$nCount_RNA) >= minCov) {
            countLOW[[contador]] <- min(Seurat_list_Deng[[contador]]$nCount_RNA)
        } else {
            countLOW[[contador]] <- quantile(Seurat_list_Deng[[contador]]$nCount_RNA, prob = 0.01)
        }
        countHIGH[[contador]] <- quantile(Seurat_list_Deng[[contador]]$nCount_RNA, prob = 0.99)
        featureLOW[[contador]] <- quantile(Seurat_list_Deng[[contador]]$nFeature_RNA, prob = 0.01)
        featureHIGH[[contador]] <- quantile(Seurat_list_Deng[[contador]]$nFeature_RNA, prob = 0.99)

        ## subset
        Seurat_list_Deng[[contador]] <- subset(Seurat_list_Deng[[contador]], subset = (nFeature_RNA > featureLOW[[contador]]) &
            (nFeature_RNA < featureHIGH[[contador]]) & (nCount_RNA > countLOW[[contador]]) & (nCount_RNA < countHIGH[[contador]]) & (mitoRatio < 0.20))
    }
    Pre_filter <- "DONE"
}
rm(minCov)

### Group all Products in a single metadata table
merged_metadata_df <- list()
for (cont in seq_along(Seurat_list_Deng)) {
    merged_metadata_df[[cont]] <- Seurat_list_Deng[[cont]]@meta.data
}
rm(cont)
merged_metadata_df <- bind_rows(merged_metadata_df)

### Exploratory analysis
# Values
# Order: Patient_14, Patient_15, Patient_16, Patient_18, Patient_20, Patient_21, Patient_26, Patient_27, Patient_28, Patient_33, Patient_34, Patient_37,
#        Patient_38, Patient_40, Patient_41, Patient_42, Patient_43, Patient_49, Patient_50, Patient_54, Patient_55, Patient_56, Patient_59, Patient_64

Min_feat <- c(
    500, 750, 1100, 700, 1000, 1300, 1000, 800, 700, 750, 800, 1400,
    800, 1500, 1500, 900, 650, 1600, 1500, 1200, 1500, 1000, 1750, 1500
)
Min_counts <- c(
    900, 1400, 2200, 1500, 3000, 3000, 2000, 1500, 1400, 1250, 1500, 3100,
    1250, 3000, 3000, 1500, 1000, 4000, 4000, 3000, 4000, 2200, 4500, 3000
)
Max_counts <- c(
    26000, 12500, 20000, 15000, 22500, 40000, 30000, 21000, 12000, 27500, 12500, 15500,
    22000, 35000, 40000, 22000, 15000, 42000, 42000, 38000, 45000, 22000, 35000, 44000
)
Max_mito_ratio <- c(
    0.1, 0.1, 0.08, 0.09, 0.2, 0.11, 0.09, 0.09, 0.1, 0.2, 0.1, 0.2,
    0.08, 0.1, 0.1, 0.2, 0.2, 0.12, 0.12, 0.1, 0.12, 0.13, 0.15, 0.13
)

for (i in seq_along(Seurat_list_Deng)) {
    a <- max(Seurat_list_Deng[[i]]@meta.data$nCount_RNA)
    print(paste0(i, ": ", a))
}

colores <- hue_pal()(length(Seurat_list_Deng))

# Visualize the number of cell counts per sample
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Exploratory_analysis/Pre_QC_Ncells.pdf")
merged_metadata_df %>%
    ggplot(aes(x = Product, fill = Product)) +
    geom_bar() +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 3) +
    theme_classic() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    ggtitle("Number of Cells pre-QC")
dev.off()

# Transformation into single cell experiment + Knee plots representation
counter <- 0
sce_Seurat_list_Deng <- list()
metadata_bcrank <- list()
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Exploratory_analysis/Pre_QC_Knee_plots.pdf")
for (i in Patient_list) {
    counter <- counter + 1
    sce_Seurat_list_Deng[[i]] <- as.SingleCellExperiment(Seurat_list_Deng[[i]])
    bcrank <- DropletUtils::barcodeRanks(counts(sce_Seurat_list_Deng[[i]]))

    uniq <- !duplicated(bcrank$rank) # Only showing unique points for plotting speed.

    plot(bcrank$rank[uniq], bcrank$total[uniq], log = "xy", xlab = "Rank", ylab = "Total UMI count", main = paste0(i), cex.lab = 1.2)
    abline(h = metadata(bcrank)$inflection, col = "darkgreen", lty = 2)
    abline(h = metadata(bcrank)$knee, col = "dodgerblue", lty = 2)
    abline(h = Max_counts[counter], col = "red", lty = 2)
    abline(h = Min_counts[counter], col = "darkorchid1", lty = 2)
    legend("bottomleft", legend = c("Inflection", "Knee", "Doublets_estim", "My_thresh"), col = c("darkgreen", "dodgerblue", "red", "darkorchid1"), lty = 2, cex = 1.2)

    metadata_bcrank[[i]] <- metadata(bcrank)
    print(i)
}
dev.off()
rm(i, bcrank, uniq, counter)

# Visualize the number UMIs/transcripts per cell
pt_h <- list()
for (i in seq_along(Seurat_list_Deng)) {
    pt_h[[i]] <- ggplot(data = Seurat_list_Deng[[i]]@meta.data, aes(x = nCount_RNA)) +
        geom_density(alpha = 0.2, fill = colores[i]) +
        scale_x_log10() +
        theme_classic() +
        ylab("Cell density") +

        # Adjust thresholds per dataset
        geom_vline(xintercept = Min_counts[i], col = colores[i], alpha = 0.8, linetype = "dotted") +
        geom_vline(xintercept = Max_counts[i], col = colores[i], alpha = 0.8, linetype = "dotted") +
        labs(
            title = "Number UMIs/transcripts per cell pre-QC:",
            subtitle = paste0(unique(Seurat_list_Deng[[i]]@meta.data$Product))
        )
}

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Exploratory_analysis/Pre_QC_UMIs_per_cell.pdf")
pt_h
dev.off()
rm(i, pt_h)

# Visualize the distribution of genes detected per cell
pt_h <- list()
for (i in seq_along(Seurat_list_Deng)) {
    pt_h[[i]] <- ggplot(data = Seurat_list_Deng[[i]]@meta.data, aes(x = nFeature_RNA)) +
        geom_density(alpha = 0.2, fill = colores[i]) +
        theme_classic() +
        scale_x_log10() +

        # Adjust thresholds per dataset
        geom_vline(xintercept = Min_feat[i], col = colores[i], alpha = 0.6, linetype = "dotted") +
        labs(
            title = "Number genes per cell pre-QC:",
            subtitle = paste0(unique(Seurat_list_Deng[[i]]@meta.data$Product))
        )
}

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Exploratory_analysis/Pre_QC_Genes_per_cell.pdf")
pt_h
dev.off()
rm(i, pt_h)

# Visualize the distribution of mitochondrial gene expression detected per cell
pt_h <- list()
for (i in seq_along(Seurat_list_Deng)) {
    pt_h[[i]] <- ggplot(data = Seurat_list_Deng[[i]]@meta.data, aes(x = mitoRatio)) +
        geom_density(alpha = 0.2, fill = colores[i]) +
        scale_x_log10() +
        theme_classic() +

        # Adjust thresholds per dataset
        geom_vline(xintercept = Max_mito_ratio[i], col = colores[i], alpha = 0.6, linetype = "dotted") +
        labs(
            title = "Mitocondrial ratio per cell:",
            subtitle = paste0(unique(Seurat_list_Deng[[i]]@meta.data$Product))
        )
}

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Exploratory_analysis/Pre_QC_Mito_ratio_per_cell.pdf")
pt_h
dev.off()
rm(i, pt_h)

# Visualize the distribution of complexity per cell
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Exploratory_analysis/Pre_QC_Complexity_per_cell.pdf")
merged_metadata_df %>%
    ggplot(aes(color = Product, x = log10GenesPerUMI, fill = Product)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    ggtitle("Complexity per cell")
dev.off()

# Violin plots
p_Vln <- list()
for (cont in seq_along(Seurat_list_Deng)) {
    p_Vln[[cont]] <- VlnPlot(Seurat_list_Deng[[cont]], features = c("nCount_RNA", "nFeature_RNA", "mitoRatio"), col = colores[cont], pt.size = 0, group.by = "Product")
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Exploratory_analysis/Pre_QC_Vln_plots.pdf")
p_Vln
dev.off()

rm(cont, p_Vln)

# nGenes vs nUMI graph (Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs)
p1 <- list()
for (i in seq_along(Seurat_list_Deng)) {
    p1[[i]] <- ggplot(data = Seurat_list_Deng[[i]]@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = mitoRatio)) +
        geom_point(size = 0.7) +
        scale_color_gradientn(colors = c("black", "blue", "green2", "red", "yellow")) +
        stat_smooth(method = lm) +
        scale_x_log10() +
        scale_y_log10() +
        theme_classic() +
        ggtitle(paste0(Patient_list[i])) +
        geom_vline(xintercept = Min_counts[i], col = "black") +
        geom_vline(xintercept = Max_counts[i], col = "black") +
        geom_hline(yintercept = Min_feat[i], col = "black")
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Exploratory_analysis/Pre_QC_nGenes_vs_nUMI.pdf")
p1
dev.off()

rm(i, p1)

# mitoRatio vs nUMI graph
p1 <- list()
for (i in seq_along(Seurat_list_Deng)) {
    p1[[i]] <- ggplot(data = Seurat_list_Deng[[i]]@meta.data, aes(x = nCount_RNA, y = mitoRatio)) +
        geom_point(size = 0.7) +
        theme_classic() +
        ggtitle(paste0(Patient_list[i]))
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Exploratory_analysis/Pre_QC_mitoRatio_vs_nUMI.pdf")
p1
dev.off()

rm(i, p1)

### Subset samples by determined quality control thresholds
for (i in seq_along(Seurat_list_Deng)) {
    Seurat_list_Deng[[i]] <- subset(Seurat_list_Deng[[i]],
        subset = (nFeature_RNA > Min_feat[i]) &
            (nCount_RNA > Min_counts[i]) &
            (nCount_RNA < Max_counts[i]) &
            (mitoRatio < Max_mito_ratio[i]) &
            (log10GenesPerUMI > 0.8)
    )
}
rm(i)

### Calculate Doublets
# https://github.com/chris-mcginnis-ucsf/DoubletFinder
# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Predict_doublets

Seurat_doublets_finder <- Seurat_list_Deng

sweep_res_list_Deng <- list()
sweep_stats_list_Deng <- list()
bcmvn_list_Deng <- list()
pK <- list()
BCmetric <- list()
pK_choose <- list()
nExp <- list()
per_doubl <- c(
    0.104, 0.096, 0.052, 0.036, 0.024, 0.064, 0.056, 0.040, 0.052, 0.096, 0.072, 0.044,
    0.112, 0.088, 0.080, 0.040, 0.136, 0.080, 0.080, 0.080, 0.056, 0.076, 0.068, 0.076
) # Calculated looking at 10X table

# pK Identification (no ground-truth)
for (i in seq_along(Seurat_doublets_finder)) {
    Seurat_doublets_finder[[i]] <- NormalizeData(Seurat_doublets_finder[[i]])
    Seurat_doublets_finder[[i]] <- FindVariableFeatures(Seurat_doublets_finder[[i]], selection.method = "vst", nfeatures = 2000)
    Seurat_doublets_finder[[i]] <- ScaleData(Seurat_doublets_finder[[i]])
    Seurat_doublets_finder[[i]] <- RunPCA(Seurat_doublets_finder[[i]], npcs = 20)
    Seurat_doublets_finder[[i]] <- RunUMAP(Seurat_doublets_finder[[i]], dims = 1:10)
    nExp[[i]] <- round(ncol(Seurat_doublets_finder[[i]]) * per_doubl[i])

    sweep_res_list_Deng[[i]] <- paramSweep_v3(Seurat_doublets_finder[[i]], PCs = 1:10, sct = FALSE, num.cores = 12)
    sweep_stats_list_Deng[[i]] <- summarizeSweep(sweep_res_list_Deng[[i]], GT = FALSE)
    bcmvn_list_Deng[[i]] <- find.pK(sweep_stats_list_Deng[[i]])
    dev.off()

    pK[[i]] <- as.numeric(as.character(bcmvn_list_Deng[[i]]$pK))
    BCmetric[[i]] <- bcmvn_list_Deng[[i]]$BCmetric
    pK_choose[[i]] <- as.numeric(pK[[i]][which(BCmetric[[i]] %in% max(BCmetric[[i]]))])
    print(pK_choose[[i]])
}

rm(i)

# Plot data
h1 <- list()
for (i in seq_along(Seurat_doublets_finder)) {
    df <- data.frame(x = pK[[i]], y = BCmetric[[i]])
    h1[[i]] <- ggplot(df, aes(x = x, y = y)) +
        geom_point(shape = 16, size = 3, color = "blue") +
        geom_line(color = "blue", lty = 1) +
        geom_vline(xintercept = pK_choose[[i]], color = "red", linetype = 2, size = 1) +
        geom_text(aes_(x = pK_choose[[i]], y = max(BCmetric[[i]]), label = as.character(pK_choose[[i]])),
            hjust = 1.2, vjust = -0.2, color = "red", size = 6
        ) +

        # Set theme and labels
        theme_minimal() +
        labs(title = "The BCmvn distributions", x = "pK", y = "BCmetric") +

        # Adjust axis limits and text size
        scale_x_continuous(limits = c(min(as.numeric(pK[[i]])), max(as.numeric(pK[[i]])))) +
        theme(
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 16, face = "bold")
        )
}

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Exploratory_analysis/Doublet_Finder_pK_calculation.pdf")
h1
dev.off()
rm(i, h1)

# Doublets prediction
h2 <- list()
DF_name <- list()
for (i in seq_along(Seurat_doublets_finder)) {
    Seurat_doublets_finder[[i]] <- doubletFinder_v3(Seurat_doublets_finder[[i]], pN = 0.25, pK = pK_choose[[i]], nExp = nExp[[i]], PCs = 1:10)
    # name of the DF prediction can change, so extract the correct column name.
    DF_name[[i]] <- colnames(Seurat_doublets_finder[[i]]@meta.data)[grepl("DF.classification", colnames(Seurat_doublets_finder[[i]]@meta.data))]
    h2[[i]] <- cowplot::plot_grid(
        nrow = 2, DimPlot(Seurat_doublets_finder[[i]], group.by = DF_name[[i]]) + NoAxes() + ggtitle(paste0(Patient_list[i])),
        VlnPlot(Seurat_doublets_finder[[i]], features = "nFeature_RNA", group.by = DF_name[[i]], pt.size = 0.1)
    )
}

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Exploratory_analysis/Doublet_Finder_Classification.pdf")
h2
dev.off()
rm(i, h2)

# Remove all predicted doublets from our data.
for (i in seq_along(Seurat_doublets_finder)) {
    print(Seurat_doublets_finder[[i]])

    Sys.sleep(2)

    Singlets <- rownames(Seurat_doublets_finder[[i]]@meta.data)[Seurat_doublets_finder[[i]]@meta.data[, DF_name[[i]]] == "Singlet"]
    print(length(Singlets))

    Sys.sleep(2)

    Seurat_list_Deng[[i]] <- subset(Seurat_list_Deng[[i]], cells = Singlets)
    print(Seurat_list_Deng[[i]])

    Sys.sleep(2)

    print(paste0(Patient_list[i], " DONE"))
    cat(paste0("----------------------------------------------------------------------------\n"))
}

### Subset by CAR+ cells
for(i in seq_along(Seurat_list_Deng)){
    Seurat_list_Deng[[i]] <- subset(x = Seurat_list_Deng[[i]], subset = `FMC63-28Z` > 0)
}

Seurat_list_Deng

### Re-doing of exploratory analysis graphs after QC
# Group again all Products in a single metadata table
merged_metadata_df <- list()
for (cont in seq_along(Seurat_list_Deng)) {
    merged_metadata_df[[cont]] <- Seurat_list_Deng[[cont]]@meta.data
}
rm(cont)
merged_metadata_df <- bind_rows(merged_metadata_df)

# Visualize the number of cell counts per sample
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Exploratory_analysis/Post_QC_Ncells.pdf")
merged_metadata_df %>%
    ggplot(aes(x = Product, fill = Product)) +
    geom_bar() +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 4) +
    theme_classic() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    ggtitle("Number of Cells post-QC")
dev.off()

# Visualize the number UMIs/transcripts per cell
pt_h <- list()
for (i in seq_along(Seurat_list_Deng)) {
    pt_h[[i]] <- ggplot(data = Seurat_list_Deng[[i]]@meta.data, aes(x = nCount_RNA)) +
        geom_density(alpha = 0.2, fill = colores[i]) +
        scale_x_log10() +
        theme_classic() +
        ylab("Cell density") +
        labs(
            title = "Number UMIs/transcripts per cell post-QC:",
            subtitle = paste0(unique(Seurat_list_Deng[[i]]@meta.data$Product))
        )
}

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Exploratory_analysis/Post_QC_UMIs_per_cell.pdf")
pt_h
dev.off()
rm(i, pt_h)

# Visualize the distribution of genes detected per cell
pt_h <- list()
for (i in seq_along(Seurat_list_Deng)) {
    pt_h[[i]] <- ggplot(data = Seurat_list_Deng[[i]]@meta.data, aes(x = nFeature_RNA)) +
        geom_density(alpha = 0.2, fill = colores[i]) +
        theme_classic() +
        scale_x_log10() +
        labs(
            title = "Number genes per cell post-QC:",
            subtitle = paste0(unique(Seurat_list_Deng[[i]]@meta.data$Product))
        )
}

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Exploratory_analysis/Post_QC_Genes_per_cell.pdf")
pt_h
dev.off()
rm(i, pt_h)

# Visualize the distribution of mitochondrial gene expression detected per cell
pt_h <- list()
for (i in seq_along(Seurat_list_Deng)) {
    pt_h[[i]] <- ggplot(data = Seurat_list_Deng[[i]]@meta.data, aes(x = mitoRatio)) +
        geom_density(alpha = 0.2, fill = colores[i]) +
        scale_x_log10() +
        theme_classic() +

        # Adjust thresholds per dataset
        labs(
            title = "Mitocondrial ratio per cell:",
            subtitle = paste0(unique(Seurat_list_Deng[[i]]@meta.data$Product))
        )
}

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Exploratory_analysis/Post_QC_Mito_ratio_per_cell.pdf")
pt_h
dev.off()
rm(i, pt_h)

# Violin plots
p_Vln <- list()
for (cont in seq_along(Seurat_list_Deng)) {
    p_Vln[[cont]] <- VlnPlot(Seurat_list_Deng[[cont]], features = c("nCount_RNA", "nFeature_RNA", "mitoRatio"), col = colores[cont], pt.size = 0, group.by = "Product")
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Exploratory_analysis/Post_QC_Vln_plots.pdf")
p_Vln
dev.off()

rm(cont, p_Vln)

# nGenes vs nUMI graph (Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs)
p1 <- list()
for (i in seq_along(Seurat_list_Deng)) {
    p1[[i]] <- ggplot(data = Seurat_list_Deng[[i]]@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = mitoRatio, grou)) +
        geom_point(size = 0.7) +
        scale_color_gradientn(colors = c("black", "blue", "green2", "red", "yellow")) +
        stat_smooth(method = lm) +
        scale_x_log10() +
        scale_y_log10() +
        theme_classic() +
        ggtitle(paste0(Patient_list[i]))
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Exploratory_analysis/Post_QC_nGenes_vs_nUMI.pdf")
p1
dev.off()

rm(i, p1)

setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/RDS")
saveRDS(Seurat_list_Deng, file = "PostQC_CellRanger_Deng_RDS.rds")

##### Normalization #####

# rm(list = ls())
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/RDS")
Seurat_list_Deng <- readRDS("PostQC_CellRanger_Deng_RDS.rds")
load("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Codigo/Gene_Markers_Info/cycle.rda")

for (i in seq_along(Seurat_list_Deng)) {
    Seurat_list_Deng[[i]] <- NormalizeData(Seurat_list_Deng[[i]])

    # Score cells for cell cycle
    Seurat_list_Deng[[i]] <- CellCycleScoring(Seurat_list_Deng[[i]],
        g2m.features = g2m_genes,
        s.features = s_genes
    )

    Seurat_list_Deng[[i]] <- FindVariableFeatures(Seurat_list_Deng[[i]],
        selection.method = "vst",
        nfeatures = 2000
    )

    Seurat_list_Deng[[i]] <- ScaleData(Seurat_list_Deng[[i]])
}
rm(i)

saveRDS(Seurat_list_Deng, file = "Normalized_CellRanger_Deng_RDS.rds")


##### Merge ---- Worst case scenario #####
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/RDS")
Seurat_list_Deng <- readRDS("Normalized_CellRanger_Deng_RDS.rds")

Seurat_list_Deng_2 <- merge(x = Seurat_list_Deng[[1]], y = c(
    Seurat_list_Deng[[2]], Seurat_list_Deng[[3]],
    Seurat_list_Deng[[4]], Seurat_list_Deng[[5]], Seurat_list_Deng[[6]],
    Seurat_list_Deng[[7]], Seurat_list_Deng[[8]], Seurat_list_Deng[[9]],
    Seurat_list_Deng[[10]], Seurat_list_Deng[[11]], Seurat_list_Deng[[12]],
    Seurat_list_Deng[[13]], Seurat_list_Deng[[14]], Seurat_list_Deng[[15]],
    Seurat_list_Deng[[16]], Seurat_list_Deng[[17]], Seurat_list_Deng[[18]],
    Seurat_list_Deng[[19]], Seurat_list_Deng[[20]], Seurat_list_Deng[[21]],
    Seurat_list_Deng[[22]], Seurat_list_Deng[[23]], Seurat_list_Deng[[24]]
))

Seurat_list_Deng <- Seurat_list_Deng_2
rm(Seurat_list_Deng_2)

Seurat_list_Deng <- FindVariableFeatures(Seurat_list_Deng,
    selection.method = "vst",
    nfeatures = 2000
)

Seurat_list_Deng <- ScaleData(Seurat_list_Deng)

# Run PCA
Seurat_list_Deng <- RunPCA(object = Seurat_list_Deng, reduction.name = "pca_wo_integ")

# ElbowPlot(Seurat_list_Deng, ndims = 50, reduction = "pca_wo_integ")
# dev.off()

# Run UMAP
Seurat_list_Deng <- RunUMAP(Seurat_list_Deng,
    dims = 1:30,
    reduction = "pca_wo_integ",
    reduction.name = "umap_wo_integ"
)

# Plot UMAP
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Worst_case_scenario/WO_integ_Seurat.pdf")
DimPlot(Seurat_list_Deng,
    group.by = "Product",
    reduction = "umap_wo_integ"
)
dev.off()

# Quality control metrics
p <- FeaturePlot(Seurat_list_Deng,
    reduction = "umap_wo_integ",
    features = c("nCount_RNA", "nFeature_RNA", "log10GenesPerUMI", "mitoRatio"),
    order = TRUE
)

# Seurat_list_Deng$Marker <- 0
# for (i in seq_along(rownames(Seurat_list_Deng@meta.data))) {
#     if (Seurat_list_Deng$Product[i] == "Patient_43") {
#         Seurat_list_Deng$Marker[i] <- "SI"
#     } else {
#         Seurat_list_Deng$Marker[i] <- "NO"
#     }
# }
#
# which(Seurat_list_Deng$Marker == "SI")
#
# DimPlot(Seurat_list_Deng, reduction = "umap_wo_integ", group.by = "Marker", order = TRUE)

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Deng_et_al/Plots/Worst_case_scenario/QC_metrics_WO_integ_Seurat.pdf")
p + plot_annotation(title = paste0(Seurat_list_Deng@meta.data$orig.ident %>% unique()), theme = theme(plot.title = element_text(size = 16)))
dev.off()

rm(p)

##### END OF SCRIPT #####