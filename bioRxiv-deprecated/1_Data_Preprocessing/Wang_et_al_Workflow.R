###############################################################################
###############################################################################

# Program: Wang_et_al_Workflow.R
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


###### Load data from CellRanger ######
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Count_Matrices/Cell_Ranger")
Patient_list <- gsub("_Wang_et_al", "", gsub("run_", "", dir()[dir() %>% grep(pattern = "^run_")]))
Patient_list

Seurat_list_Wang <- list()
Sum_CAR <- list()
Table_CAR <- list()
for (contador in Patient_list) {
    cell_ranger_dir <- paste0("./run_", contador, "_Wang_et_al/outs/filtered_feature_bc_matrix")
    cell_ranger_data <- Read10X(data.dir = cell_ranger_dir)
    Seurat_list_Wang[[contador]] <- CreateSeuratObject(counts = cell_ranger_data, project = "Wang_et_al", min.features = 200, min.cells = 3)
    Seurat_list_Wang[[contador]]$Product <- contador
    Seurat_list_Wang[[contador]] <- RenameCells(Seurat_list_Wang[[contador]], new.names = gsub("-1", "", colnames(Seurat_list_Wang[[contador]])))
    Sum_CAR[[contador]] <- sum(Seurat_list_Wang[[contador]]$RNA@counts[which(rownames(Seurat_list_Wang[[contador]]) == "IFNB1-5LTR"), ] > 0)
    Table_CAR[[contador]] <- table(Seurat_list_Wang[[contador]]$RNA@counts[which(rownames(Seurat_list_Wang[[contador]]) == "IFNB1-5LTR"), Seurat_list_Wang[[contador]]$RNA@counts[which(rownames(Seurat_list_Wang[[contador]]) == "IFNB1-5LTR"), ] > 0])

    print(Seurat_list_Wang[[contador]])
    print(Sum_CAR[[contador]])
    print(paste0(contador, " DONE"))
}

rm(contador, cell_ranger_dir, cell_ranger_data)


##### Quality control #####

### Add more metadata: novelty score and mito ratio
for (contador in seq_along(Seurat_list_Wang)) {
    Seurat_list_Wang[[contador]]$log10GenesPerUMI <- log10(Seurat_list_Wang[[contador]]$nFeature_RNA) / log10(Seurat_list_Wang[[contador]]$nCount_RNA)

    Seurat_list_Wang[[contador]]$mitoRatio <- (PercentageFeatureSet(object = Seurat_list_Wang[[contador]], pattern = "^MT-")) / 100

    Seurat_list_Wang[[contador]]$log10GenesPerUMI <- log10(Seurat_list_Wang[[contador]]$nFeature_RNA) / log10(Seurat_list_Wang[[contador]]$nCount_RNA)
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
    for (contador in seq_along(Seurat_list_Wang)) {
        if (min(Seurat_list_Wang[[contador]]$nCount_RNA) >= minCov) {
            countLOW[[contador]] <- min(Seurat_list_Wang[[contador]]$nCount_RNA)
        } else {
            countLOW[[contador]] <- quantile(Seurat_list_Wang[[contador]]$nCount_RNA, prob = 0.01)
        }
        countHIGH[[contador]] <- quantile(Seurat_list_Wang[[contador]]$nCount_RNA, prob = 0.99)
        featureLOW[[contador]] <- quantile(Seurat_list_Wang[[contador]]$nFeature_RNA, prob = 0.01)
        featureHIGH[[contador]] <- quantile(Seurat_list_Wang[[contador]]$nFeature_RNA, prob = 0.99)

        ## subset
        Seurat_list_Wang[[contador]] <- subset(Seurat_list_Wang[[contador]], subset = (nFeature_RNA > featureLOW[[contador]]) &
            (nFeature_RNA < featureHIGH[[contador]]) & (nCount_RNA > countLOW[[contador]]) & (nCount_RNA < countHIGH[[contador]]) & (mitoRatio < 0.20))
    }
    Pre_filter <- "DONE"
}
rm(minCov)

### Group all Products in a single metadata table
merged_metadata_df <- list()
for (cont in seq_along(Seurat_list_Wang)) {
    merged_metadata_df[[cont]] <- Seurat_list_Wang[[cont]]@meta.data
}
rm(cont)
merged_metadata_df <- bind_rows(merged_metadata_df)

### Exploratory analysis
# Values
# Order: Product_Donor_1, Product_Donor_2, Product_Donor_3, Stimulated_Product_Donor_1, Stimulated_Product_Donor_2, Stimulated_Product_Donor_3

Min_feat <- c(1400, 1400, 1200, 925, 1050, 1000)
Min_counts <- c(3500, 3500, 3300, 2800, 2200, 2800)
Max_counts <- c(17000, 20000, 18000, 25000, 27000, 30000)
Max_mito_ratio <- c(0.09, 0.09, 0.09, 0.09, 0.09, 0.09)

for (i in seq_along(Seurat_list_Wang)) {
    a <- max(Seurat_list_Wang[[i]]@meta.data$nCount_RNA)
    print(paste0(i, ": ", a))
}

colores <- hue_pal()(length(Seurat_list_Wang))

# Visualize the number of cell counts per sample
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Exploratory_analysis/Pre_QC_Ncells.pdf")
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
    ggtitle("Number of Cells pre-QC")
dev.off()

# Transformation into single cell experiment + Knee plots representation
counter <- 0
sce_Seurat_list_Wang <- list()
metadata_bcrank <- list()
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Exploratory_analysis/Pre_QC_Knee_plots.pdf")
for (i in Patient_list) {
    counter <- counter + 1
    sce_Seurat_list_Wang[[i]] <- as.SingleCellExperiment(Seurat_list_Wang[[i]])
    bcrank <- DropletUtils::barcodeRanks(counts(sce_Seurat_list_Wang[[i]]))

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
for (i in seq_along(Seurat_list_Wang)) {
    pt_h[[i]] <- ggplot(data = Seurat_list_Wang[[i]]@meta.data, aes(x = nCount_RNA)) +
        geom_density(alpha = 0.2, fill = colores[i]) +
        scale_x_log10() +
        theme_classic() +
        ylab("Cell density") +

        # Adjust thresholds per dataset
        geom_vline(xintercept = Min_counts[i], col = colores[i], alpha = 0.8, linetype = "dotted") +
        geom_vline(xintercept = Max_counts[i], col = colores[i], alpha = 0.8, linetype = "dotted") +
        labs(
            title = "Number UMIs/transcripts per cell pre-QC:",
            subtitle = paste0(unique(Seurat_list_Wang[[i]]@meta.data$Product))
        )
}

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Exploratory_analysis/Pre_QC_UMIs_per_cell.pdf")
pt_h
dev.off()
rm(i, pt_h)

# Visualize the distribution of genes detected per cell
pt_h <- list()
for (i in seq_along(Seurat_list_Wang)) {
    pt_h[[i]] <- ggplot(data = Seurat_list_Wang[[i]]@meta.data, aes(x = nFeature_RNA)) +
        geom_density(alpha = 0.2, fill = colores[i]) +
        theme_classic() +
        scale_x_log10() +

        # Adjust thresholds per dataset
        geom_vline(xintercept = Min_feat[i], col = colores[i], alpha = 0.6, linetype = "dotted") +
        labs(
            title = "Number genes per cell pre-QC:",
            subtitle = paste0(unique(Seurat_list_Wang[[i]]@meta.data$Product))
        )
}

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Exploratory_analysis/Pre_QC_Genes_per_cell.pdf")
pt_h
dev.off()
rm(i, pt_h)

# Visualize the distribution of mitochondrial gene expression detected per cell
pt_h <- list()
for (i in seq_along(Seurat_list_Wang)) {
    pt_h[[i]] <- ggplot(data = Seurat_list_Wang[[i]]@meta.data, aes(x = mitoRatio)) +
        geom_density(alpha = 0.2, fill = colores[i]) +
        scale_x_log10() +
        theme_classic() +

        # Adjust thresholds per dataset
        geom_vline(xintercept = Max_mito_ratio[i], col = colores[i], alpha = 0.6, linetype = "dotted") +
        labs(
            title = "Mitocondrial ratio per cell:",
            subtitle = paste0(unique(Seurat_list_Wang[[i]]@meta.data$Product))
        )
}

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Exploratory_analysis/Pre_QC_Mito_ratio_per_cell.pdf")
pt_h
dev.off()
rm(i, pt_h)

# Visualize the distribution of complexity per cell
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Exploratory_analysis/Pre_QC_Complexity_per_cell.pdf")
merged_metadata_df %>%
    ggplot(aes(color = Product, x = log10GenesPerUMI, fill = Product)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    ggtitle("Complexity per cell")
dev.off()

# Violin plots
p_Vln <- list()
for (cont in seq_along(Seurat_list_Wang)) {
    p_Vln[[cont]] <- VlnPlot(Seurat_list_Wang[[cont]], features = c("nCount_RNA", "nFeature_RNA", "mitoRatio"), col = colores[cont], group.by = "Product", pt.size = 0)
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Exploratory_analysis/Pre_QC_Vln_plots.pdf")
p_Vln
dev.off()

rm(cont, p_Vln)

# nGenes vs nUMI graph (Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs)
p1 <- list()
for (i in seq_along(Seurat_list_Wang)) {
    p1[[i]] <- ggplot(data = Seurat_list_Wang[[i]]@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = mitoRatio)) +
        geom_point(size = 0.7) +
        scale_color_gradientn(colors = c("black", "blue", "green2", "red", "yellow")) +
        stat_smooth(method = lm) +
        scale_x_log10() +
        scale_y_log10() +
        theme_classic() +
        ggtitle(paste0(Patient_list[i])) +
        geom_vline(xintercept = Min_counts[i], col = "black") +
        geom_hline(yintercept = Min_feat[i], col = "black")
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Exploratory_analysis/Pre_QC_nGenes_vs_nUMI.pdf")
p1
dev.off()

rm(i, p1)

# mitoRatio vs nUMI graph
p1 <- list()
for (i in seq_along(Seurat_list_Wang)) {
    p1[[i]] <- ggplot(data = Seurat_list_Wang[[i]]@meta.data, aes(x = nCount_RNA, y = mitoRatio)) +
        geom_point(size = 0.7) +
        theme_classic() +
        ggtitle(paste0(Patient_list[i]))
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Exploratory_analysis/Pre_QC_mitoRatio_vs_nUMI.pdf")
p1
dev.off()

rm(i, p1)

### Subset samples by determined quality control thresholds
for (i in seq_along(Seurat_list_Wang)) {
    Seurat_list_Wang[[i]] <- subset(Seurat_list_Wang[[i]],
        subset = (nFeature_RNA > Min_feat[i]) &
            (nCount_RNA > Min_counts[i]) &
            (nCount_RNA < Max_counts[i]) &
            (mitoRatio < Max_mito_ratio[i])
    )
}
rm(i)

### Calculate Doublets
# https://github.com/chris-mcginnis-ucsf/DoubletFinder
# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Predict_doublets

Seurat_doublets_finder <- Seurat_list_Wang

sweep_res_list_Wang <- list()
sweep_stats_list_Wang <- list()
bcmvn_list_Wang <- list()
pK <- list()
BCmetric <- list()
pK_choose <- list()
nExp <- list()
per_doubl <- c(0.064, 0.068, 0.064, 0.068, 0.076, 0.072) # Calculated looking at 10X table

# pK Identification (no ground-truth)
for (i in seq_along(Seurat_doublets_finder)) {
    Seurat_doublets_finder[[i]] <- NormalizeData(Seurat_doublets_finder[[i]])
    Seurat_doublets_finder[[i]] <- FindVariableFeatures(Seurat_doublets_finder[[i]], selection.method = "vst", nfeatures = 2000)
    Seurat_doublets_finder[[i]] <- ScaleData(Seurat_doublets_finder[[i]])
    Seurat_doublets_finder[[i]] <- RunPCA(Seurat_doublets_finder[[i]], npcs = 20)
    Seurat_doublets_finder[[i]] <- RunUMAP(Seurat_doublets_finder[[i]], dims = 1:10)
    nExp[[i]] <- round(ncol(Seurat_doublets_finder[[i]]) * per_doubl[i])

    sweep_res_list_Wang[[i]] <- paramSweep_v3(Seurat_doublets_finder[[i]], PCs = 1:10, sct = FALSE, num.cores = 12)
    sweep_stats_list_Wang[[i]] <- summarizeSweep(sweep_res_list_Wang[[i]], GT = FALSE)
    bcmvn_list_Wang[[i]] <- find.pK(sweep_stats_list_Wang[[i]])
    dev.off()

    pK[[i]] <- as.numeric(as.character(bcmvn_list_Wang[[i]]$pK))
    BCmetric[[i]] <- bcmvn_list_Wang[[i]]$BCmetric
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

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Exploratory_analysis/Doublet_Finder_pK_calculation.pdf")
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

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Exploratory_analysis/Doublet_Finder_Classification.pdf")
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

    Seurat_list_Wang[[i]] <- subset(Seurat_list_Wang[[i]], cells = Singlets)
    print(Seurat_list_Wang[[i]])

    Sys.sleep(2)

    print(paste0(Patient_list[i], " DONE"))
    cat(paste0("----------------------------------------------------------------------------\n"))
}

### Subset by CAR+ cells
for (i in seq_along(Seurat_list_Wang)){
    Seurat_list_Wang[[i]] <- subset(x = Seurat_list_Wang[[i]], subset = `IFNB1-5LTR` > 0)
}

### Re-doing of exploratory analysis graphs after QC
# Group again all Products in a single metadata table
merged_metadata_df <- list()
for (cont in seq_along(Seurat_list_Wang)) {
    merged_metadata_df[[cont]] <- Seurat_list_Wang[[cont]]@meta.data
}
rm(cont)
merged_metadata_df <- bind_rows(merged_metadata_df)

# Visualize the number of cell counts per sample
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Exploratory_analysis/Post_QC_Ncells.pdf")
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
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Exploratory_analysis/Post_QC_UMIs_per_cell.pdf")
merged_metadata_df %>%
    ggplot(aes(color = Product, x = nCount_RNA, fill = Product)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    ylab("Cell density") +
    ggtitle("Number UMIs/transcripts per cell post-QC")
dev.off()

# Visualize the distribution of genes detected per cell
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Exploratory_analysis/Post_QC_Genes_per_cell.pdf")
merged_metadata_df %>%
    ggplot(aes(color = Product, x = nFeature_RNA, fill = Product)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    scale_x_log10() +
    ggtitle("Genes per cell post-QC")
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per cell
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Exploratory_analysis/Post_QC_Mito_ratio_per_cell.pdf")
merged_metadata_df %>%
    ggplot(aes(color = Product, x = mitoRatio, fill = Product)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    ggtitle("Mitocondrial ratio per cell post-QC")
dev.off()

# Violin plots
p_Vln <- list()
for (cont in seq_along(Seurat_list_Wang)) {
    p_Vln[[cont]] <- VlnPlot(Seurat_list_Wang[[cont]], features = c("nCount_RNA", "nFeature_RNA", "mitoRatio"), col = colores[cont], group.by = "Product", pt.size = 0)
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Exploratory_analysis/Post_QC_Vln_plots.pdf")
p_Vln
dev.off()

rm(cont, p_Vln)

# nGenes vs nUMI graph (Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs)
p1 <- list()
for (i in seq_along(Seurat_list_Wang)) {
    p1[[i]] <- ggplot(data = Seurat_list_Wang[[i]]@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = mitoRatio, grou)) +
        geom_point(size = 0.7) +
        scale_color_gradientn(colors = c("black", "blue", "green2", "red", "yellow")) +
        stat_smooth(method = lm) +
        scale_x_log10() +
        scale_y_log10() +
        theme_classic() +
        ggtitle(paste0(Patient_list[i]))
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Exploratory_analysis/Post_QC_nGenes_vs_nUMI.pdf")
p1
dev.off()

rm(i, p1)

setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/RDS")
saveRDS(Seurat_list_Wang, file = "PostQC_CellRanger_Wang_RDS.rds")


##### Normalization #####

# rm(list = ls())
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/RDS")
Seurat_list_Wang <- readRDS("PostQC_CellRanger_Wang_RDS.rds")
load("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Codigo/Gene_Markers_Info/cycle.rda")

for (i in seq_along(Seurat_list_Wang)) {
    Seurat_list_Wang[[i]] <- NormalizeData(Seurat_list_Wang[[i]])

    # Score cells for cell cycle
    Seurat_list_Wang[[i]] <- CellCycleScoring(Seurat_list_Wang[[i]],
        g2m.features = g2m_genes,
        s.features = s_genes)

    Seurat_list_Wang[[i]] <- FindVariableFeatures(Seurat_list_Wang[[i]],
        selection.method = "vst",
        nfeatures = 2000)

    Seurat_list_Wang[[i]] <- ScaleData(Seurat_list_Wang[[i]])
}
rm(i)

# Mark this strange population of cells with lower than 1000 features, for further analysis
for (i in seq_along(Seurat_list_Wang)) {
    Seurat_list_Wang[[i]]$Marked_cells <- (Seurat_list_Wang[[i]]$nFeature_RNA < 1000)
}

saveRDS(Seurat_list_Wang, file = "Normalized_CellRanger_Wang_RDS.rds")


##### Merge ---- Worst case scenario #####
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/RDS")
Seurat_list_Wang <- readRDS("Normalized_CellRanger_Wang_RDS.rds")

Seurat_list_Wang_2 <- merge(x = Seurat_list_Wang[[1]], y = c(Seurat_list_Wang[[2]], Seurat_list_Wang[[3]],
    Seurat_list_Wang[[4]], Seurat_list_Wang[[5]], Seurat_list_Wang[[6]]
))

Seurat_list_Wang <- Seurat_list_Wang_2
rm(Seurat_list_Wang_2)

Seurat_list_Wang <- FindVariableFeatures(Seurat_list_Wang,
    selection.method = "vst",
    nfeatures = 2000
)

Seurat_list_Wang <- ScaleData(Seurat_list_Wang)

# Run PCA
Seurat_list_Wang <- RunPCA(object = Seurat_list_Wang, reduction.name = "pca_wo_integ")

# ElbowPlot(Seurat_list_Wang, ndims = 50, reduction = "pca_wo_integ")
# dev.off()

# Run UMAP
Seurat_list_Wang <- RunUMAP(Seurat_list_Wang,
    dims = 1:30,
    reduction = "pca_wo_integ",
    reduction.name = "umap_wo_integ"
)

# Plot UMAP
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Worst_case_scenario/WO_integ_Seurat.pdf")
DimPlot(Seurat_list_Wang,
    group.by = "Product",
    reduction = "umap_wo_integ"
)
dev.off()

# Quality control metrics
p <- FeaturePlot(Seurat_list_Wang,
    reduction = "umap_wo_integ",
    features = c("nCount_RNA", "nFeature_RNA", "log10GenesPerUMI", "mitoRatio"),
    order = TRUE
)

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Wang_et_al/Plots/Worst_case_scenario/QC_metrics_WO_integ_Seurat.pdf")
p + plot_annotation(title = paste0(Seurat_list_Wang@meta.data$orig.ident %>% unique()), theme = theme(plot.title = element_text(size = 16)))
dev.off()

rm(p)

##### END OF SCRIPT #####