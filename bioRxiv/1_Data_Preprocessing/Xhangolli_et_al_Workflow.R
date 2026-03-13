###############################################################################
###############################################################################

# Program: Xhangolli_et_al_Workflow.R
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
library(data.table)
library(Matrix)
library(harmony)
library(VennDiagram)
library(patchwork)


##### Set seed #####
set.seed(2504)

###### Load data from Dropseq ######
Seurat_list_Xhangolli <- list()

setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Count_Matrices/Dropseq")
Patient_list <- dir()
Patient_list

for (contador in Patient_list) {
    Dropseq_dir <- paste0("./", contador, "/")
    tt <- fread(paste0(Dropseq_dir, contador, "_dge_matrix.txt.gz"))
    Dropseq_table <- Matrix(as.matrix(tt[, -1, with=FALSE]), sparse=TRUE)
    rownames(Dropseq_table) <- tt$GENE
    is_cell <- DropletUtils::emptyDrops(Dropseq_table)$FDR <= 0.001
    cell_counts <- Dropseq_table[, which(is_cell), drop = FALSE]
    Seurat_list_Xhangolli[[contador]] <- CreateSeuratObject(counts = cell_counts, project = "Xhangolli_et_al", min.features = 200, min.cells = 3)
    Seurat_list_Xhangolli[[contador]]$Product <- contador
    print(Seurat_list_Xhangolli[[contador]])
    print(paste0(contador, " DONE"))
}

rm(contador, Dropseq_dir, Dropseq_table, cell_counts, is_cell)


##### Quality control #####

# Add more metadata: novelty score and mito ratio
for (contador in seq_along(Seurat_list_Xhangolli)) {
    Seurat_list_Xhangolli[[contador]]$log10GenesPerUMI <- log10(Seurat_list_Xhangolli[[contador]]$nFeature_RNA) / log10(Seurat_list_Xhangolli[[contador]]$nCount_RNA)

    Seurat_list_Xhangolli[[contador]]$mitoRatio <- (PercentageFeatureSet(object = Seurat_list_Xhangolli[[contador]], pattern = "^MT-")) / 100

    Seurat_list_Xhangolli[[contador]]$log10GenesPerUMI <- log10(Seurat_list_Xhangolli[[contador]]$nFeature_RNA) / log10(Seurat_list_Xhangolli[[contador]]$nCount_RNA)
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
    for (contador in seq_along(Seurat_list_Xhangolli)) {
        if (min(Seurat_list_Xhangolli[[contador]]$nCount_RNA) >= minCov) {
            countLOW[[contador]] <- min(Seurat_list_Xhangolli[[contador]]$nCount_RNA)
        } else {
            countLOW[[contador]] <- quantile(Seurat_list_Xhangolli[[contador]]$nCount_RNA, prob = 0.01)
        }
        countHIGH[[contador]] <- quantile(Seurat_list_Xhangolli[[contador]]$nCount_RNA, prob = 0.99)
        featureLOW[[contador]] <- quantile(Seurat_list_Xhangolli[[contador]]$nFeature_RNA, prob = 0.01)
        featureHIGH[[contador]] <- quantile(Seurat_list_Xhangolli[[contador]]$nFeature_RNA, prob = 0.99)

        ## subset
        Seurat_list_Xhangolli[[contador]] <- subset(Seurat_list_Xhangolli[[contador]], subset = (nFeature_RNA > featureLOW[[contador]]) &
            (nFeature_RNA < featureHIGH[[contador]]) & (nCount_RNA > countLOW[[contador]]) & (nCount_RNA < countHIGH[[contador]]) & (mitoRatio < 0.20))
    }
    Pre_filter <- "DONE"
}
rm(minCov)

### Group all Products in a single metadata table
merged_metadata_df <- list()
for (cont in seq_along(Seurat_list_Xhangolli)) {
    merged_metadata_df[[cont]] <- Seurat_list_Xhangolli[[cont]]@meta.data
}
rm(cont)
merged_metadata_df <- bind_rows(merged_metadata_df)

### Exploratory analysis
# Values
# Order: Control, Raji_stim_1, Raji_stim_2

Min_feat <- c(300, 275, 275)
Max_feat <- c(5000, 4750, 4750)
Min_counts <- c(350, 350, 350)
Max_counts <- c(25000, 22500, 22500)
Max_mito_ratio <- c(0.14, 0.10, 0.10)

colores <- hue_pal()(length(Seurat_list_Xhangolli))

# Visualize the number of cell counts per sample
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Exploratory_analysis/Pre_QC_Ncells.pdf")
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
sce_Seurat_list_Xhangolli <- list()
metadata_bcrank <- list()
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Exploratory_analysis/Pre_QC_Knee_plots.pdf")
for (i in Patient_list) {
    counter <- counter + 1
    sce_Seurat_list_Xhangolli[[i]] <- as.SingleCellExperiment(Seurat_list_Xhangolli[[i]])
    bcrank <- DropletUtils::barcodeRanks(counts(sce_Seurat_list_Xhangolli[[i]]))

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
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Exploratory_analysis/Pre_QC_UMIs_per_cell.pdf")
merged_metadata_df %>%
    ggplot(aes(color = Product, x = nCount_RNA, fill = Product)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    ylab("Cell density") +

    # Adjust thresholds per dataset
    geom_vline(xintercept = Min_counts[1], col = colores[1], alpha = 0.6, linetype = "dotted") +
    geom_vline(xintercept = Max_counts[1], col = colores[1], alpha = 0.6, linetype = "dotted") +
    geom_vline(xintercept = Min_counts[2], col = colores[2], alpha = 0.6, linetype = "dotted") +
    geom_vline(xintercept = Max_counts[2], col = colores[2], alpha = 0.6, linetype = "dotted") +
    geom_vline(xintercept = Min_counts[3], col = colores[3], alpha = 0.6, linetype = "dotted") +
    geom_vline(xintercept = Max_counts[3], col = colores[3], alpha = 0.6, linetype = "dotted") +
    ggtitle("Number UMIs/transcripts per cell")
dev.off()

# Visualize the distribution of genes detected per cell
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Exploratory_analysis/Pre_QC_Genes_per_cell.pdf")
merged_metadata_df %>%
    ggplot(aes(color = Product, x = nFeature_RNA, fill = Product)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    scale_x_log10() +

    # Adjust thresholds per dataset
    geom_vline(xintercept = Min_feat[1], col = colores[1], alpha = 0.6, linetype = "dotted") +
    geom_vline(xintercept = Max_feat[1], col = colores[1], alpha = 0.6, linetype = "dotted") +
    geom_vline(xintercept = Min_feat[2], col = colores[2], alpha = 0.6, linetype = "dotted") +
    geom_vline(xintercept = Max_feat[2], col = colores[2], alpha = 0.6, linetype = "dotted") +
    geom_vline(xintercept = Min_feat[3], col = colores[3], alpha = 0.6, linetype = "dotted") +
    geom_vline(xintercept = Max_feat[3], col = colores[3], alpha = 0.6, linetype = "dotted") +
    ggtitle("Genes per cell")
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per cell
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Exploratory_analysis/Pre_QC_Mito_ratio_per_cell.pdf")
merged_metadata_df %>%
    ggplot(aes(color = Product, x = mitoRatio, fill = Product)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +

    # Adjust thresholds per dataset
    geom_vline(xintercept = Max_mito_ratio[1], col = colores[1], alpha = 0.6, linetype = "dotted") +
    geom_vline(xintercept = Max_mito_ratio[2], col = colores[2], alpha = 0.6, linetype = "dotted") +
    geom_vline(xintercept = Max_mito_ratio[3], col = colores[3], alpha = 0.6, linetype = "dotted") +
    ggtitle("Mitocondrial ratio per cell")
dev.off()

# Visualize the distribution of complexity per cell
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Exploratory_analysis/Pre_QC_Complexity_per_cell.pdf")
merged_metadata_df %>%
    ggplot(aes(color = Product, x = log10GenesPerUMI, fill = Product)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    ggtitle("Complexity per cell")
dev.off()

# Violin plots
p_Vln <- list()
for (cont in seq_along(Seurat_list_Xhangolli)) {
    p_Vln[[cont]] <- VlnPlot(Seurat_list_Xhangolli[[cont]], features = c("nCount_RNA", "nFeature_RNA", "mitoRatio"), col = colores[cont], group.by = "Product", pt.size = 0)
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Exploratory_analysis/Pre_QC_Vln_plots.pdf")
p_Vln
dev.off()

rm(cont, p_Vln)

# nGenes vs nUMI graph (Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs)
p1 <- list()
for (i in seq_along(Seurat_list_Xhangolli)) {
    p1[[i]] <- ggplot(data = Seurat_list_Xhangolli[[i]]@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = mitoRatio)) +
        geom_point(size = 0.7) +
        scale_color_gradientn(colors = c("black", "blue", "green2", "red", "yellow")) +
        stat_smooth(method = lm) +
        scale_x_log10() +
        scale_y_log10() +
        theme_classic() +
        ggtitle(paste0(Patient_list[i])) +
        geom_vline(xintercept = Min_counts[i], col = "black") +
        geom_vline(xintercept = Max_counts[i], col = "black") +
        geom_hline(yintercept = Min_feat[i], col = "black") +
        geom_hline(yintercept = Max_feat[i], col = "black")
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Exploratory_analysis/Pre_QC_nGenes_vs_nUMI.pdf")
p1
dev.off()

rm(i, p1)

# mitoRatio vs nUMI graph
p1 <- list()
for (i in seq_along(Seurat_list_Xhangolli)) {
    p1[[i]] <- ggplot(data = Seurat_list_Xhangolli[[i]]@meta.data, aes(x = nCount_RNA, y = mitoRatio)) +
        geom_point(size = 0.7) +
        theme_classic() +
        ggtitle(paste0(Patient_list[i]))
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Exploratory_analysis/Pre_QC_mitoRatio_vs_nUMI.pdf")
p1
dev.off()

rm(i, p1)

### Subset samples by determined quality control thresholds
for (i in seq_along(Seurat_list_Xhangolli)) {
    Seurat_list_Xhangolli[[i]] <- subset(Seurat_list_Xhangolli[[i]],
        subset = (nFeature_RNA > Min_feat[i]) &
            (nFeature_RNA < Max_feat[i]) &
            (mitoRatio < Max_mito_ratio[i])
    )
}
rm(i)

### Calculate Doublets
# https://github.com/chris-mcginnis-ucsf/DoubletFinder
# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Predict_doublets

Seurat_doublets_finder <- Seurat_list_Xhangolli

sweep_res_list_Xhangolli <- list()
sweep_stats_list_Xhangolli <- list()
bcmvn_list_Xhangolli <- list()
pK <- list()
BCmetric <- list()
pK_choose <- list()
nExp <- list()
per_doubl <- c(0.02, 0.02, 0.02) # scFTD seq method -- demonstrated a doublet rate of less than 2%

# pK Identification (no ground-truth)
for (i in seq_along(Seurat_doublets_finder)) {
    Seurat_doublets_finder[[i]] <- NormalizeData(Seurat_doublets_finder[[i]])
    Seurat_doublets_finder[[i]] <- FindVariableFeatures(Seurat_doublets_finder[[i]], selection.method = "vst", nfeatures = 2000)
    Seurat_doublets_finder[[i]] <- ScaleData(Seurat_doublets_finder[[i]])
    Seurat_doublets_finder[[i]] <- RunPCA(Seurat_doublets_finder[[i]], npcs = 20)
    Seurat_doublets_finder[[i]] <- RunUMAP(Seurat_doublets_finder[[i]], dims = 1:10)
    nExp[[i]] <- round(ncol(Seurat_doublets_finder[[i]]) * per_doubl[i])

    sweep_res_list_Xhangolli[[i]] <- paramSweep_v3(Seurat_doublets_finder[[i]], PCs = 1:10, sct = FALSE, num.cores = 12)
    sweep_stats_list_Xhangolli[[i]] <- summarizeSweep(sweep_res_list_Xhangolli[[i]], GT = FALSE)
    bcmvn_list_Xhangolli[[i]] <- find.pK(sweep_stats_list_Xhangolli[[i]])
    dev.off()

    pK[[i]] <- as.numeric(as.character(bcmvn_list_Xhangolli[[i]]$pK))
    BCmetric[[i]] <- bcmvn_list_Xhangolli[[i]]$BCmetric
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

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Exploratory_analysis/Doublet_Finder_pK_calculation.pdf")
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

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Exploratory_analysis/Doublet_Finder_Classification.pdf")
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

    Seurat_list_Xhangolli[[i]] <- subset(Seurat_list_Xhangolli[[i]], cells = Singlets)
    print(Seurat_list_Xhangolli[[i]])

    Sys.sleep(2)

    print(paste0(Patient_list[i], " DONE"))
    cat(paste0("----------------------------------------------------------------------------\n"))
}

### Re-doing of exploratory analysis graphs after QC
# Group again all Products in a single metadata table
merged_metadata_df <- list()
for (cont in seq_along(Seurat_list_Xhangolli)) {
    merged_metadata_df[[cont]] <- Seurat_list_Xhangolli[[cont]]@meta.data
}
rm(cont)
merged_metadata_df <- bind_rows(merged_metadata_df)

# Visualize the number of cell counts per sample
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Exploratory_analysis/Post_QC_Ncells.pdf")
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
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Exploratory_analysis/Post_QC_UMIs_per_cell.pdf")
merged_metadata_df %>%
    ggplot(aes(color = Product, x = nCount_RNA, fill = Product)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    ylab("Cell density") +
    ggtitle("Number UMIs/transcripts per cell post-QC")
dev.off()

# Visualize the distribution of genes detected per cell
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Exploratory_analysis/Post_QC_Genes_per_cell.pdf")
merged_metadata_df %>%
    ggplot(aes(color = Product, x = nFeature_RNA, fill = Product)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    scale_x_log10() +
    ggtitle("Genes per cell post-QC")
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per cell
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Exploratory_analysis/Post_QC_Mito_ratio_per_cell.pdf")
merged_metadata_df %>%
    ggplot(aes(color = Product, x = mitoRatio, fill = Product)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    ggtitle("Mitocondrial ratio per cell post-QC")
dev.off()

# Violin plots
p_Vln <- list()
for (cont in seq_along(Seurat_list_Xhangolli)) {
    p_Vln[[cont]] <- VlnPlot(Seurat_list_Xhangolli[[cont]], features = c("nCount_RNA", "nFeature_RNA", "mitoRatio"), col = colores[cont], group.by = "Product", pt.size = 0)
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Exploratory_analysis/Post_QC_Vln_plots.pdf")
p_Vln
dev.off()

rm(cont, p_Vln)

# nGenes vs nUMI graph (Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs)
p1 <- list()
for (i in seq_along(Seurat_list_Xhangolli)) {
    p1[[i]] <- ggplot(data = Seurat_list_Xhangolli[[i]]@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = mitoRatio, grou)) +
        geom_point(size = 0.7) +
        scale_color_gradientn(colors = c("black", "blue", "green2", "red", "yellow")) +
        stat_smooth(method = lm) +
        scale_x_log10() +
        scale_y_log10() +
        theme_classic() +
        ggtitle(paste0(Patient_list[i]))
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Exploratory_analysis/Post_QC_nGenes_vs_nUMI.pdf")
p1
dev.off()

rm(i, p1)

setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/RDS")
saveRDS(Seurat_list_Xhangolli, file = "PostQC_CellRanger_Xhangolli_RDS.rds")

##### Normalization #####

# rm(list = ls())
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/RDS")
Seurat_list_Xhangolli <- readRDS("PostQC_CellRanger_Xhangolli_RDS.rds")
load("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Codigo/Gene_Markers_Info/cycle.rda")

for (i in seq_along(Seurat_list_Xhangolli)) {
    Seurat_list_Xhangolli[[i]] <- NormalizeData(Seurat_list_Xhangolli[[i]])

    # Score cells for cell cycle
    Seurat_list_Xhangolli[[i]] <- CellCycleScoring(Seurat_list_Xhangolli[[i]],
        g2m.features = g2m_genes,
        s.features = s_genes
    )
}
rm(i)


##### Discarding residual cancer cells described in the original article
Cancerous_cell_detection <- Seurat_list_Xhangolli
p1 <- list()
markers <- list()
for (i in seq_along(Cancerous_cell_detection)) {
    Cancerous_cell_detection[[i]] <- ScaleData(object = Cancerous_cell_detection[[i]])
    Cancerous_cell_detection[[i]] <- RunPCA(object = Cancerous_cell_detection[[i]])
    Cancerous_cell_detection[[i]] <- FindNeighbors(object = Cancerous_cell_detection[[i]])
    Cancerous_cell_detection[[i]] <- FindClusters(object = Cancerous_cell_detection[[i]])
    Cancerous_cell_detection[[i]] <- RunUMAP(object = Cancerous_cell_detection[[i]], dims = 1:10)
    markers[[i]] <- FindAllMarkers(
        object = Cancerous_cell_detection[[i]],
        only.pos = TRUE,
        logfc.threshold = 0.25
    )
    p1[[i]] <- cowplot::plot_grid(
        nrow = 2, DimPlot(object = Cancerous_cell_detection[[i]], reduction = "umap"), FeaturePlot(object = Cancerous_cell_detection[[i]], reduction = "umap", features = c("IGHM", "CD3D", "CD19"), pt.size = 0.4, order = TRUE, label = TRUE)
    )
}

setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Cancer_cells_Removal")
pdf("Cancer_cell_individual_detection.pdf")
p1
dev.off()

Cancer_cells_Raji_stim_1 <- paste0("Raji_stim_1_", names(Idents(Cancerous_cell_detection[[2]]))[Idents(Cancerous_cell_detection[[2]]) == 6])
Cancer_cells_Raji_stim_2 <- paste0("Raji_stim_2_", names(Idents(Cancerous_cell_detection[[3]]))[Idents(Cancerous_cell_detection[[3]]) == 8])


## Merge and integrate with Harmony - To see this cancerous population together
Seurat_merged <- merge(x = Cancerous_cell_detection[[1]], y = c(
    Cancerous_cell_detection[[2]], Cancerous_cell_detection[[3]]),
    add.cell.ids = c("Control", "Raji_stim_1", "Raji_stim_2"))

Seurat_merged <- FindVariableFeatures(Seurat_merged,
        selection.method = "vst",
        nfeatures = 2000
)

Seurat_merged <- ScaleData(Seurat_merged)

Seurat_merged <- RunPCA(object = Seurat_merged)

Seurat_merged <- RunHarmony(Seurat_merged, "Product")

Seurat_merged <- RunUMAP(Seurat_merged,
        dims = 1:50,
        reduction = "harmony",
        reduction.name = "umap_harmony"
)

Seurat_merged <- FindNeighbors(object= Seurat_merged, reduction = "harmony", dims = 1:20)

Seurat_merged <- FindClusters(Seurat_merged)

pdf("Cancer_cell_integrated_detection.pdf")
cowplot::plot_grid(
        nrow = 2, DimPlot(object = Seurat_merged, reduction = "umap_harmony"), FeaturePlot(object = Seurat_merged, reduction = "umap_harmony", features = c("IGHM", "CD3D", "CD19"), pt.size = 0.4, order = TRUE, label = TRUE)
    )
dev.off()

Cancer_cells_integrated <- names(Idents(Seurat_merged))[Idents(Seurat_merged) == 10]
Cancer_cells_no_integrated <- c(Cancer_cells_Raji_stim_1, Cancer_cells_Raji_stim_2)

intersect(Cancer_cells_integrated, Cancer_cells_no_integrated)

length(Cancer_cells_integrated)
length(Cancer_cells_no_integrated)
length(intersect(Cancer_cells_integrated, Cancer_cells_no_integrated))

venn.diagram(x= list(Cancer_cells_integrated, Cancer_cells_no_integrated), filename = 'Venn_diagramm.png', category.names = c("Integrated" , "No_integrated"), output=TRUE)

discard_cells <- unique(c(Cancer_cells_integrated, Cancer_cells_no_integrated))
length(discard_cells)

Seurat_merged <- subset(Seurat_merged, cells = discard_cells, invert = TRUE)

table(Seurat_merged@meta.data$Product)

## Re-scale the data
Seurat_merged <- ScaleData(Seurat_merged)

Seurat_merged <- RunPCA(object = Seurat_merged, reduction.name = "pca_2")

Seurat_merged <- RunHarmony(Seurat_merged, "Product", reduction = "pca_2", reduction.save = "harmony_2")

Seurat_merged <- RunUMAP(Seurat_merged,
        dims = 1:50,
        reduction = "harmony_2",
        reduction.name = "umap_harmony_2"
)

Seurat_merged <- FindNeighbors(object= Seurat_merged, reduction = "harmony_2", dims = 1:20)

Seurat_merged <- FindClusters(Seurat_merged)

# DimPlot(Seurat_merged, reduction = "umap_harmony_2", group.by = "Product", pt.size = .4)

pdf("UMAP_post_cancer_cells_removal.pdf")
cowplot::plot_grid(
        nrow = 2, DimPlot(object = Seurat_merged, reduction = "umap_harmony_2"), FeaturePlot(object = Seurat_merged, reduction = "umap_harmony_2", features = c("IGHM", "CD3D", "CD19"), pt.size = 0.4, order = TRUE, label = TRUE)
    )
dev.off()

Seurat_list_Xhangolli <- NULL
Seurat_list_Xhangolli <- SplitObject(Seurat_merged, split.by = "Product")

for(i in seq_along(Seurat_list_Xhangolli)){
    Seurat_list_Xhangolli[[i]]@meta.data$RNA_snn_res.0.8 <- NULL
    Seurat_list_Xhangolli[[i]]@meta.data$seurat_clusters <- NULL
}

setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/RDS")
saveRDS(Seurat_list_Xhangolli, file = "Normalized_CellRanger_Xhangolli_RDS.rds")


##### Merge ---- Worst case scenario #####
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/RDS")
Seurat_list_Xhangolli <- readRDS("Normalized_CellRanger_Xhangolli_RDS.rds")

Seurat_list_Xhangolli_2 <- merge(x = Seurat_list_Xhangolli[[1]], y = c(Seurat_list_Xhangolli[[2]], Seurat_list_Xhangolli[[3]]))

Seurat_list_Xhangolli <- Seurat_list_Xhangolli_2
rm(Seurat_list_Xhangolli_2)

Seurat_list_Xhangolli <- FindVariableFeatures(Seurat_list_Xhangolli,
    selection.method = "vst",
    nfeatures = 2000
)

Seurat_list_Xhangolli <- ScaleData(Seurat_list_Xhangolli)

# Run PCA
Seurat_list_Xhangolli <- RunPCA(object = Seurat_list_Xhangolli, reduction.name = "pca_wo_integ")

# ElbowPlot(Seurat_list_Xhangolli, ndims = 50, reduction = "pca_wo_integ")
# dev.off()

# Run UMAP
Seurat_list_Xhangolli <- RunUMAP(Seurat_list_Xhangolli,
    dims = 1:25,
    reduction = "pca_wo_integ",
    reduction.name = "umap_wo_integ"
)

# Plot UMAP
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Worst_case_scenario/WO_integ_Seurat.pdf")
DimPlot(Seurat_list_Xhangolli,
    group.by = "Product",
    reduction = "umap_wo_integ"
)
dev.off()

# Quality control metrics
p <- FeaturePlot(Seurat_list_Xhangolli,
    reduction = "umap_wo_integ",
    features = c("nCount_RNA", "nFeature_RNA", "log10GenesPerUMI", "mitoRatio"),
    order = TRUE
)

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Xhangolli_et_al/Plots/Worst_case_scenario/QC_metrics_WO_integ_Seurat.pdf")
p + plot_annotation(title = paste0(Seurat_list_Xhangolli@meta.data$orig.ident %>% unique()), theme = theme(plot.title = element_text(size = 16)))
dev.off()

rm(p)

##### END OF SCRIPT #####