###############################################################################
###############################################################################

# Program: Haradvala_et_al_Workflow.R
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
if (FALSE) {
    setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Count_Matrices/Cell_Ranger/")
    Patient_list <- gsub("run_", "", dir()[dir() %>% grep(pattern = "run_")])

    if (FALSE) {
        Seurat_list_Haradvala <- list()
        for (contador in Patient_list) {
            cell_ranger_dir <- paste0("./run_", contador, "/outs/filtered_feature_bc_matrix")
            cell_ranger_data <- Read10X(data.dir = cell_ranger_dir)
            Seurat_list_Haradvala[[contador]] <- CreateSeuratObject(counts = cell_ranger_data, project = "Haradvala_et_al", min.features = 200, min.cells = 3)
            Seurat_list_Haradvala[[contador]]$Product <- contador

            print(Seurat_list_Haradvala[[contador]])
            print(paste0(contador, " DONE"))
        }

        rm(contador, cell_ranger_dir, cell_ranger_data)

        saveRDS(Seurat_list_Haradvala, "Seurat_list_Haradvala.RDS")
    } else {
        Seurat_list_Haradvala <- readRDS("Seurat_list_Haradvala.RDS")
    }


    ### Re-Identification of all samples
    for (contador in Patient_list) {
        Seurat_list_Haradvala[[contador]] <- RenameCells(Seurat_list_Haradvala[[contador]], new.names = gsub("-1", "", colnames(Seurat_list_Haradvala[[contador]])))
    }

    setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Datasets/Haradhvala_et_al/Downloaded_data/Authors_Metadata_Reduced")
    Filtered_reduced_authors_meta <- read.csv("Haradvala_filtered_metadata.csv", row.names = 1)
    Filtered_reduced_authors_meta_splited <- split(Filtered_reduced_authors_meta, Filtered_reduced_authors_meta$channel)

    # I only want to keep the multiplexed ones
    Filtered_reduced_authors_meta_splited <- Filtered_reduced_authors_meta_splited[1:32]

    for (contador in names(Filtered_reduced_authors_meta_splited)) {
        rownames(Filtered_reduced_authors_meta_splited[[contador]]) <- gsub("-[0-9]+-[0-9]+$", "", rownames(Filtered_reduced_authors_meta_splited[[contador]]))
    }

    ## Sample pool data:
    # A1: Axi-R-17 Infusion, Tisa-R-22 Infusion, Tisa-N-23 Infusion
    # A2: Axi-R-15 Infusion, Axi-N-18 Infusion, Tisa-N-24 Infusion, Tisa-N-31 Infusion
    # A3: Axi-N-16 Infusion, Tisa-R-26 Infusion, Tisa-N-28 Infusion
    # A4: Axi-R-15 Baseline, Axi-R-17 Baseline, Tisa-R-22 Baseline, Tisa-N-23 Baseline
    # B1: Axi-R-13 Infusion, Axi-R-19 Infusion, Tisa-N-25 Baseline, Tisa-R-30 Infusion
    # B2: Axi-N-10 Infusion, Axi-R-11 Infusion, Tisa-N-29 Infusion, Tisa-R-32 Infusion
    # B3: Axi-N-07 Infusion, Axi-R-08 Infusion, Axi-N-09 Infusion, Tisa-N-27 Infusion
    # B4: Axi-R-13 Baseline, Axi-R-19 Baseline, Tisa-N-25 Infusion, Tisa-R-30 Baseline
    # B5: Axi-R-08 Baseline, Axi-N-09 Baseline, Axi-N-10 Baseline, Axi-R-11 Baseline
    # B6: Axi-N-06 Baseline, Axi-N-18 Baseline, Tisa-N-24 Baseline, Tisa-N-31 Baseline
    # C1: Axi-N-16 D7, Axi-R-17 D7, Tisa-N-28 D7
    # C2: Axi-N-06 D7, Tisa-R-26 D7, Tisa-N-31 D7
    # E1: Axi-N-07 D7, Axi-R-15 D7, Tisa-N-27 D7, Tisa-R-32 D7
    # E2: Axi-R-11 D7, Axi-N-18 D7, Tisa-R-22 D7, Tisa-R-30 D7
    # G1: Axi-N-10 D7, Axi-R-19 D7, Tisa-N-23 D7, Tisa-N-25 D7
    # G2: Axi-R-08 D7, Axi-N-09 D7, Tisa-N-24 D7, Tisa-N-29 D7

    Seurat_list_Haradvala_reintification <- Seurat_list_Haradvala[1:32]

    # Loop through each element in the list

    for (i in seq_along(Seurat_list_Haradvala_reintification)) {
        # Extract necessary columns
        table1 <- Seurat_list_Haradvala_reintification[[i]]@meta.data
        table2 <- Filtered_reduced_authors_meta_splited[[i]]

        # Match row names
        matching_rows <- intersect(rownames(table1), rownames(table2))

        # Add barcode and timepoint columns from table2 to table1 based on shared row names
        table1$barcode <- NA # Initialize a barcode column in table1
        table1$timepoint <- NA # Initialize a timepoint column in table1

        # Only add values from table2 to table1 where row names match
        table1[matching_rows, "barcode"] <- table2[matching_rows, "barcode"]
        table1[matching_rows, "timepoint"] <- table2[matching_rows, "timepoint"]

        # Returning variables back
        Seurat_list_Haradvala_reintification[[i]]@meta.data <- table1
        Filtered_reduced_authors_meta_splited[[i]] <- table2

        # Remove rows with NA in the 'barcode' column
        Seurat_list_Haradvala_reintification[[i]] <- subset(Seurat_list_Haradvala_reintification[[i]], subset = barcode != "NA")

        # Merge columns
        Seurat_list_Haradvala_reintification[[i]]$barcode_timepoint <- paste(Seurat_list_Haradvala_reintification[[i]]$barcode, Seurat_list_Haradvala_reintification[[i]]$timepoint, sep = "_")

        # Rename the cell using renameCell function
        new_cell_name <- paste(colnames(Seurat_list_Haradvala_reintification[[i]]), i, sep = "-")
        Seurat_list_Haradvala_reintification[[i]] <- RenameCells(Seurat_list_Haradvala_reintification[[i]], new.name = new_cell_name)

        # DONE
        print(paste0("DONE nº", i))
    }

    # Merge in one object to split after bt the newly created column
    Seurat_list_Haradvala_reintification_merged <- merge(
        x = Seurat_list_Haradvala_reintification[[1]],
        y = c(
            Seurat_list_Haradvala_reintification[[2]],
            Seurat_list_Haradvala_reintification[[3]],
            Seurat_list_Haradvala_reintification[[4]],
            Seurat_list_Haradvala_reintification[[5]],
            Seurat_list_Haradvala_reintification[[6]],
            Seurat_list_Haradvala_reintification[[7]],
            Seurat_list_Haradvala_reintification[[8]],
            Seurat_list_Haradvala_reintification[[9]],
            Seurat_list_Haradvala_reintification[[10]],
            Seurat_list_Haradvala_reintification[[11]],
            Seurat_list_Haradvala_reintification[[12]],
            Seurat_list_Haradvala_reintification[[13]],
            Seurat_list_Haradvala_reintification[[14]],
            Seurat_list_Haradvala_reintification[[15]],
            Seurat_list_Haradvala_reintification[[16]],
            Seurat_list_Haradvala_reintification[[17]],
            Seurat_list_Haradvala_reintification[[18]],
            Seurat_list_Haradvala_reintification[[19]],
            Seurat_list_Haradvala_reintification[[20]],
            Seurat_list_Haradvala_reintification[[21]],
            Seurat_list_Haradvala_reintification[[22]],
            Seurat_list_Haradvala_reintification[[23]],
            Seurat_list_Haradvala_reintification[[24]],
            Seurat_list_Haradvala_reintification[[25]],
            Seurat_list_Haradvala_reintification[[26]],
            Seurat_list_Haradvala_reintification[[27]],
            Seurat_list_Haradvala_reintification[[28]],
            Seurat_list_Haradvala_reintification[[29]],
            Seurat_list_Haradvala_reintification[[30]],
            Seurat_list_Haradvala_reintification[[31]],
            Seurat_list_Haradvala_reintification[[32]]
        )
    )

    # Split by the newly created column
    Seurat_list_Haradvala_final <- list()
    Seurat_list_Haradvala_final <- c(Seurat_list_Haradvala_final, SplitObject(Seurat_list_Haradvala_reintification_merged, split.by = "barcode_timepoint"))

    names(Seurat_list_Haradvala_final)

    Seurat_list_Haradvala_final_filtered_list <- Seurat_list_Haradvala_final[!grepl("Baseline", names(Seurat_list_Haradvala_final))]
    Seurat_list_Haradvala_final_filtered_list <- Seurat_list_Haradvala_final_filtered_list[!grepl("D7", names(Seurat_list_Haradvala_final_filtered_list))]
    names(Seurat_list_Haradvala_final_filtered_list)

    ## Change the names with these matchings
    # "Tisa-R-22_Infusion" - Har_Pat22_IP
    # "Tisa-N-23_Infusion" - Har_Pat23_IP
    # "Axi-R-17_Infusion" - Har_Pat17_IP
    # "Tisa-N-24_Infusion" - Har_Pat24_IP
    # "Tisa-N-31_Infusion" - Har_Pat31_IP
    # "Axi-R-15_Infusion" - Har_Pat15_IP
    # "Axi-N-18_Infusion" - Har_Pat18_IP
    # "Tisa-R-26_Infusion" - Har_Pat26_IP
    # "Tisa-N-28_Infusion" - Har_Pat28_IP
    # "Axi-N-16_Infusion" - Har_Pat16_IP
    # "Tisa-R-30_Infusion" - Har_Pat30_IP
    # "Axi-R-19_Infusion" - Har_Pat19_IP
    # "Axi-R-13_Infusion" - Har_Pat13_IP
    # "Tisa-R-32_Infusion" - Har_Pat32_IP
    # "Tisa-N-29_Infusion" - Har_Pat29_IP
    # "Axi-N-10_Infusion" - Har_Pat10_IP
    # "Axi-R-11_Infusion" - Har_Pat11_IP
    # "Tisa-N-27_Infusion" - Har_Pat27_IP
    # "Axi-N-09_Infusion" - Har_Pat9_IP
    # "Axi-R-08_Infusion" - Har_Pat8_IP
    # "Axi-N-07_Infusion" - Har_Pat7_IP
    # "Tisa-N-25_Infusion" - Har_Pat25_IP

    # Original names in Seurat object
    original_names <- names(Seurat_list_Haradvala_final_filtered_list)

    new_names_list <- c(
        "Har_Pat22_IP", "Har_Pat23_IP", "Har_Pat17_IP", "Har_Pat24_IP", "Har_Pat31_IP",
        "Har_Pat15_IP", "Har_Pat18_IP", "Har_Pat26_IP", "Har_Pat28_IP", "Har_Pat16_IP",
        "Har_Pat30_IP", "Har_Pat19_IP", "Har_Pat13_IP", "Har_Pat32_IP", "Har_Pat29_IP",
        "Har_Pat10_IP", "Har_Pat11_IP", "Har_Pat27_IP", "Har_Pat9_IP", "Har_Pat8_IP",
        "Har_Pat7_IP", "Har_Pat25_IP"
    )

    names(Seurat_list_Haradvala_final_filtered_list) <- new_names_list

    # Now the names in Seurat object are changed
    print(original_names)
    print(names(Seurat_list_Haradvala_final_filtered_list))

    # Add the other elements to the object
    Seurat_list_Haradvala_bis <- c(Seurat_list_Haradvala[33:76], Seurat_list_Haradvala_final_filtered_list)

    ## Order to have the same order as in the excel metadata - This is the order:
    # Har_Pat1_IP
    # Har_Pat2_D7
    # Har_Pat2_IP
    # Har_Pat3_IP
    # Har_Pat4_IP
    # Har_Pat4_D7
    # Har_Pat5_IP
    # Har_Pat6_D7
    # Har_Pat7_D7
    # Har_Pat7_IP
    # Har_Pat8_IP
    # Har_Pat8_D7
    # Har_Pat9_IP
    # Har_Pat9_D7
    # Har_Pat10_IP
    # Har_Pat10_D7
    # Har_Pat11_D7
    # Har_Pat11_IP
    # Har_Pat12_D7
    # Har_Pat12_IP
    # Har_Pat12_D14
    # Har_Pat13_IP
    # Har_Pat13_D7
    # Har_Pat14_D7
    # Har_Pat14_D14
    # Har_Pat14_IP
    # Har_Pat15_D7
    # Har_Pat15_IP
    # Har_Pat16_D7
    # Har_Pat16_IP
    # Har_Pat17_IP
    # Har_Pat17_D7
    # Har_Pat18_D7
    # Har_Pat18_IP
    # Har_Pat19_D7
    # Har_Pat19_IP
    # Har_Pat20_D7
    # Har_Pat20_IP
    # Har_Pat20_D14
    # Har_Pat21_D7
    # Har_Pat21_D14
    # Har_Pat21_IP
    # Har_Pat22_D7
    # Har_Pat22_IP
    # Har_Pat23_D7
    # Har_Pat23_IP
    # Har_Pat24_IP
    # Har_Pat24_D7
    # Har_Pat25_IP
    # Har_Pat25_D7
    # Har_Pat26_IP
    # Har_Pat26_D7
    # Har_Pat27_D7
    # Har_Pat27_IP
    # Har_Pat28_IP
    # Har_Pat28_D7
    # Har_Pat29_IP
    # Har_Pat29_IP_retreat
    # Har_Pat29_D7
    # Har_Pat29_D7_retreat
    # Har_Pat30_D7
    # Har_Pat30_IP
    # Har_Pat31_D7
    # Har_Pat31_IP
    # Har_Pat32_IP
    # Har_Pat32_D7

    # Define the desired order (names of Seurat objects)
    new_order <- c(
        "Har_Pat1_IP", "Har_Pat2_D7", "Har_Pat2_IP", "Har_Pat3_IP", "Har_Pat4_IP",
        "Har_Pat4_D7", "Har_Pat5_IP", "Har_Pat6_D7", "Har_Pat7_D7", "Har_Pat7_IP",
        "Har_Pat8_IP", "Har_Pat8_D7", "Har_Pat9_IP", "Har_Pat9_D7", "Har_Pat10_IP",
        "Har_Pat10_D7", "Har_Pat11_D7", "Har_Pat11_IP", "Har_Pat12_D7", "Har_Pat12_IP",
        "Har_Pat12_D14", "Har_Pat13_IP", "Har_Pat13_D7", "Har_Pat14_D7", "Har_Pat14_D14",
        "Har_Pat14_IP", "Har_Pat15_D7", "Har_Pat15_IP", "Har_Pat16_D7", "Har_Pat16_IP",
        "Har_Pat17_IP", "Har_Pat17_D7", "Har_Pat18_D7", "Har_Pat18_IP", "Har_Pat19_D7",
        "Har_Pat19_IP", "Har_Pat20_D7", "Har_Pat20_IP", "Har_Pat20_D14", "Har_Pat21_D7",
        "Har_Pat21_D14", "Har_Pat21_IP", "Har_Pat22_D7", "Har_Pat22_IP", "Har_Pat23_D7",
        "Har_Pat23_IP", "Har_Pat24_IP", "Har_Pat24_D7", "Har_Pat25_IP", "Har_Pat25_D7",
        "Har_Pat26_IP", "Har_Pat26_D7", "Har_Pat27_D7", "Har_Pat27_IP", "Har_Pat28_IP",
        "Har_Pat28_D7", "Har_Pat29_IP", "Har_Pat29_IP_retreat", "Har_Pat29_D7",
        "Har_Pat29_D7_retreat", "Har_Pat30_D7", "Har_Pat30_IP", "Har_Pat31_D7",
        "Har_Pat31_IP", "Har_Pat32_IP", "Har_Pat32_D7"
    )

    # Reorder the list of Seurat objects based on the new order
    reordered_seurat_list <- Seurat_list_Haradvala_bis[new_order]

    for (i in seq_along(reordered_seurat_list)) {
        reordered_seurat_list[[i]]@meta.data$Product <- names(reordered_seurat_list)[i]
    }

    Seurat_list_Haradvala <- NULL
    Seurat_list_Haradvala <- reordered_seurat_list

    # Save the object to save time with previous steps (Remember to put and IF:False to not execute all of this chunk of code)
    setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/")
    saveRDS(Seurat_list_Haradvala, "Seurat_list_Haradvala_inicio.RDS")
}

##### Quality control #####

## Read object created before
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/")
Seurat_list_Haradvala <- readRDS("Seurat_list_Haradvala_inicio.RDS")

## Count number of cells with CAR+ per sample to get the number of inicial cells pre-QC
Sum_CAR <- list()
Sum_Kymriah <- list()
Sum_Yescarta <- list()

for (contador in seq_along(Seurat_list_Haradvala)) {
    Sum_CAR[[contador]] <- sum(Seurat_list_Haradvala[[contador]]$RNA@counts[which(rownames(Seurat_list_Haradvala[[contador]]) %in% c("Kymriah", "Yescarta")), ] > 0)

    Sum_Kymriah[[contador]] <- sum(Seurat_list_Haradvala[[contador]]$RNA@counts[which(rownames(Seurat_list_Haradvala[[contador]]) == "Kymriah"), ] > 0)
    Sum_Yescarta[[contador]] <- sum(Seurat_list_Haradvala[[contador]]$RNA@counts[which(rownames(Seurat_list_Haradvala[[contador]]) == "Yescarta"), ] > 0)

    print(names(Seurat_list_Haradvala)[contador])
    Sys.sleep(0.5)
    print(Seurat_list_Haradvala[[contador]])
    print(paste0("Total CAR sum: ", Sum_CAR[[contador]]))
    print(paste0("Total Kymriah: ", Sum_Kymriah[[contador]]))
    print(paste0("Total Yescarta: ", Sum_Yescarta[[contador]]))
    Sys.sleep(0.5)
}

### Add more metadata: novelty score and mito ratio
load("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Gene_Markers_Info/cycle.rda")
for (contador in seq_along(Seurat_list_Haradvala)) {
    Seurat_list_Haradvala[[contador]]$log10GenesPerUMI <- log10(Seurat_list_Haradvala[[contador]]$nFeature_RNA) / log10(Seurat_list_Haradvala[[contador]]$nCount_RNA)

    Seurat_list_Haradvala[[contador]]$mitoRatio <- (PercentageFeatureSet(object = Seurat_list_Haradvala[[contador]], pattern = "^MT-")) / 100

    Seurat_list_Haradvala[[contador]]$log10GenesPerUMI <- log10(Seurat_list_Haradvala[[contador]]$nFeature_RNA) / log10(Seurat_list_Haradvala[[contador]]$nCount_RNA)

}

### Pre-filter
minCov <- 1000 # if a sample has a good coverage, then it doesn't set a lower thresold for nCount.
countLOW <- list()
countHIGH <- list()
featureLOW <- list()
featureHIGH <- list()

if (exists("Pre_filter")) {
    print("Prefiltering DONE")
} else {
    for (contador in seq_along(Seurat_list_Haradvala)) {
        if (min(Seurat_list_Haradvala[[contador]]$nCount_RNA) >= minCov) {
            countLOW[[contador]] <- min(Seurat_list_Haradvala[[contador]]$nCount_RNA)
        } else {
            countLOW[[contador]] <- quantile(Seurat_list_Haradvala[[contador]]$nCount_RNA, prob = 0.01)
        }
        countHIGH[[contador]] <- quantile(Seurat_list_Haradvala[[contador]]$nCount_RNA, prob = 0.99)
        featureLOW[[contador]] <- quantile(Seurat_list_Haradvala[[contador]]$nFeature_RNA, prob = 0.01)
        featureHIGH[[contador]] <- quantile(Seurat_list_Haradvala[[contador]]$nFeature_RNA, prob = 0.99)

        ## subset
        Seurat_list_Haradvala[[contador]] <- subset(Seurat_list_Haradvala[[contador]], subset = (nFeature_RNA > featureLOW[[contador]]) &
            (nFeature_RNA < featureHIGH[[contador]]) & (nCount_RNA > countLOW[[contador]]) & (nCount_RNA < countHIGH[[contador]]) & (mitoRatio < 0.20))
    }
    Pre_filter <- "DONE"
}
rm(minCov)

### Group all Products in a single metadata table
merged_metadata_df <- list()
for (cont in seq_along(Seurat_list_Haradvala)) {
    merged_metadata_df[[cont]] <- Seurat_list_Haradvala[[cont]]@meta.data
}
rm(cont)
merged_metadata_df <- bind_rows(merged_metadata_df)

Patient_list <- names(Seurat_list_Haradvala)

### Exploratory analysis
# Values
# Order: Har_Pat1_IP, Har_Pat2_D7, Har_Pat2_IP, Har_Pat3_IP, Har_Pat4_IP, Har_Pat4_D7, Har_Pat5_IP, Har_Pat6_D7, Har_Pat7_D7, Har_Pat7_IP,
#        Har_Pat8_IP, Har_Pat8_D7, Har_Pat9_IP, Har_Pat9_D7, Har_Pat10_IP, Har_Pat10_D7, Har_Pat11_D7, Har_Pat11_IP, Har_Pat12_D7, Har_Pat12_IP,
#        Har_Pat12_D14, Har_Pat13_IP, Har_Pat13_D7, Har_Pat14_D7, Har_Pat14_D14, Har_Pat14_IP, Har_Pat15_D7, Har_Pat15_IP, Har_Pat16_D7, Har_Pat16_IP,
#        Har_Pat17_IP, Har_Pat17_D7, Har_Pat18_D7, Har_Pat18_IP, Har_Pat19_D7, Har_Pat19_IP, Har_Pat20_D7, Har_Pat20_IP, Har_Pat20_D14, Har_Pat21_D7,
#        Har_Pat21_D14, Har_Pat21_IP, Har_Pat22_D7, Har_Pat22_IP, Har_Pat23_D7, Har_Pat23_IP, Har_Pat24_IP, Har_Pat24_D7, Har_Pat25_IP, Har_Pat25_D7,
#        Har_Pat26_IP, Har_Pat26_D7, Har_Pat27_D7, Har_Pat27_IP, Har_Pat28_IP, Har_Pat28_D7, Har_Pat29_IP, Har_Pat29_IP_retreat, Har_Pat29_D7, Har_Pat29_D7_retreat,
#        Har_Pat30_D7, Har_Pat30_IP, Har_Pat31_D7, Har_Pat31_IP, Har_Pat32_IP, Har_Pat32_D7

Min_feat <- c(
    500, 1000, 1000, 1000, 1000, 800, 1000, 1000, 1250, 1000,
    900, 1500, 1000, 1000, 650, 1600, 1500, 1000, 800, 2000,
    800, 1000, 1000, 700, 500, 900, 1000, 500, 1000, 500,
    500, 1500, 1500, 500, 1300, 1600, 750, 1200, 700, 600,
    800, 500, 1100, 800, 1500, 700, 600, 1000, 2000, 1700,
    700, 1500, 1500, 1500, 650, 1600, 1500, 1200, 1500, 1000,
    1000, 1500, 1400, 500, 1500, 1000
)
Min_counts <- c(
    700, 2000, 1700, 2500, 2200, 1000, 1500, 1800, 3000, 3000,
    1500, 3000, 2500, 3000, 1500, 3000, 3000, 2300, 1500, 5000,
    1800, 2000, 2000, 1000, 1000, 1500, 2100, 800, 2700, 750,
    800, 3000, 3000, 900, 3600, 3000, 1000, 3000, 1500, 1000,
    1000, 750, 3000, 1200, 3000, 1300, 1000, 2800, 3000, 4000,
    1000, 3000, 3000, 3000, 1000, 3000, 3000, 3000, 3000, 3500,
    3000, 3000, 3000, 900, 3000, 1600
)
Max_counts <- c(
    13250, 25000, 32000, 17000, 27000, 10000, 28000, 15000, 35000, 27500,
    22000, 27000, 22000, 15000, 26000, 32000, 47000, 28000, 10000, 46000,
    10000, 37000, 25000, 10000, 4000, 19000, 22500, 8500, 18000, 9500,
    5000, 23000, 10000, 10000, 17000, 52000, 10000, 30000, 8000, 7000,
    5000, 6000, 25000, 11000, 15000, 10000, 10000, 11000, 48000, 25000,
    7000, 25000, 12000, 32000, 8000, 28000, 32000, 27000, 30000, 24000,
    20000, 40000, 26500, 7500, 35000, 17500
)
Max_mito_ratio <- rep(0.1, 66)

Seurat_list_Haradvala[[32]]$nCount_RNA %>% max()

for (i in seq_along(Seurat_list_Haradvala)) {
    a <- max(Seurat_list_Haradvala[[i]]@meta.data$nCount_RNA)
    print(paste0(i, ": ", a))
}

colores <- hue_pal()(length(Seurat_list_Haradvala))

# Visualize the number of cell counts per sample
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Exploratory_analysis/Pre_QC_Ncells.pdf")
merged_metadata_df %>%
    ggplot(aes(x = Product, fill = Product)) +
    geom_bar() +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 3) +
    theme_classic() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none" # Remove the legend
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    ggtitle("Number of Cells pre-QC")
dev.off()

# Transformation into single cell experiment + Knee plots representation
counter <- 0
sce_Seurat_list_Haradvala <- list()
metadata_bcrank <- list()
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Exploratory_analysis/Pre_QC_Knee_plots.pdf")
for (i in Patient_list) {
    counter <- counter + 1
    sce_Seurat_list_Haradvala[[i]] <- as.SingleCellExperiment(Seurat_list_Haradvala[[i]])
    bcrank <- DropletUtils::barcodeRanks(counts(sce_Seurat_list_Haradvala[[i]]))

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
for (i in seq_along(Seurat_list_Haradvala)) {
    pt_h[[i]] <- ggplot(data = Seurat_list_Haradvala[[i]]@meta.data, aes(x = nCount_RNA)) +
        geom_density(alpha = 0.2, fill = colores[i]) +
        scale_x_log10() +
        theme_classic() +
        ylab("Cell density") +

        # Adjust thresholds per dataset
        geom_vline(xintercept = Min_counts[i], col = colores[i], alpha = 0.8, linetype = "dotted") +
        geom_vline(xintercept = Max_counts[i], col = colores[i], alpha = 0.8, linetype = "dotted") +
        labs(
            title = "Number UMIs/transcripts per cell pre-QC:",
            subtitle = paste0(unique(Seurat_list_Haradvala[[i]]@meta.data$Product))
        )
}

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Exploratory_analysis/Pre_QC_UMIs_per_cell.pdf")
pt_h
dev.off()
rm(i, pt_h)

# Visualize the distribution of genes detected per cell
pt_h <- list()
for (i in seq_along(Seurat_list_Haradvala)) {
    pt_h[[i]] <- ggplot(data = Seurat_list_Haradvala[[i]]@meta.data, aes(x = nFeature_RNA)) +
        geom_density(alpha = 0.2, fill = colores[i]) +
        theme_classic() +
        scale_x_log10() +

        # Adjust thresholds per dataset
        geom_vline(xintercept = Min_feat[i], col = colores[i], alpha = 0.6, linetype = "dotted") +
        labs(
            title = "Number genes per cell pre-QC:",
            subtitle = paste0(unique(Seurat_list_Haradvala[[i]]@meta.data$Product))
        )
}

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Exploratory_analysis/Pre_QC_Genes_per_cell.pdf")
pt_h
dev.off()
rm(i, pt_h)

# Visualize the distribution of mitochondrial gene expression detected per cell
pt_h <- list()
for (i in seq_along(Seurat_list_Haradvala)) {
    pt_h[[i]] <- ggplot(data = Seurat_list_Haradvala[[i]]@meta.data, aes(x = mitoRatio)) +
        geom_density(alpha = 0.2, fill = colores[i]) +
        scale_x_log10() +
        theme_classic() +

        # Adjust thresholds per dataset
        geom_vline(xintercept = Max_mito_ratio[i], col = colores[i], alpha = 0.6, linetype = "dotted") +
        labs(
            title = "Mitocondrial ratio per cell:",
            subtitle = paste0(unique(Seurat_list_Haradvala[[i]]@meta.data$Product))
        )
}

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Exploratory_analysis/Pre_QC_Mito_ratio_per_cell.pdf")
pt_h
dev.off()
rm(i, pt_h)

# Visualize the distribution of complexity per cell
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Exploratory_analysis/Pre_QC_Complexity_per_cell.pdf")
merged_metadata_df %>%
    ggplot(aes(color = Product, x = log10GenesPerUMI, fill = Product)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    theme(
        legend.position = "none" # Remove the legend
    ) +
    ggtitle("Complexity per cell")
dev.off()

# Violin plots
p_Vln <- list()
for (cont in seq_along(Seurat_list_Haradvala)) {
    p_Vln[[cont]] <- VlnPlot(Seurat_list_Haradvala[[cont]], features = c("nCount_RNA", "nFeature_RNA", "mitoRatio"), col = colores[cont], pt.size = 0, group.by = "Product")
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Exploratory_analysis/Pre_QC_Vln_plots.pdf")
p_Vln
dev.off()

rm(cont, p_Vln)

# nGenes vs nUMI graph (Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs)
p1 <- list()
for (i in seq_along(Seurat_list_Haradvala)) {
    p1[[i]] <- ggplot(data = Seurat_list_Haradvala[[i]]@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = mitoRatio)) +
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
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Exploratory_analysis/Pre_QC_nGenes_vs_nUMI.pdf")
p1
dev.off()

rm(i, p1)

# mitoRatio vs nUMI graph
p1 <- list()
for (i in seq_along(Seurat_list_Haradvala)) {
    p1[[i]] <- ggplot(data = Seurat_list_Haradvala[[i]]@meta.data, aes(x = nCount_RNA, y = mitoRatio)) +
        geom_point(size = 0.7) +
        theme_classic() +
        ggtitle(paste0(Patient_list[i]))
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Exploratory_analysis/Pre_QC_mitoRatio_vs_nUMI.pdf")
p1
dev.off()

rm(i, p1)

### Subset samples by determined quality control thresholds
for (i in seq_along(Seurat_list_Haradvala)) {
    Seurat_list_Haradvala[[i]] <- subset(Seurat_list_Haradvala[[i]],
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

Seurat_doublets_finder <- Seurat_list_Haradvala

sweep_res_list_Haradvala <- list()
sweep_stats_list_Haradvala <- list()
bcmvn_list_Haradvala <- list()
pK <- list()
BCmetric <- list()
pK_choose <- list()
nExp <- list()

per_doubl <- c(
    0.096, 0.024, 0.02, 0.068, 0.008, 0.068, 0.044, 0.008, 0.024, 0.008,
    0.012, 0.008, 0.04, 0.008, 0.016, 0.004, 0.004, 0.004, 0.056, 0.004,
    0.072, 0.004, 0.004, 0.072, 0.16, 0.048, 0.008, 0.008, 0.012, 0.008,
    0.004, 0.008, 0.008, 0.004, 0.004, 0.004, 0.044, 0.024, 0.072, 0.104,
    0.088, 0.088, 0.04, 0.032, 0.016, 0.032, 0.036, 0.036, 0.04, 0.004,
    0.06, 0.012, 0.044, 0.08, 0.02, 0.016, 0.032, 0.032, 0.012, 0.012,
    0.008, 0.016, 0.036, 0.024, 0.072, 0.04
) # Calculated looking at 10X table

# pK Identification (no ground-truth)
for (i in seq_along(Seurat_doublets_finder)) {
    Seurat_doublets_finder[[i]] <- NormalizeData(Seurat_doublets_finder[[i]])
    Seurat_doublets_finder[[i]] <- FindVariableFeatures(Seurat_doublets_finder[[i]], selection.method = "vst", nfeatures = 2000)
    Seurat_doublets_finder[[i]] <- ScaleData(Seurat_doublets_finder[[i]])
    Seurat_doublets_finder[[i]] <- RunPCA(Seurat_doublets_finder[[i]], npcs = 20)
    Seurat_doublets_finder[[i]] <- RunUMAP(Seurat_doublets_finder[[i]], dims = 1:10)
    nExp[[i]] <- round(ncol(Seurat_doublets_finder[[i]]) * per_doubl[i])

    sweep_res_list_Haradvala[[i]] <- paramSweep_v3(Seurat_doublets_finder[[i]], PCs = 1:10, sct = FALSE, num.cores = 12)
    sweep_stats_list_Haradvala[[i]] <- summarizeSweep(sweep_res_list_Haradvala[[i]], GT = FALSE)
    bcmvn_list_Haradvala[[i]] <- find.pK(sweep_stats_list_Haradvala[[i]])
    dev.off()

    pK[[i]] <- as.numeric(as.character(bcmvn_list_Haradvala[[i]]$pK))
    BCmetric[[i]] <- bcmvn_list_Haradvala[[i]]$BCmetric
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

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Exploratory_analysis/Doublet_Finder_pK_calculation.pdf")
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

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Exploratory_analysis/Doublet_Finder_Classification.pdf")
h2
dev.off()
rm(i, h2)

# Remove all predicted doublets from our data.
for (i in seq_along(Seurat_doublets_finder)) {
    print(Seurat_doublets_finder[[i]])

    Sys.sleep(0.5)

    Singlets <- rownames(Seurat_doublets_finder[[i]]@meta.data)[Seurat_doublets_finder[[i]]@meta.data[, DF_name[[i]]] == "Singlet"]
    print(length(Singlets))

    Sys.sleep(0.5)

    Seurat_list_Haradvala[[i]] <- subset(Seurat_list_Haradvala[[i]], cells = Singlets)
    print(Seurat_list_Haradvala[[i]])

    Sys.sleep(0.5)

    print(paste0(Patient_list[i], " DONE"))
    cat(paste0("----------------------------------------------------------------------------\n"))
}

### Score cells for cell cycle
Seurat_list_Haradvala_fix <- Seurat_list_Haradvala

for(contador in seq_along(Seurat_list_Haradvala_fix)){
    Seurat_list_Haradvala_fix[[contador]] <- NormalizeData(Seurat_list_Haradvala_fix[[contador]])

    # Score cells for cell cycle
    Seurat_list_Haradvala_fix[[contador]] <- CellCycleScoring(Seurat_list_Haradvala_fix[[contador]],
        g2m.features = g2m_genes,
        s.features = s_genes
    )
}

for(contador in seq_along(Seurat_list_Haradvala_fix)){
    Seurat_list_Haradvala[[contador]]@meta.data <- Seurat_list_Haradvala_fix[[contador]]@meta.data
}

rm(Seurat_list_Haradvala_fix)


### Subset by CAR+ cells
names_to_remove <- list()
for (i in seq_along(Seurat_list_Haradvala)) {
    if (any(rownames(Seurat_list_Haradvala[[i]]) == "Kymriah") && !any(rownames(Seurat_list_Haradvala[[i]]) == "Yescarta")) {
        Seurat_list_Haradvala[[i]] <- subset(x = Seurat_list_Haradvala[[i]], subset = Kymriah > 0)
        print("Caso 1")
    } else if (!any(rownames(Seurat_list_Haradvala[[i]]) == "Kymriah") && any(rownames(Seurat_list_Haradvala[[i]]) == "Yescarta")) {
        Seurat_list_Haradvala[[i]] <- subset(x = Seurat_list_Haradvala[[i]], subset = Yescarta > 0)
        print("Caso 2")
    } else if (any(rownames(Seurat_list_Haradvala[[i]]) == "Kymriah") && any(rownames(Seurat_list_Haradvala[[i]]) == "Yescarta")) {
        Seurat_list_Haradvala[[i]] <- subset(x = Seurat_list_Haradvala[[i]], subset = Kymriah > 0 | Yescarta > 0)
        print("Caso 3")
    } else {
        print("Caso 4")
        print(paste0(names(Seurat_list_Haradvala)[i], " es CAR-"))
        names_to_remove <- c(names_to_remove, names(Seurat_list_Haradvala)[i])
        Sys.sleep(2)
    }
}

Seurat_list_Haradvala %>% length()
Seurat_list_Haradvala <- Seurat_list_Haradvala[!(names(Seurat_list_Haradvala) %in% names_to_remove)]
Seurat_list_Haradvala %>% length()

### Re-doing of exploratory analysis graphs after QC
# Group again all Products in a single metadata table
merged_metadata_df <- list()
for (cont in seq_along(Seurat_list_Haradvala)) {
    merged_metadata_df[[cont]] <- Seurat_list_Haradvala[[cont]]@meta.data
}
rm(cont)
merged_metadata_df <- bind_rows(merged_metadata_df)

# Visualize the number of cell counts per sample
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Exploratory_analysis/Post_QC_Ncells.pdf")
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
for (i in seq_along(Seurat_list_Haradvala)) {
    pt_h[[i]] <- ggplot(data = Seurat_list_Haradvala[[i]]@meta.data, aes(x = nCount_RNA)) +
        geom_density(alpha = 0.2, fill = colores[i]) +
        scale_x_log10() +
        theme_classic() +
        ylab("Cell density") +
        labs(
            title = "Number UMIs/transcripts per cell post-QC:",
            subtitle = paste0(unique(Seurat_list_Haradvala[[i]]@meta.data$Product))
        )
}

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Exploratory_analysis/Post_QC_UMIs_per_cell.pdf")
pt_h
dev.off()
rm(i, pt_h)

# Visualize the distribution of genes detected per cell
pt_h <- list()
for (i in seq_along(Seurat_list_Haradvala)) {
    pt_h[[i]] <- ggplot(data = Seurat_list_Haradvala[[i]]@meta.data, aes(x = nFeature_RNA)) +
        geom_density(alpha = 0.2, fill = colores[i]) +
        theme_classic() +
        scale_x_log10() +
        labs(
            title = "Number genes per cell post-QC:",
            subtitle = paste0(unique(Seurat_list_Haradvala[[i]]@meta.data$Product))
        )
}

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Exploratory_analysis/Post_QC_Genes_per_cell.pdf")
pt_h
dev.off()
rm(i, pt_h)

# Visualize the distribution of mitochondrial gene expression detected per cell
pt_h <- list()
for (i in seq_along(Seurat_list_Haradvala)) {
    pt_h[[i]] <- ggplot(data = Seurat_list_Haradvala[[i]]@meta.data, aes(x = mitoRatio)) +
        geom_density(alpha = 0.2, fill = colores[i]) +
        scale_x_log10() +
        theme_classic() +

        # Adjust thresholds per dataset
        labs(
            title = "Mitocondrial ratio per cell:",
            subtitle = paste0(unique(Seurat_list_Haradvala[[i]]@meta.data$Product))
        )
}

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Exploratory_analysis/Post_QC_Mito_ratio_per_cell.pdf")
pt_h
dev.off()
rm(i, pt_h)

# Violin plots
p_Vln <- list()
for (cont in seq_along(Seurat_list_Haradvala)) {
    p_Vln[[cont]] <- VlnPlot(Seurat_list_Haradvala[[cont]], features = c("nCount_RNA", "nFeature_RNA", "mitoRatio"), col = colores[cont], pt.size = 0, group.by = "Product")
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Exploratory_analysis/Post_QC_Vln_plots.pdf")
p_Vln
dev.off()

rm(cont, p_Vln)

# nGenes vs nUMI graph (Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs)
p1 <- list()
for (i in seq_along(Seurat_list_Haradvala)) {
    p1[[i]] <- ggplot(data = Seurat_list_Haradvala[[i]]@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = mitoRatio, grou)) +
        geom_point(size = 0.7) +
        scale_color_gradientn(colors = c("black", "blue", "green2", "red", "yellow")) +
        stat_smooth(method = lm) +
        scale_x_log10() +
        scale_y_log10() +
        theme_classic() +
        ggtitle(paste0(Patient_list[i]))
}
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Exploratory_analysis/Post_QC_nGenes_vs_nUMI.pdf")
p1
dev.off()

rm(i, p1)

setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/RDS")
saveRDS(Seurat_list_Haradvala, file = "PostQC_CellRanger_Haradvala_RDS.rds")

##### Normalization #####

# rm(list = ls())
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/RDS")
Seurat_list_Haradvala <- readRDS("PostQC_CellRanger_Haradvala_RDS.rds")

for (i in seq_along(Seurat_list_Haradvala)) {
    Seurat_list_Haradvala[[i]] <- NormalizeData(Seurat_list_Haradvala[[i]])

    Seurat_list_Haradvala[[i]] <- FindVariableFeatures(Seurat_list_Haradvala[[i]],
        selection.method = "vst",
        nfeatures = 2000
    )

    Seurat_list_Haradvala[[i]] <- ScaleData(Seurat_list_Haradvala[[i]])
}
rm(i)

saveRDS(Seurat_list_Haradvala, file = "Normalized_CellRanger_Haradvala_RDS.rds")


##### Merge ---- Worst case scenario #####
setwd("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/RDS")
Seurat_list_Haradvala <- readRDS("Normalized_CellRanger_Haradvala_RDS.rds")

Seurat_list_Haradvala_2 <- merge(x = Seurat_list_Haradvala[[1]], y = do.call(c, Seurat_list_Haradvala[2:65]))

Seurat_list_Haradvala <- Seurat_list_Haradvala_2
rm(Seurat_list_Haradvala_2)

Seurat_list_Haradvala <- FindVariableFeatures(Seurat_list_Haradvala,
    selection.method = "vst",
    nfeatures = 2000
)

Seurat_list_Haradvala <- ScaleData(Seurat_list_Haradvala)

# Run PCA
Seurat_list_Haradvala <- RunPCA(object = Seurat_list_Haradvala, reduction.name = "pca_wo_integ")

# ElbowPlot(Seurat_list_Haradvala, ndims = 50, reduction = "pca_wo_integ")
# dev.off()

# Run UMAP
Seurat_list_Haradvala <- RunUMAP(Seurat_list_Haradvala,
    dims = 1:30,
    reduction = "pca_wo_integ",
    reduction.name = "umap_wo_integ"
)

# Plot UMAP
pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Worst_case_scenario/WO_integ_Seurat.pdf")
DimPlot(Seurat_list_Haradvala,
    group.by = "Product",
    reduction = "umap_wo_integ"
) + NoLegend()
dev.off()

# Quality control metrics
p <- FeaturePlot(Seurat_list_Haradvala,
    reduction = "umap_wo_integ",
    features = c("nCount_RNA", "nFeature_RNA", "log10GenesPerUMI", "mitoRatio"),
    order = TRUE
)

pdf("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Haradvala_et_al/Plots/Worst_case_scenario/QC_metrics_WO_integ_Seurat.pdf")
p + plot_annotation(title = paste0(Seurat_list_Haradvala@meta.data$orig.ident %>% unique()), theme = theme(plot.title = element_text(size = 16)))
dev.off()

rm(p)

##### END OF SCRIPT #####