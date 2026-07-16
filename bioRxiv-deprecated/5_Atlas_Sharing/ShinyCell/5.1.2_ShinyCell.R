###############################################################################
###############################################################################

# Program: 5.1.2_ShinyCell.R
# Author: Sergio Cámara Peña
# Date: 01/10/2025
# Version: FINAL

###############################################################################
###############################################################################

library(ShinyCell)

setwd("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Shiny_app")

inpFile <- "Python_scVI_adata_big_V4_state4_Normalized_for_shiny.h5ad"

scConf <- createConfig(inpFile)

citation = "Reference pending - manuscript in preparation."

makeShinyApp(
    inpFile, scConf,
    gene.mapping = TRUE,
    shiny.title = "Functional CAR+ T Cell Atlas",
    shiny.dir = "shinyAtlas/",
    shiny.footnotes = citation,
    default.gene1 = "CD4",
    default.gene2 = "CD8A",
    default.multigene = c(
        "CD3D", "CD4", "CD8A", "NKG7", "GNLY",
        "FOXP3", "IL2RA", "CD69", "CD27", "TCF7",
        "CCR7", "SELL"
    )
)

##### END OF SCRIPT #####