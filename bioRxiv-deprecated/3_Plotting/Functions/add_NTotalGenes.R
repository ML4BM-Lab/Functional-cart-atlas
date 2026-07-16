# Function to add NTotalGenes column to res.zenith
add_NTotalGenes <- function(res.zenith, go.gs) {
  # Initialize an empty vector to store NTotalGenes values
  NTotalGenes <- numeric(nrow(res.zenith))
  
  # Loop through each row in res.zenith and retrieve NTotalGenes for each Geneset
  for (i in 1:nrow(res.zenith)) {
    # Get the pathway name from the Geneset column
    target_pathway <- res.zenith$Geneset[i]
    
    # Check if the pathway exists in the GeneSetCollection (go.gs)
    if (target_pathway %in% names(go.gs)) {
      # Extract the gene set for the pathway
      matching_geneset <- go.gs[[target_pathway]]
      
      # Get the total number of genes in the gene set
      NTotalGenes[i] <- length(geneIds(matching_geneset))
    } else {
      # If the pathway isn't found, set NTotalGenes to NA
      NTotalGenes[i] <- NA
    }
  }
  
  # Add the NTotalGenes vector as a new column to res.zenith
  res.zenith$NTotalGenes <- NTotalGenes
  
  # Return the updated res.zenith with the new column
  return(res.zenith)
}
