library(GSEABase)

# Function to extract gene sets from a GeneSetCollection based on pathway prefixes
extract_gene_sets <- function(gsc, pathway_prefixes) {
  if (!inherits(gsc, "GeneSetCollection")) {
    stop("The input must be a GeneSetCollection object.")
  }
  
  # Extract pathway IDs and their associated gene sets
  gene_sets <- setNames(
    lapply(gsc, function(gs) geneIds(gs)),  # Extract gene symbols as a list
    sapply(gsc, setName)  # Use pathway IDs as names
  )
  
  # Filter pathways whose IDs start with any of the given prefixes
  selected_gene_sets <- gene_sets[names(gene_sets) %in% grep(paste0("^", paste(pathway_prefixes, collapse = "|")), names(gene_sets), value = TRUE)]
  
  return(selected_gene_sets)
}