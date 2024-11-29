
#' Get all protein coding genes data
#'
#' This function gets the df from EBI of all protein coding genes.
#'
#' @return A dataframe with protein coding gene data from EBI database.
#' @export
get_protein_coding_genes <- function() {

  # Use fread for faster and memory-efficient reading
  protein_coding_genes <- data.table::fread(
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/locus_types/gene_with_protein_product.txt",
    sep = "\t",     # Define tab-separated format
    header = TRUE,  # File contains headers
    select = NULL   # You can specify column names to load, if needed
  )

  # Convert to data.frame only if necessary
  return(as.data.frame(protein_coding_genes))
}

