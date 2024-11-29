
#' Helper Function 1 - get organ columns
#'
#' This function performs a specific task. Gets the organ of interest columns
#' from the mouse or human gene expression data.
#' https://apps.kaessmannlab.org/evodevoapp/
#' mouse IDs are different. add script to match to human orthologs.
#'
#' @param organism A string of either 'mouse' or 'human'.
#' @param organ A string of the organ of interest.
#' @return A data frame with the columns of the organ of interest of either mouse or human.
#' @export
get_organ <- function(organism, organ) {
  # Validate input
  if (!organism %in% c("mouse", "human")) {
    stop("Unsupported organism: choose 'mouse' or 'human'")
  }

  # Select and process the data
  organ <- as.character(organ)
  gene_data <- get(organism) %>%  # Dynamically access the dataset (mouse or human)
    dplyr::select('Names', starts_with(organ)) %>%
    dplyr::rename(ensembl_gene_id = Names)

  return(gene_data)
}



#' Main function - get gene expression of the organ and organism of interest
#'
#' This is the main functions. Gets the gene expression RPKM of the organ and
#' organism of interest. Either mouse or human.
#'
#' @param protein_coding_genes A df with at least one column of hgnc_id and ensembl_gene_id
#' @param organism A string of either 'mouse' or 'human'.
#' @param organ A vector of strings of the organs of interest.
#' @return A data frame with the columns of the organ of interest of either mouse or human of gene expression.
#' @export
get_gene_expression <- function(protein_coding_genes, organism, vector_of_organs) {

  results <- data.frame()

  # Loop through the vector of organs and merge results
  for (string in vector_of_organs) {
    result <- get_organ(organism, string)
    if (nrow(results) == 0) {
      results <- result
    } else {
      results <- results %>%
        full_join(result, by = "ensembl_gene_id")
    }
  }

  # Merge the combined organ data with protein_coding_genes
  organ_hgnc <- protein_coding_genes %>%
    left_join(results, by = 'ensembl_gene_id') %>%
    distinct()

  return(organ_hgnc)
}
