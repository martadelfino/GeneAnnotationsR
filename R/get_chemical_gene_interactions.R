
#' Helper Function 1 - downlaod chemical gene interactions
#'
#'
#'
#' @return Df of the raw chemical gene interactions data
#' @export
download_chemical_gene_interactions <- function() {

  url <- "https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz"
  temp <- tempfile()

  # Download the gzipped file
  download.file(url, temp)

  # Uncompress and read the TSV file into a data frame
  df <- read.delim(gzfile(temp), skip = 27, header = FALSE, sep = ',')

  df <- df[-2, ]

  colnames(df) <- as.character(df[1, ])
  df <- df[-1, ]

  # Rename the column A to remove the '#'
  colnames(df)[colnames(df) == '# ChemicalName'] <- 'ChemicalName'

  # Remove temporary file and run garbage collection
  #unlink(temp)
  #gc()

  return(df)
}


#' Main Function 1 - gets chemical gene interactions for each gene
#'
#'
#' @param protein_coding_genes Needs to be a df with hgnc_id and entrez_id columns
#' @return Df of the chemical gene interactions data joined by entrez_id
#' @export
get_chemical_gene_interactions <- function(protein_coding_genes) {
  # only using Human at the moment, but could use mouse too
  chemical_gene_interactions <- download_chemical_gene_interactions()

  print(chemical_gene_interactions)

  chemical_gene_interactions <- chemical_gene_interactions %>%
    dplyr::rename('entrez_id' = 'GeneID') %>%
    dplyr::select(-'ChemicalName', -'CasRN', -'GeneSymbol',
                  -'PubMedIDs', -'Organism') %>%
    dplyr::filter(OrganismID == '9606') %>%
    dplyr::mutate(entrez_id = as.character(entrez_id)) %>%
    dplyr::select(-'OrganismID')

  #print(chemical_gene_interactions)

  hgnc_chemical_gene_interactions <- protein_coding_genes %>%
    dplyr::mutate(entrez_id = as.character(entrez_id)) %>%
    left_join(chemical_gene_interactions, by = 'entrez_id')

  # Remove unused data and force garbage collection
 # rm(chemical_gene_interactions)
  #gc()

  return(hgnc_chemical_gene_interactions)
}
