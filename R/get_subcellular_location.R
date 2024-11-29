
#' Helper Function 1 - download subcellular locations from protein atlas
#'
#'
#'
#' @return A dataframe with subcellular locations from protein atlas
#' @export
download_subcell_location_file <- function() {
  # Download and unzip the file, then read the TSV directly
  url <- "https://www.proteinatlas.org/download/tsv/subcellular_location.tsv.zip"
  temp <- tempfile()

  # Download the zip file
  download.file(url, temp)

  # Unzip and read the TSV file into a data frame
  df <- read_tsv(unz(temp, "subcellular_location.tsv"))

  return(df)
}


#' Main Function - get subcellular locations
#'
#' keeping all reliability levels.
#'
#' @return A dataframe with subcellular locations mapped to hgnc_id
#' @export
get_subcellular_location <- function(protein_coding_genes) {

  # Subcellular location file
  subcellular <- download_subcell_location_file()

  subcellular <- subcellular %>%
    dplyr::rename(ensembl_gene_id = Gene)

  hgnc_subcellular <- protein_coding_genes %>%
    dplyr::select(hgnc_id, ensembl_gene_id) %>%
    left_join(subcellular, by = 'ensembl_gene_id') %>%
    dplyr::select(-`Gene name`, -`ensembl_gene_id`, -`Reliability`, -`GO id`)

  return(hgnc_subcellular)
}

