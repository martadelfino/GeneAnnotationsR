
#' Helper Function 1 - download protein expression
#'
#'
#'
#' @return A dataframe with protein expression from protein atlas
#' @export
download_normal_prot_expr <- function() {
  # Download and unzip the file, then read the TSV directly
  url <- "https://www.proteinatlas.org/download/tsv/normal_ihc_data.tsv.zip"
  temp <- tempfile()

  # Download the zip file
  download.file(url, temp)

  # Unzip and read the TSV file into a data frame
  df <- read_tsv(unz(temp, "normal_ihc_data.tsv"))

  return(df)
}


#' Main Function - get protein expression
#'
#'
#'
#' @return A dataframe with protein expression mapped to hgnc_id
#' @export
get_protein_expression <- function(protein_coding_genes) {

  protein_expression <- download_normal_prot_expr()

  protein_expression <- protein_expression %>%
    dplyr::filter(Level != "Not detected" & Level != 'Not representative') %>%
    dplyr::select(-'Gene name') %>%
    dplyr::rename('ensembl_gene_id' = 'Gene')

  hgnc_protein_expression <- protein_coding_genes %>%
    dplyr::left_join(protein_expression, by = 'ensembl_gene_id')

  return(hgnc_protein_expression)
}

