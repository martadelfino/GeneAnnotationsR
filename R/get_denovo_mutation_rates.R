
#' Helper Function 1 - download mutation supplementary file xls
#'
#'
#'
#' @return A temp xls file
#' @export
download_mutation_supp_file <- function() {
  # Create a temporary file name with the .xls extension
  temp_file <- tempfile(fileext = ".xls")

  # Download the file from the URL
  download.file(
    url = "https://pmc.ncbi.nlm.nih.gov/articles/instance/4222185/bin/NIHMS612569-supplement-3.xls",
    destfile = temp_file,
    mode = "wb"
  )

  return(temp_file)
}

#' Helper Function 2 - read mutation supplementary file xls
#'
#'
#'
#' @return A df from the xls mutations file
#' @export
read_sheet <- function(temp_file, sheet_number) {
  # Read the specified sheet from the downloaded file
  df <- read_excel(temp_file, sheet = sheet_number)

  # Convert to a regular data frame if desired
  df <- as.data.frame(df)

  return(df)
}


#' Helper Function 3 - get sheet 2 in mutation supplementary file xls
#'
#'
#'
#' @return A df from the xls mutations file from sheet 2 only
#' @export
access_denovo_mutation_rates_sheet <- function() {
  # Download the mutation rate supplement file
  temp_file <- download_mutation_supp_file()

  # Read the second sheet
  df2 <- read_sheet(temp_file, 2)

  # Clean up the temporary file
  unlink(temp_file)

  return(df2)
}


#' Main Function - get denovo mutation rates
#'
#'
#'
#' @return A dataframe with protein expression mapped to hgnc_id
#' @export
get_denovo_mutation_rates <- function(protein_coding_genes) {

  pcg_for_symbolcheck <- get_protein_coding_genes()

  pcg <- protein_coding_genes %>%
    dplyr::select(hgnc_id, symbol) %>%
    dplyr::rename(gene_symbol = symbol)

  # Get the data from the second sheet
  denovo <- access_denovo_mutation_rates_sheet()
  # pull gene symbols
  symbols <- dplyr::pull(denovo, gene)
  # Check the gene symbols
  denovo_symbols_check <- hgnc.checker(symbols, pcg_for_symbolcheck)
  # Merge the data with the gene symbols
  denovo_with_symbcheck <- denovo_symbols_check %>%
    inner_join(denovo, by = c("gene_symbol" = "gene"))

  # merging denovo symbol checked and pcg
  denovo_pcg <- pcg %>%
    left_join(denovo_with_symbcheck, by = c("hgnc_id" = "hgnc_id")) %>%
    dplyr::select(-`gene_symbol.x`, -`gene_symbol.y`, -`refseq_accession`, -`type`,
                  -`transcript`)

  return(denovo_pcg)
}
