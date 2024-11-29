
#' Check HGNC Gene Symbols
#'
#' This function checks HGNC gene symbols against a database of protein-coding genes
#' and returns matching information.
#'
#' @param gene_symbols A vector of gene symbols to check.
#' @param file_path A file path to the HGNC data file. If the file does not exist, it will be downloaded.
#' @param download_url URL to download the HGNC data file if it is not present.
#' @return A data frame with columns `hgnc_id`, `gene_symbol`, and `type`.
#' @importFrom dplyr select mutate filter bind_rows
#' @importFrom readr read_delim
#' @importFrom tidyr separate_rows
#' @importFrom utils download.file
#' @export
hgnc_checker <- function(gene_symbols,
                         file_path = "gene_with_protein_product.txt",
                         download_url = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/locus_types/gene_with_protein_product.txt") {
  # Ensure the file is available, download if necessary
  if (!file.exists(file_path)) {
    message("File not found. Downloading the file...")
    tryCatch({
      utils::download.file(download_url, file_path, mode = "wb")
      message("File downloaded successfully.")
    }, error = function(e) {
      stop("Error downloading the file: ", e$message)
    })
  }

  # Load the data file
  tryCatch({
    gene_file <- readr::read_delim(file_path, delim = "\t", col_names = TRUE)
  }, error = function(e) {
    stop("Error reading the file: ", e$message)
  })

  # Define helper functions
  check.approved <- function(input_genes, database) {
    database %>%
      dplyr::select(hgnc_id, symbol) %>%
      dplyr::filter(symbol != "", !is.na(symbol), symbol %in% input_genes) %>%
      dplyr::rename(gene_symbol = symbol) %>%
      dplyr::mutate(type = "approved_symbol")
  }

  check.synonyms <- function(input_genes, database) {
    database %>%
      dplyr::select(hgnc_id, alias_symbol) %>%
      tidyr::separate_rows(alias_symbol, sep = "\\|") %>%
      dplyr::mutate(gene_symbol = trimws(alias_symbol)) %>%
      dplyr::filter(gene_symbol %in% input_genes) %>%
      dplyr::mutate(type = "synonym_symbol")
  }

  check.previous <- function(input_genes, database) {
    database %>%
      dplyr::select(hgnc_id, prev_symbol) %>%
      tidyr::separate_rows(prev_symbol, sep = "\\|") %>%
      dplyr::mutate(gene_symbol = trimws(prev_symbol)) %>%
      dplyr::filter(gene_symbol %in% input_genes) %>%
      dplyr::mutate(type = "previous_symbol")
  }

  # Check gene symbols
  genes <- trimws(gene_symbols)
  hgnc_approved <- check.approved(genes, gene_file)
  hgnc_synonyms <- check.synonyms(genes[!genes %in% hgnc_approved$gene_symbol], gene_file)
  hgnc_previous <- check.previous(
    genes[!genes %in% c(hgnc_approved$gene_symbol, hgnc_synonyms$gene_symbol)],
    gene_file
  )
  genes_not_found <- genes[!genes %in% c(hgnc_approved$gene_symbol, hgnc_synonyms$gene_symbol, hgnc_previous$gene_symbol)]
  hgnc_notfound <- data.frame(
    hgnc_id = rep("-", length(genes_not_found)),
    gene_symbol = genes_not_found,
    type = "notfound_proteincoding_symbol"
  )

  # Combine results
  results <- dplyr::bind_rows(hgnc_approved, hgnc_synonyms, hgnc_previous, hgnc_notfound)

  # Handle duplicates
  duplicates_symbol <- results %>%
    dplyr::filter(duplicated(gene_symbol)) %>%
    dplyr::pull(gene_symbol)
  if (length(duplicates_symbol)) {
    results <- results %>%
      dplyr::filter(!gene_symbol %in% duplicates_symbol) %>%
      dplyr::bind_rows(data.frame(hgnc_id = "-", gene_symbol = duplicates_symbol, type = "ambiguous_symbol"))
  }

  duplicates_id <- results %>%
    dplyr::filter(duplicated(hgnc_id), hgnc_id != "-") %>%
    dplyr::pull(hgnc_id)
  if (length(duplicates_id)) {
    results <- results %>%
      dplyr::filter(!hgnc_id %in% duplicates_id) %>%
      dplyr::bind_rows(results %>%
                         dplyr::filter(hgnc_id %in% duplicates_id) %>%
                         dplyr::mutate(hgnc_id = "-", type = "ambiguous_symbol"))
  }

  return(results)
}
