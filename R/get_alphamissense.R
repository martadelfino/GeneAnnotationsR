#' Helper Function 1 - read alphamissense file from zenodo
#'
#' This function performs a specific task. Reads the alphamissense gene hg38
#' scores file from zenodo.
#'
#' @return Dataframe of alphamissense gene hg38 scores.
#' @export
read_alphamissense_gene_hg38 <- function() {
  url <- "https://zenodo.org/records/8208688/files/AlphaMissense_gene_hg38.tsv.gz"
  temp_file <- tempfile()
  download.file(url, temp_file)
  alphamissense <- readr::read_tsv(temp_file, skip = 3, col_names = TRUE)

  alphamissense <- alphamissense %>%
    dplyr::mutate(ensembl_transcript_id = str_remove(transcript_id, "\\..*")) %>%
    dplyr::select(-transcript_id)

  return(alphamissense)
}



#' Helper Function 2 - Maps alphamissense ENST to ENSG and HGNC IDs
#'
#' This function performs a specific task. Uses biomart to map alphamissense
#' ENST ids to ENSG and HGNC IDs.
#'
#' @param alphamissense A df that is output of read_alphamissense_gene_hg38().
#' It is a df with two columns 1. ENST and 2. mean_am_pathogenicity.
#' @return Dataframe of alphamissense gene hg38 scores with their mapped ENSG,
#' HGNC IDs.
#' @export
map_alphamissense_ids <- function(alphamissense) {

  enst <- as.vector(alphamissense$ensembl_transcript_id)
  ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attributes <- c('ensembl_transcript_id','ensembl_gene_id', 'hgnc_id')
  mapping <- biomaRt::getBM(attributes = attributes,
                            filters= 'ensembl_transcript_id',
                            values = enst,
                            mart = ensembl)

  alphamissense_mapping <- alphamissense %>%
    left_join(mapping, by = join_by(ensembl_transcript_id))

  return(alphamissense_mapping)
}




#' Main function - get alphamissense gene scores
#'
#' This is the main function. Gets alphamissense gene sequence scores from
#' Zenodo and adds HGNC IDs
#'
#' @param genes A df with at least one column with hgnc_id
#' @return A df with hgnc_id and the max alphamissense score for that gene
#' @export
get_alphamissense <- function(genes) {
  #source('./R/get_protein_coding_genes.R')
  protein_coding_genes <- get_protein_coding_genes() %>%
    dplyr::select(hgnc_id, ensembl_gene_id, refseq_accession, mane_select) %>%
    separate(mane_select, into = c('ensembl_transcript_id', 'ncbi_refseq'),
             sep ="\\|") %>%
    dplyr::mutate(ensembl_transcript_id = str_remove(ensembl_transcript_id, "\\..*"))

  alphamissense <- read_alphamissense_gene_hg38()
  alphamissense_mapping <- map_alphamissense_ids(alphamissense)

  hgnc_alphamissense <- protein_coding_genes %>%
    left_join(alphamissense_mapping, by = 'hgnc_id')

  hgnc_alphamissense_max <- hgnc_alphamissense %>%
    dplyr::select(hgnc_id, mean_am_pathogenicity) %>%
    group_by(hgnc_id) %>%
    slice_max(order_by = mean_am_pathogenicity, n = 1) %>%
    ungroup() %>%
    distinct()

  genes_alphamissense <- genes %>%
    left_join(hgnc_alphamissense_max, by = join_by(hgnc_id))

  return(genes_alphamissense)
}






#source('./R/get_protein_coding_genes.R')

#genes <- get_protein_coding_genes() %>%
 # dplyr::select(hgnc_id)
#genes5 <- head(genes, 5)

#test <- get_alphamissense(genes5)



