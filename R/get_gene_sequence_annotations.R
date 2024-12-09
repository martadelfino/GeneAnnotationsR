
#' Helper function 1 - access biomart structure page
#'
#' This is a helper function. Accesses the structure feature page to get various length annotations.
#'
#' @return A df with ensemble_gene_id, and the various length annotations
#' @export
access_structure_page <- function() {

  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://www.ensembl.org")
  results1 <- biomaRt::getBM(attributes = c('ensembl_gene_id',
                                            'start_position',
                                            'end_position',
                                            'transcript_length',
                                            'cds_length',
                                            '5_utr_start',
                                            '5_utr_end',
                                            '3_utr_start',
                                            '3_utr_end'),
                             filters = 'transcript_is_canonical', # transcript_is_canonical
                             values = TRUE,
                             mart = ensembl)

  return(results1)

}


#' Helper function 2 - access biomart feature page
#'
#' This is a helper function. Accesses the biomart feature page to get %GC.
#'
#' @return A df with ensemble_gene_id, and the %GC
#' @export
access_feature_page <- function() {

  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
  results2 <- biomaRt::getBM(attributes = c('ensembl_gene_id',
                                            'percentage_gene_gc_content',
                                            'version'),
                             filters = 'transcript_is_canonical',
                             values = TRUE,
                             mart = ensembl)

  return(results2)
}


#' Main function - get gene sequence annotations from biomart
#'
#' This is the main function. Gets gene sequence annotations from biomart
#'
#' @param protein_coding_genes A df with at least hgnc_id and ensembl_gene_id
#' @return A df with hgnc_id, ensemble_gene_id, and the gene sequence annotations
#' Specifically: gene length, transcript length, cds length, %GC, and version of the ensembl entry
#' @export
get_gene_sequence_annotations <- function(protein_coding_genes, keep_version = FALSE) {

  biomart1 <- access_structure_page()
  biomart2 <- access_feature_page()

  results <- biomart1 %>%
    left_join(biomart2, by = 'ensembl_gene_id')

  results_mut <- results %>%
    dplyr::mutate(gene_length = end_position - start_position + 1) %>%
    dplyr::select(-end_position, -start_position) %>%
    dplyr::mutate(`3_utr_length` = `3_utr_end` - `3_utr_start` + 1) %>%
    dplyr::select(-`3_utr_end`, -`3_utr_start`) %>%
    dplyr::mutate(`5_utr_length` = `5_utr_end` - `5_utr_start` + 1) %>%
    dplyr::select(-`5_utr_end`, -`5_utr_start`) %>%
    dplyr::mutate(`3_utr_length` = replace_na(`3_utr_length`, 0),
                  `5_utr_length` = replace_na(`5_utr_length`, 0))

  results_summarised <- results_mut %>%
    dplyr::group_by(ensembl_gene_id) %>%
    dplyr::mutate(`3_utr_length` = sum(`3_utr_length`, na.rm = TRUE),
                  `5_utr_length` = sum(`5_utr_length`, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct()

  # Conditionally keep the 'version' column
  if (!keep_version) {
    results_summarised <- results_summarised %>%
      dplyr::select(-version)
  }

  results_hgnc <- protein_coding_genes %>%
    left_join(results_summarised, by = 'ensembl_gene_id')

  return(results_hgnc)
}
