
#' Main function - get Exomiser hiphive NDD scores
#'
#' This is a main function. Read the RData file in the data directory. This
#' function only works for NDD genes.
#'
#' @param protein_coding_genes A df of entrez_id and hgnc_id.
#' @return A df with hgnc_id and the NDD hiphive scores.
#' @export
get_exomiser_hiphive_ndd <- function(protein_coding_genes) {

  #load('./data/NDD_HiPhive_Prioritiser_genes.RData')

  hiphive1 <- hiphive %>%
    dplyr::filter(ENTREZ_GENE_ID != -1)

  hgnc_entrez <- protein_coding_genes %>%
    dplyr::select(hgnc_id, entrez_id)

  hgnc_hiphive <- hgnc_entrez %>%
    left_join(hiphive1,
              by = join_by(entrez_id == ENTREZ_GENE_ID))

  hgnc_hiphive1 <- hgnc_hiphive %>%
    dplyr::group_by(hgnc_id) %>%
    dplyr::filter(EXOMISER_GENE_PHENO_SCORE == max(EXOMISER_GENE_PHENO_SCORE)) %>%
    dplyr::select(hgnc_id, EXOMISER_GENE_PHENO_SCORE) %>%
  distinct()

  hgnc_hiphive2 <- protein_coding_genes %>%
    left_join(hgnc_hiphive1) %>%
    dplyr::select(-entrez_id)

  return(hgnc_hiphive2)
}




