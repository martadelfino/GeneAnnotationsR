
#' Main function - get Exomiser hiphive NDD scores
#'
#' This is a main function. Read the RData file in the data directory. This
#' function only works for NDD genes.
#'
#' @param protein_coding_genes A df of entrez_id and hgnc_id
#' @return A df with hgnc_id and the NDD hiphive scores.
#' @export
get_exomiser_pheno_score <- function(protein_coding_genes) {

  #load('./data/NDD_Exomiser_Pheno_Score.RData')

  hiphive1 <- exomiser_pheno_score %>%
    dplyr::filter(ENTREZ_GENE_ID != -1)

  hgnc_entrez <- protein_coding_genes %>%
    dplyr::select(hgnc_id, entrez_id)

  hgnc_hiphive <- hgnc_entrez %>%
    left_join(hiphive1,
              by = join_by(entrez_id == ENTREZ_GENE_ID))

  hgnc_hiphive1 <- hgnc_hiphive %>%
    dplyr::group_by(hgnc_id) %>%
    dplyr::summarize(
      max_EXOMISER_GENE_PHENO_SCORE = max(EXOMISER_GENE_PHENO_SCORE, na.rm = TRUE),
      max_HUMAN_PHENO_SCORE = max(HUMAN_PHENO_SCORE, na.rm = TRUE),
      max_MOUSE_PHENO_SCORE = max(MOUSE_PHENO_SCORE, na.rm = TRUE),
      max_FISH_PHENO_SCORE = max(FISH_PHENO_SCORE, na.rm = TRUE),
      max_MOUSE_FISH_PHENO_SCORE = max(MOUSE_PHENO_SCORE, FISH_PHENO_SCORE, na.rm = TRUE)) %>%
    distinct()

  hgnc_hiphive1[hgnc_hiphive1 == -Inf] <- NA

  hgnc_hiphive2 <- protein_coding_genes %>%
    left_join(hgnc_hiphive1) %>%
    dplyr::select(-entrez_id)

  return(hgnc_hiphive2)
}



