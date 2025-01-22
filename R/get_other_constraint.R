
#' Main function - get Shet, DOMINO, and SCoNes constraint
#'
#' This is the main function. Reads the RData file with the Shet, DOMINO, and
#' SCoNes scores as these are static.
#'
#' @param protein_coding_genes A df with at least one column with hgnc_id.
#' @return A df with hgnc_id and the Shet, DOMINO, and SCoNes scores.
#' @export
get_other_constraint <- function(protein_coding_genes) {

  #load('./data/gene_constraint_metrics.rdata')

  metrics <- gene_constraint_metrics %>%
    dplyr::select(hgnc_id, starts_with('shet'), domino_score, scones_score)

  genes_metrics <- protein_coding_genes %>%
    left_join(metrics, by = 'hgnc_id') %>%
    dplyr::select(-`ensembl_gene_id`, -`symbol`)

  return(genes_metrics)
}
