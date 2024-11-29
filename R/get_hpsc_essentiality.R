
#' Main function - get hPSC essentiality.
#'
#' This is the main function. Gets the hPSC essentiality scores from RData file.
#' https://www.cell.com/cell-reports/fulltext/S2211-1247(19)30212-8
#'
#' @param protein_coding_genes A data frame with HGNC ids.
#' @return Data frame of hPSC essentiality scores.
#' @export
get_hpsc_essentiality <- function(protein_coding_genes) {

  #load('./data/cell_essentiality.rdata')
  hpsc_essentiality <- left_join(protein_coding_genes, cell_essentiality, by = "hgnc_id") %>%
    dplyr::select(hgnc_id, h1_mef_BF, h1_laminin_BF)

  return(hpsc_essentiality)

}
