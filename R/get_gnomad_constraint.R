
#' Main function - get Gnomad constraint scores
#'
#' This is the main function. Gets Gnomad contraints scores from v4.1
#'
#' @param protein_coding_genes A df with at least hgnc_id and ensembl_gene_id
#' @return A data frame with hgnc_id, ensemble_gene_id, and the gnomad constraint scores
#' @export
get_gnomad_constraint <- function(protein_coding_genes) {
  # latest gnomad release
  gnomad <- read_delim('https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv') %>%
    dplyr::rename('ensembl_gene_id' = 'gene_id') %>%
    dplyr::filter(canonical == 'TRUE')
  # joining with hgnc_ids
  gnomad_genes <- protein_coding_genes %>%
    left_join(gnomad, by = 'ensembl_gene_id') %>%
    dplyr::select(-`ensembl_gene_id`, -`mane_select`, -`gene`, -`transcript`,
                  -`canonical`, -`symbol`)

  return(gnomad_genes)
}
