
#' Helper function - Downloads GO IDs
#'
#'
#' @param protein_coding_genes A dataframe with protein coding genes. Must have at least hgnc_id and entrez_id columns.
#' @return A dataframe with the GO ids of the protein coding genes.
#' @export
download_go_ids <- function(protein_coding_genes) {

  full_genes_entrez <- protein_coding_genes %>%
    dplyr::select(hgnc_id, entrez_id) %>%
    mutate(entrez_id = as.character(entrez_id))

  full_go_anot <- toTable(org.Hs.egGO) %>%
    filter(gene_id %in% full_genes_entrez$entrez_id) %>%
    dplyr::select(gene_id, go_id, Ontology) %>%
    distinct()

  return(full_go_anot)
}


#' Helper function - Downloads GO terms
#'
#'
#' @param protein_coding_genes A dataframe with protein coding genes. Must have at least hgnc_id and entrez_id columns.
#' @param full_go_anot A dataframe with GO IDs from the download_go_ids function.
#' @return A dataframe with the GO terms of the protein coding genes.
#' @export
download_go_terms <- function(protein_coding_genes, full_go_anot) {

  full_genes_entrez <- protein_coding_genes %>%
    dplyr::select(hgnc_id, entrez_id) %>%
    mutate(entrez_id = as.character(entrez_id))

  keys <- unique(full_go_anot$go_id)

  full_go_terms <- AnnotationDbi::select(GO.db, keys = keys, keytype="GOID", columns=c("TERM") ) %>%
    dplyr::rename(go_id = GOID,
                  go_term = TERM)

  full_go_anot_term <- full_go_anot %>%
    left_join(full_go_terms) %>%
    dplyr::group_by(gene_id, Ontology) %>%
    summarise(go_ids = paste0(go_id, collapse = "|"),
              go_terms = paste0(go_term, collapse = "|")) %>%
    pivot_wider(names_from = Ontology,
                values_from = c("go_ids","go_terms")) %>%
    dplyr::select(gene_id, go_ids_BP, go_terms_BP,
                  go_ids_MF, go_terms_MF,
                  go_ids_CC, go_terms_CC) #%>%
  #replace(is.na(.),"-")

  full_list_go_anot <- full_genes_entrez %>%
    left_join(full_go_anot_term, by = c("entrez_id" = "gene_id")) %>%
    dplyr::select(-entrez_id) #%>%
  #replace(is.na(.),"-")

  return(full_list_go_anot)
}


#' Main function - gets GO terms and IDs
#'
#'
#' @param protein_coding_genes A dataframe with protein coding genes. Must have at least hgnc_id and entrez_id columns.
#' @return A dataframe with the GO terms + ids of the protein coding genes.
#' @export
get_go <- function(protein_coding_genes) {

  full_go_anot <- download_go_ids(protein_coding_genes)
  full_list_go_anot <- download_go_terms(protein_coding_genes, full_go_anot)

  return(full_list_go_anot)
}

