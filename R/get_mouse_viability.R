
#' Helper Function 1 - get mouse protein coding genes
#'
#' This function performs a specific task. Getting mouse protein coding genes.
#'
#' @return Large character. Of MGI ids.
#' @export
get_mouse_protein_coding_genes <- function() {
  # Get whole set of mouse protein coding genes
  mouse_genes <- read_delim("https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt") %>%
    dplyr::select("MGI Accession ID","Feature Type") %>%
    dplyr::filter(`Feature Type` == "protein coding gene")%>%
    distinct()
  mouse_proteincoding_genes = unique(mouse_genes$`MGI Accession ID`)

  return(mouse_proteincoding_genes)
}


#' Helper Function 2 - get IMPC viability.
#'
#' This function performs a specific task. Getting IMPC viability.
#'
#' @return Dataframe of IMPC viability.
#' \describe{
#'   \item{mgi_id}
#'   \item{viability_impc}
#' }
#' @export
get_impc_viability <- function() {
  via_impc <- read_delim("http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/viability.csv.gz", delim = ",") %>%
    dplyr::filter(is.na(Comment)) %>%
    rename(gene_symbol_mm  = "Gene Symbol",
           mgi_id = "Gene Accession Id",
           viability_impc = "Viability Phenotype HOMs/HEMIs") %>%
    dplyr::select(mgi_id,viability_impc) %>%
    distinct()

  return(via_impc)
}


#' Helper Function 3 - get MGI viability.
#'
#' This function performs a specific task. Getting MGI viability.
#'
#' @param mouse_proteincoding_genes Output of get_mouse_protein_coding_genes()
#' @param lethal_terms lethal terms.
#' @return Dataframe of MGI viability.
#' \describe{
#'   \item{mgi_id}
#'   \item{viability_mgi}
#' }
#' @export
get_mgi_viability <- function(mouse_proteincoding_genes, lethal_terms) {
  via_mgi <- readr::read_delim("https://www.informatics.jax.org/downloads/reports/MGI_GenePheno.rpt",  col_names = FALSE) %>%
    dplyr::select(7,5) %>%
    distinct() %>%
    rename(mgi_id = X7, mp_term = X5)  %>%
    dplyr::filter(mgi_id %in% mouse_proteincoding_genes) %>%
    dplyr::mutate(mp_term_lethal = ifelse(mp_term %in% lethal_terms, "y","n")) %>% # lethal_terms is a global variable. will need to be changed to a function output
    dplyr::select(mgi_id, mp_term_lethal) %>%
    distinct() %>%
    arrange(mgi_id, mp_term_lethal) %>%
    group_by(mgi_id) %>%
    summarise(mgi_lethal_term = paste0(unique(mp_term_lethal), collapse = "|")) %>%
    mutate(viability_mgi = ifelse(mgi_lethal_term == "n","viable","lethal")) %>%
    dplyr::select(mgi_id, viability_mgi) %>%
    distinct()

  return(via_mgi)
}


#' Main Function
#'
#' This is the main function which uses helper functions to obtain the combined mouse viability.
#'
#' @param protein_coding_genes Output of get_protein_coding_genes()
#' @return A data frame with the combined mouse viability.
#' @export
get_mouse_viability <- function(protein_coding_genes) {
  # Select only gene symbol and mgi id columns
  hgnc <- protein_coding_genes %>%
    dplyr::rename(gene_symbol = symbol, mgi_id = mgd_id) %>%
    dplyr::select(hgnc_id, mgi_id)%>%
    separate_rows(mgi_id, sep = "\\|") %>%
    dplyr::filter(!is.na(mgi_id))

  # Select unique hgnc ids and mgi ids
  hgnc_dups <- unique(hgnc$hgnc_id[duplicated(hgnc$hgnc_id)])
  mgi_dups <- unique(hgnc$mgi_id[duplicated(hgnc$mgi_id)])

  # Orthologs. Keep only one2one relationships (strict criteria)
  one2one <- hgnc %>%
    dplyr::filter(!mgi_id %in% mgi_dups) %>%
    dplyr::filter(!hgnc_id %in% hgnc_dups)

  # Get whole set of mouse protein coding genes and keep unique ones
  mouse_proteincoding_genes <- get_mouse_protein_coding_genes()

  ### viability data impc
  ### lethal genes are easy to retrieve, there is a viability report
  ### three viability outcomes: lethal, subviable and viable
  via_impc <- get_impc_viability()

  ### get lethal phenotypes from the MGI based on the set of lethal
  ### terms imported through the RData object ("lethal_terms")
  via_mgi <- get_mgi_viability(mouse_proteincoding_genes, lethal_terms)

  ## this file contains single gene knockout viability outcomes from
  ### the impc and mgi resources
  viability_impc_mgi <- one2one %>%
    left_join(via_impc) %>%
    left_join(via_mgi) %>%
    replace(is.na(.),"-") %>%
    dplyr::filter(viability_impc != "-" | viability_mgi != "-")

  return(viability_impc_mgi)
}
