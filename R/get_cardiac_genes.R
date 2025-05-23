
#' Main Function 1 - df GEL panels cardiac
#'
#'
#'
#' @return Df of GEL cardiac panels.
#' @export
get_gel_panelapp_cardiac<- function(protein_coding_genes) {

  panel_urls <- c(
    "https://panelapp.genomicsengland.co.uk/panels/212/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/134/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/652/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/47/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/49/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/238/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/13/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/842/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/843/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/214/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/277/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/76/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/45/download/01234/"
  )

  # Load and clean each panel, then bind them together
  panelapp <- map_dfr(panel_urls, load_and_clean_panel) %>%
    dplyr::group_by(hgnc_id) %>%
    dplyr::summarise(
      status_gel = paste0(unique(status_gel), collapse = "|"),
      moi_gel = paste0(unique(moi_gel), collapse = "|"),
      phenotype_gel = paste0(unique(phenotype_gel), collapse = "|"),
      level4_gel = paste0(unique(level4_gel), collapse = "|"),
      level3_gel = paste0(unique(level3_gel), collapse = "|")
    ) %>%
    dplyr::mutate(
      select_gene_confidence_gel = ifelse(grepl(3, status_gel), '1', NA_real_),
      select_gene_moi_gel = ifelse(grepl("MONOALLELIC, autosomal or pseudoautosomal|BOTH monoallelic and biallelic", moi_gel), '1', NA_real_),
      select_gene_moi_gel_ad = ifelse(grepl("MONOALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_ar = ifelse(grepl("BIALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_adar = ifelse(grepl("MONOALLELIC, autosomal|BIALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_all = ifelse(select_gene_confidence_gel == '1', '1', NA_real_)
    ) #%>%
  #replace(is.na(.), "-")

  protein_coding_genes_panelapp <- protein_coding_genes %>%
    left_join(panelapp, by = 'hgnc_id')

  return(protein_coding_genes_panelapp)
}


#' Main Function 2 - df GEL panels cardiac amber
#'
#'
#'
#' @return Df of GEL cardiac panels amber.
#' @export
get_gel_panelapp_amber_cardiac <- function(protein_coding_genes) {

  panel_urls <- c(
    "https://panelapp.genomicsengland.co.uk/panels/212/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/134/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/652/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/47/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/49/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/238/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/13/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/842/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/843/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/214/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/277/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/76/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/45/download/01234/"
  )

  # Load and clean each panel, then bind them together
  panelapp <- map_dfr(panel_urls, load_and_clean_panel) %>%
    dplyr::group_by(hgnc_id) %>%
    dplyr::summarise(
      status_gel = paste0(unique(status_gel), collapse = "|"),
      moi_gel = paste0(unique(moi_gel), collapse = "|"),
      phenotype_gel = paste0(unique(phenotype_gel), collapse = "|"),
      level4_gel = paste0(unique(level4_gel), collapse = "|"),
      level3_gel = paste0(unique(level3_gel), collapse = "|")
    ) %>%
    dplyr::mutate(
      select_gene_confidence_gel = ifelse(grepl(2, status_gel), '1', NA_real_),
      select_gene_moi_gel = ifelse(grepl("MONOALLELIC, autosomal or pseudoautosomal|BOTH monoallelic and biallelic", moi_gel), '1', NA_real_),
      select_gene_moi_gel_ad = ifelse(grepl("MONOALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_ar = ifelse(grepl("BIALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_adar = ifelse(grepl("MONOALLELIC, autosomal|BIALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_all = ifelse(select_gene_confidence_gel == '1', '1', NA_real_)
    ) #%>%
  #replace(is.na(.), "-")

  protein_coding_genes_panelapp <- protein_coding_genes %>%
    left_join(panelapp, by = 'hgnc_id')

  return(protein_coding_genes_panelapp)
}


#' Main Function 3 - df GEL panels cardiac red
#'
#'
#'
#' @return Df of GEL cardiac panels red.
#' @export
get_gel_panelapp_red_cardiac <- function(protein_coding_genes) {

  panel_urls <- c(
    "https://panelapp.genomicsengland.co.uk/panels/212/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/134/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/652/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/47/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/49/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/238/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/13/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/842/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/843/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/214/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/277/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/76/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/45/download/01234/"
  )

  # Load and clean each panel, then bind them together
  panelapp <- map_dfr(panel_urls, load_and_clean_panel) %>%
    dplyr::group_by(hgnc_id) %>%
    dplyr::summarise(
      status_gel = paste0(unique(status_gel), collapse = "|"),
      moi_gel = paste0(unique(moi_gel), collapse = "|"),
      phenotype_gel = paste0(unique(phenotype_gel), collapse = "|"),
      level4_gel = paste0(unique(level4_gel), collapse = "|"),
      level3_gel = paste0(unique(level3_gel), collapse = "|")
    ) %>%
    dplyr::mutate(
      select_gene_confidence_gel = ifelse(grepl(1, status_gel), '1', NA_real_),
      select_gene_moi_gel = ifelse(grepl("MONOALLELIC, autosomal or pseudoautosomal|BOTH monoallelic and biallelic", moi_gel), '1', NA_real_),
      select_gene_moi_gel_ad = ifelse(grepl("MONOALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_ar = ifelse(grepl("BIALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_adar = ifelse(grepl("MONOALLELIC, autosomal|BIALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_all = ifelse(select_gene_confidence_gel == '1', '1', NA_real_)
    ) #%>%
  #replace(is.na(.), "-")

  protein_coding_genes_panelapp <- protein_coding_genes %>%
    left_join(panelapp, by = 'hgnc_id')

  return(protein_coding_genes_panelapp)
}


#' Main Function 4 - df of Australia cardiac panels
#'
#'
#'
#' @return Df of austalia cardiac panels.
#' @export
get_australia_panelapp_cardiac <- function(protein_coding_genes) {

  panel_urls <- c(
    "https://panelapp-aus.org/panels/76/download/01234/",
    "https://panelapp-aus.org/panels/111/download/01234/",
    "https://panelapp-aus.org/panels/95/download/01234/",
    "https://panelapp-aus.org/panels/3270/download/01234/",
    "https://panelapp-aus.org/panels/48/download/01234/",
    "https://panelapp-aus.org/panels/60/download/01234/",
    "https://panelapp-aus.org/panels/92/download/01234/",
    "https://panelapp-aus.org/panels/183/download/01234/",
    "https://panelapp-aus.org/panels/131/download/01234/"
  )

  # Load and clean each panel, then bind them together
  panelapp <- map_dfr(panel_urls, load_and_clean_panel) %>%
    dplyr::group_by(hgnc_id) %>%
    dplyr::summarise(
      status_gel = paste0(unique(status_gel), collapse = "|"),
      moi_gel = paste0(unique(moi_gel), collapse = "|"),
      phenotype_gel = paste0(unique(phenotype_gel), collapse = "|"),
      level4_gel = paste0(unique(level4_gel), collapse = "|"),
      level3_gel = paste0(unique(level3_gel), collapse = "|")
    ) %>%
    dplyr::mutate(
      select_gene_confidence_gel = ifelse(grepl(3, status_gel), '1', NA_real_),
      select_gene_moi_gel = ifelse(grepl("MONOALLELIC, autosomal or pseudoautosomal|BOTH monoallelic and biallelic", moi_gel), '1', NA_real_),
      select_gene_moi_gel_ad = ifelse(grepl("MONOALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_ar = ifelse(grepl("BIALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_adar = ifelse(grepl("MONOALLELIC, autosomal|BIALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_all = ifelse(select_gene_confidence_gel == '1', '1', NA_real_)
    ) #%>%
  #replace(is.na(.), "-")

  protein_coding_genes_panelapp <- protein_coding_genes %>%
    left_join(panelapp, by = 'hgnc_id')

  return(protein_coding_genes_panelapp)
}


#' Main Function 5 - df of Australia cardiac panels amber
#'
#'
#'
#' @return Df of austalia cardiac panels amber.
#' @export
get_australia_panelapp_amber_cardiac <- function(protein_coding_genes) {

  panel_urls <- c(
    "https://panelapp-aus.org/panels/76/download/01234/",
    "https://panelapp-aus.org/panels/111/download/01234/",
    "https://panelapp-aus.org/panels/95/download/01234/",
    "https://panelapp-aus.org/panels/3270/download/01234/",
    "https://panelapp-aus.org/panels/48/download/01234/",
    "https://panelapp-aus.org/panels/60/download/01234/",
    "https://panelapp-aus.org/panels/92/download/01234/",
    "https://panelapp-aus.org/panels/183/download/01234/",
    "https://panelapp-aus.org/panels/131/download/01234/"
  )

  # Load and clean each panel, then bind them together
  panelapp <- map_dfr(panel_urls, load_and_clean_panel) %>%
    dplyr::group_by(hgnc_id) %>%
    dplyr::summarise(
      status_gel = paste0(unique(status_gel), collapse = "|"),
      moi_gel = paste0(unique(moi_gel), collapse = "|"),
      phenotype_gel = paste0(unique(phenotype_gel), collapse = "|"),
      level4_gel = paste0(unique(level4_gel), collapse = "|"),
      level3_gel = paste0(unique(level3_gel), collapse = "|")
    ) %>%
    dplyr::mutate(
      select_gene_confidence_gel = ifelse(grepl(2, status_gel), '1', NA_real_),
      select_gene_moi_gel = ifelse(grepl("MONOALLELIC, autosomal or pseudoautosomal|BOTH monoallelic and biallelic", moi_gel), '1', NA_real_),
      select_gene_moi_gel_ad = ifelse(grepl("MONOALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_ar = ifelse(grepl("BIALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_adar = ifelse(grepl("MONOALLELIC, autosomal|BIALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_all = ifelse(select_gene_confidence_gel == '1', '1', NA_real_)
    ) #%>%
  #replace(is.na(.), "-")

  protein_coding_genes_panelapp <- protein_coding_genes %>%
    left_join(panelapp, by = 'hgnc_id')

  return(protein_coding_genes_panelapp)
}


#' Main Function 6 - df of Australia cardiac panels red
#'
#'
#'
#' @return Df of austalia cardiac panels red.
#' @export
get_australia_panelapp_red_cardiac <- function(protein_coding_genes) {

  panel_urls <- c(
    "https://panelapp-aus.org/panels/76/download/01234/",
    "https://panelapp-aus.org/panels/111/download/01234/",
    "https://panelapp-aus.org/panels/95/download/01234/",
    "https://panelapp-aus.org/panels/3270/download/01234/",
    "https://panelapp-aus.org/panels/48/download/01234/",
    "https://panelapp-aus.org/panels/60/download/01234/",
    "https://panelapp-aus.org/panels/92/download/01234/",
    "https://panelapp-aus.org/panels/183/download/01234/",
    "https://panelapp-aus.org/panels/131/download/01234/"
  )

  # Load and clean each panel, then bind them together
  panelapp <- map_dfr(panel_urls, load_and_clean_panel) %>%
    dplyr::group_by(hgnc_id) %>%
    dplyr::summarise(
      status_gel = paste0(unique(status_gel), collapse = "|"),
      moi_gel = paste0(unique(moi_gel), collapse = "|"),
      phenotype_gel = paste0(unique(phenotype_gel), collapse = "|"),
      level4_gel = paste0(unique(level4_gel), collapse = "|"),
      level3_gel = paste0(unique(level3_gel), collapse = "|")
    ) %>%
    dplyr::mutate(
      select_gene_confidence_gel = ifelse(grepl(1, status_gel), '1', NA_real_),
      select_gene_moi_gel = ifelse(grepl("MONOALLELIC, autosomal or pseudoautosomal|BOTH monoallelic and biallelic", moi_gel), '1', NA_real_),
      select_gene_moi_gel_ad = ifelse(grepl("MONOALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_ar = ifelse(grepl("BIALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_adar = ifelse(grepl("MONOALLELIC, autosomal|BIALLELIC, autosomal|BOTH monoallelic and biallelic", moi_gel) & select_gene_confidence_gel == '1', '1', NA_real_),
      select_gene_moi_gel_all = ifelse(select_gene_confidence_gel == '1', '1', NA_real_)
    ) #%>%
  #replace(is.na(.), "-")

  protein_coding_genes_panelapp <- protein_coding_genes %>%
    left_join(panelapp, by = 'hgnc_id')

  return(protein_coding_genes_panelapp)
}
