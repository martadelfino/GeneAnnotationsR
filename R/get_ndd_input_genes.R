
#### GENETREK GENE LIST COLUMNS:
#	ACMG SF v3.2 updated 2023/06/22
# Candidate NDD genes v1
# Candidate NDD genes v2
# Candidate NDD genes v3
# Candidate NDD genes v4
# Candidate NDD genes v5
# Chromatin Remodling Genes
# DBD Tier1 & AR genes (updated 2022-09-01)
# DBD Tier234 genes (updated 2022-09-01)
# DDg2p candidate|Cognition|brain genes (updated 2024-01-10)
# DDg2p definitive|Cognition|brain genes (updated 2024-01-10)
# Decreased Head Circumference HP_0040195 (HPO release 2024-01-16)
# Early onset or syndromic epilepsy v4 - Genomics England - 2023-03-22
# Epifactors v2.0 (updated 2022-07-07)
# Epilepsy genes - EPPAN (Mayo Clinic laboratories) - 2023-05-30
# Epilepsy genes tier 1 doi.org/10.1016/j.ejpn.2022.12.005
# Epilepsy genes tier 2 doi.org/10.1016/j.ejpn.2022.12.005
# Epilepsy HPO_1250 (HPO release 2024-01-16)
# gnomad.v4.0.constraint_metrics - lof.oe_ci.upper < 0.6
# GTEX (2017-06-05_v8 release)
# Haploinsufficient genes n= 2,987 (pHaplo >= 0.86)
# High Confidence Epilepsy Genes
# High Confidence NDD genes v1
# High Confidence NDD genes v2
# High Confidence NDD genes v3
# Increased Head Circumference HP_0040194 (HPO release 2024-01-16)
# pLI superior to 0.9 (gnomad.v2.1.1)
# postbirth BrainSpan
# prebirth BrainSpan
# SFARI_1 genes (updated 2024-01-16)
# SFARI_23S genes (updated 2024-01-16)
# SPARK genes (updated september 2022
# SynGO v1.2 (release 2023-12-01
# SysNDD candidate genes (updated 2023-11-09)
# SysNDD definitive genes (updated 2023-11-09)
# Triplosensitive genes n=1,559 (pTriplo >= 0.94)
# AUTISM|DD categories from Fu et al. Nat Genet 2022 (PMID: 35982160)
# AUTISM|DD categories from Wang et al. PNAS 2022 (PMID: 36350923)
# High-confidence autism genes from Zhou et al. Nat Genet 2022 (PMID: 35982159)
# DDg2p Inheritance Cognition|brain genes (updated 2024-01-10)
# SysNDD Inheritance (updated 2023-11-09)
# Inheritance type combined v3
# EAGLE Score from SFARI website
# pLI
# LOEUF
# MOEUF
# NCI|CI categories from SPARK-iWES-v1
# Autism-associated genes or constrained genes
# Prevalence autism
# Prevalence controls
# Number of carriers among autistic individuals
# Autism odds ratio
# p-value
# Confidence interval lower bound
# Confidence interval upper bound




#' Helper Function 1 - download genetrek file
#'
#' This function performs a specific task. Getting the latest genetrek file.
#'
#' @return Df of genetrek table.
#' @export
download_genetrek <- function() {

  # URL of the TSV file
  url <- "https://genetrek.pasteur.fr/downloadAllData?filetype=tsv"

  df <- read_tsv(url)

  return(df)
}


#' Main Function 2 - cleans genetrek file and gives you results for the genes specified
#'
#'
#'
#' @return Df of genetrek for specified genes.
#' @export
get_genetrek <- function(protein_coding_genes, columns = NULL) {

  genetrek <- download_genetrek()

  # Conditionally select columns if specified
  if (!is.null(columns)) {
    genetrek <- genetrek %>%
      dplyr::select(!!!rlang::syms(columns))
  }
  #print(genetrek)
  genetrek <- genetrek %>%
    dplyr::select(-`NCBI ID`, -`Ensembl ID`, `Gene`) %>%
    dplyr::rename(hgnc_id = `HGNC ID`) %>%
    dplyr::mutate(hgnc_id = ifelse(is.na(hgnc_id), NA,
                                   paste0("HGNC:", hgnc_id)))

  hgnc_genetrek <- protein_coding_genes %>%
    left_join(genetrek, by = 'hgnc_id') #%>%
    #dplyr::filter(`Gene type` == 'protein-coding')

  return(hgnc_genetrek)
}


#' Helper Function 1 - load and clean panels
#'
#' This function performs a specific task. Getting the latest genetrek file.
#'
#' @return cleaning of the panel.
#' @export
load_and_clean_panel <- function(url) {
  read_delim(url, delim = "\t", col_names = TRUE,
             col_types = cols(.default = "c")) %>%
    dplyr::filter(`Entity type` == "gene") %>%
    dplyr::rename(
      hgnc_id = HGNC,
      status_gel = GEL_Status,
      moi_gel = Model_Of_Inheritance,
      phenotype_gel = Phenotypes,
      level4_gel = Level4,
      level3_gel = Level3
    ) %>%
    dplyr::mutate(
      phenotype_gel = gsub("\\|", ",", phenotype_gel),
      phenotype_gel = gsub("\\t", ",", phenotype_gel)
    ) %>% #replace(is.na(.), "-") %>%
    dplyr::select(hgnc_id, status_gel, moi_gel, phenotype_gel, level4_gel,
                  level3_gel)
}


#' Main Function 2 - df of GEL panels NDD
#'
#'
#'
#' @return A df of the GEL panels for NDD.
#' @export
get_gel_panelapp <- function(protein_coding_genes) {

  panel_urls <- c(
    "https://panelapp.genomicsengland.co.uk/panels/285/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/197/download/01234",
    "https://panelapp.genomicsengland.co.uk/panels/78/download/01234/",
    "https://panelapp.genomicsengland.co.uk/panels/96/download/01234/"
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





#' Main Function 3 - df of Australia Panels NDD
#'
#'
#'
#' @return Df of the Australia panels NDD
#' @export
get_australia_panelapp <- function(protein_coding_genes) {

  panel_urls <- c(
    "https://panelapp.agha.umccr.org/panels/20/download/01234/",
    "https://panelapp.agha.umccr.org/panels/122/download/01234/",
    "https://panelapp.agha.umccr.org/panels/250/download/01234/",
    "https://panelapp.agha.umccr.org/panels/3136/download/01234/"
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

