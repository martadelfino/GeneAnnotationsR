
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
#' @return A df of the GEL panels for NDD. Green confidence.
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
#' @return Df of the Australia panels NDD. Green confidence.
#' @export
get_australia_panelapp <- function(protein_coding_genes) {

  panel_urls <- c(
    "https://panelapp-aus.org/panels/20/download/01234/",
    "https://panelapp-aus.org/panels/112/download/01234/",
    "https://panelapp-aus.org/panels/250/download/01234/",
    "https://panelapp-aus.org/panels/3136/download/01234/"
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


#' Main Function 4 - amber conf gel panelapp
#'
#'
#'
#' @return A df of the GEL panels for NDD. Amber confidence.
#' @export
get_gel_panelapp_amber <- function(protein_coding_genes) {

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


#' Main Function 5 - red conf gel panelapp
#'
#'
#'
#' @return A df of the GEL panels for NDD. Red confidence.
#' @export
get_gel_panelapp_red <- function(protein_coding_genes) {

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


#' Main Function 6 - amber conf aus panelapp
#'
#'
#'
#' @return Df of the Australia panels NDD. Amber confidence.
#' @export
get_australia_panelapp_amber <- function(protein_coding_genes) {

  panel_urls <- c(
    "https://panelapp-aus.org/panels/20/download/01234/",
    "https://panelapp-aus.org/panels/112/download/01234/",
    "https://panelapp-aus.org/panels/250/download/01234/",
    "https://panelapp-aus.org/panels/3136/download/01234/"
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

#' Main Function 7 - red conf aus panelapp
#'
#'
#'
#' @return Df of the Australia panels NDD. Red confidence.
#' @export
get_australia_panelapp_red <- function(protein_coding_genes) {

  panel_urls <- c(
    "https://panelapp-aus.org/panels/20/download/01234/",
    "https://panelapp-aus.org/panels/112/download/01234/",
    "https://panelapp-aus.org/panels/250/download/01234/",
    "https://panelapp-aus.org/panels/3136/download/01234/"
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



#' Main Function 8 - df of SysNDD definitive
#'
#'
#'
#' @return df of SysNDD. Definitive confidence.
#' @export
get_sysndd <- function(protein_coding_genes) {

  df <- sysndd_genes %>%
    dplyr::filter(entities_category == 'Definitive') %>%
    distinct()

  df_agg <- df %>%
    dplyr::group_by(hgnc_id) %>%
    dplyr::summarise(entities_inheritance_filter = paste0(unique(entities_inheritance_filter), collapse = "|"))

  df3 <- df_agg %>%
    dplyr::mutate(
      select_gene_confidence_sysndd = '1',
      select_gene_moi_sysndd_ad = ifelse(grepl("Autosomal dominant", entities_inheritance_filter) & select_gene_confidence_sysndd == '1', '1', NA),
      select_gene_moi_sysndd_ar = ifelse(grepl("Autosomal recessive", entities_inheritance_filter) & select_gene_confidence_sysndd == '1', '1', NA),
      select_gene_moi_sysndd_adar = ifelse(grepl("Autosomal dominant|Autosomal recessive", entities_inheritance_filter) & select_gene_confidence_sysndd == '1', '1', NA),
      select_gene_moi_sysndd_all = ifelse(select_gene_confidence_sysndd == '1','1', NA))

  protein_coding_genes_sysndd <- protein_coding_genes %>%
    left_join(df3, by = 'hgnc_id')

  return(protein_coding_genes_sysndd)
  }


#' Main Function 9 - df of SysNDD moderate
#'
#'
#'
#' @return df of SysNDD. Moderate confidence.
#' @export
get_sysndd_moderate <- function(protein_coding_genes) {

  df <- sysndd_genes %>%
    dplyr::filter(entities_category == 'Moderate') %>%
    distinct()

  df_agg <- df %>%
    dplyr::group_by(hgnc_id) %>%
    dplyr::summarise(entities_inheritance_filter = paste0(unique(entities_inheritance_filter), collapse = "|"))

  df3 <- df_agg %>%
    dplyr::mutate(
      select_gene_confidence_sysndd = '1',
      select_gene_moi_sysndd_ad = ifelse(grepl("Autosomal dominant", entities_inheritance_filter) & select_gene_confidence_sysndd == '1', '1', NA),
      select_gene_moi_sysndd_ar = ifelse(grepl("Autosomal recessive", entities_inheritance_filter) & select_gene_confidence_sysndd == '1', '1', NA),
      select_gene_moi_sysndd_adar = ifelse(grepl("Autosomal dominant|Autosomal recessive", entities_inheritance_filter) & select_gene_confidence_sysndd == '1', '1', NA),
      select_gene_moi_sysndd_all = ifelse(select_gene_confidence_sysndd == '1','1', NA))

  protein_coding_genes_sysndd <- protein_coding_genes %>%
    left_join(df3, by = 'hgnc_id')

  return(protein_coding_genes_sysndd)
}


#' Main Function 10 - df of SysNDD limited
#'
#'
#'
#' @return df of SysNDD. Limited confidence.
#' @export
get_sysndd_limited <- function(protein_coding_genes) {

  df <- sysndd_genes %>%
    dplyr::filter(entities_category == 'Limited') %>%
    distinct()

  df_agg <- df %>%
    dplyr::group_by(hgnc_id) %>%
    dplyr::summarise(entities_inheritance_filter = paste0(unique(entities_inheritance_filter), collapse = "|"))

  df3 <- df_agg %>%
    dplyr::mutate(
      select_gene_confidence_sysndd = '1',
      select_gene_moi_sysndd_ad = ifelse(grepl("Autosomal dominant", entities_inheritance_filter) & select_gene_confidence_sysndd == '1', '1', NA),
      select_gene_moi_sysndd_ar = ifelse(grepl("Autosomal recessive", entities_inheritance_filter) & select_gene_confidence_sysndd == '1', '1', NA),
      select_gene_moi_sysndd_adar = ifelse(grepl("Autosomal dominant|Autosomal recessive", entities_inheritance_filter) & select_gene_confidence_sysndd == '1', '1', NA),
      select_gene_moi_sysndd_all = ifelse(select_gene_confidence_sysndd == '1','1', NA))

  protein_coding_genes_sysndd <- protein_coding_genes %>%
    left_join(df3, by = 'hgnc_id')

  return(protein_coding_genes_sysndd)
}


#### ddg2p genes now

#' Main Function 11 - df of ddg2p
#'
#'
#'
#' @return df of ddg2p. Definitive and strong confidence. Brain/Cognition only.
#' @export
get_ddg2p <- function(protein_coding_genes) {

  df <- ddg2p_genes %>%
    dplyr::select(`hgnc id`, `confidence category`, `allelic requirement`, `organ specificity list`) %>%
    dplyr::rename(hgnc_id = `hgnc id`,
                  confidence = `confidence category`,
                  allelic_requirement = `allelic requirement`,
                  organ_specificity = `organ specificity list`) %>%
    dplyr::mutate(hgnc_id = paste0("HGNC:", hgnc_id)) %>%
    dplyr::filter(grepl("Brain/Cognition", organ_specificity)) %>%
    dplyr::select(hgnc_id, confidence, allelic_requirement) %>%
    dplyr::filter(confidence == 'definitive' | confidence == 'strong' ) %>%
    distinct()

  df_agg <- df %>%
    dplyr::group_by(hgnc_id) %>%
    dplyr::summarise(allelic_requirement = paste0(unique(allelic_requirement), collapse = "|"))

  df2 <- df_agg %>%
    dplyr::mutate(
      select_gene_confidence_ddg2p = '1',
      select_gene_moi_ddg2p_ad = ifelse(grepl("monoallelic_autosomal", allelic_requirement) & select_gene_confidence_ddg2p == '1', '1', NA),
      select_gene_moi_ddg2p_ar = ifelse(grepl("biallelic_autosomal", allelic_requirement) & select_gene_confidence_ddg2p == '1', '1', NA),
      select_gene_moi_ddg2p_adar = ifelse(grepl("monoallelic_autosomal|biallelic_autosomal", allelic_requirement) & select_gene_confidence_ddg2p == '1','1', NA),
      select_gene_moi_ddg2p_all = ifelse(select_gene_confidence_ddg2p == '1','1', NA))

  protein_coding_genes_ddg2p <- protein_coding_genes %>%
    left_join(df2, by = 'hgnc_id')

  return(protein_coding_genes_ddg2p)
}


#' Main Function 12 - df of ddg2p moderate
#'
#'
#'
#' @return df of ddg2p. Moderate confidence. Brain/Cognition only.
#' @export
get_ddg2p_moderate <- function(protein_coding_genes) {

  df <- ddg2p_genes %>%
    dplyr::select(`hgnc id`, `confidence category`, `allelic requirement`, `organ specificity list`) %>%
    dplyr::rename(hgnc_id = `hgnc id`,
                  confidence = `confidence category`,
                  allelic_requirement = `allelic requirement`,
                  organ_specificity = `organ specificity list`) %>%
    dplyr::mutate(hgnc_id = paste0("HGNC:", hgnc_id)) %>%
    dplyr::filter(grepl("Brain/Cognition", organ_specificity)) %>%
    dplyr::select(hgnc_id, confidence, allelic_requirement) %>%
    dplyr::filter(confidence == 'moderate') %>%
    distinct()

  df_agg <- df %>%
    dplyr::group_by(hgnc_id) %>%
    dplyr::summarise(allelic_requirement = paste0(unique(allelic_requirement), collapse = "|"))

  df2 <- df_agg %>%
    dplyr::mutate(
      select_gene_confidence_ddg2p = '1',
      select_gene_moi_ddg2p_ad = ifelse(grepl("monoallelic_autosomal", allelic_requirement) & select_gene_confidence_ddg2p == '1', '1', NA),
      select_gene_moi_ddg2p_ar = ifelse(grepl("biallelic_autosomal", allelic_requirement) & select_gene_confidence_ddg2p == '1', '1', NA),
      select_gene_moi_ddg2p_adar = ifelse(grepl("monoallelic_autosomal|biallelic_autosomal", allelic_requirement) & select_gene_confidence_ddg2p == '1','1', NA),
      select_gene_moi_ddg2p_all = ifelse(select_gene_confidence_ddg2p == '1','1', NA))

  protein_coding_genes_ddg2p <- protein_coding_genes %>%
    left_join(df2, by = 'hgnc_id')

  return(protein_coding_genes_ddg2p)
}


#' Main Function 11 - df of ddg2p limited
#'
#'
#'
#' @return df of ddg2p. Limited confidence. Brain/Cognition only.
#' @export
get_ddg2p_limited <- function(protein_coding_genes) {

  df <- ddg2p_genes %>%
    dplyr::select(`hgnc id`, `confidence category`, `allelic requirement`, `organ specificity list`) %>%
    dplyr::rename(hgnc_id = `hgnc id`,
                  confidence = `confidence category`,
                  allelic_requirement = `allelic requirement`,
                  organ_specificity = `organ specificity list`) %>%
    dplyr::mutate(hgnc_id = paste0("HGNC:", hgnc_id)) %>%
    dplyr::filter(grepl("Brain/Cognition", organ_specificity)) %>%
    dplyr::select(hgnc_id, confidence, allelic_requirement) %>%
    dplyr::filter(confidence == 'limited') %>%
    distinct()

  df_agg <- df %>%
    dplyr::group_by(hgnc_id) %>%
    dplyr::summarise(allelic_requirement = paste0(unique(allelic_requirement), collapse = "|"))

  df2 <- df_agg %>%
    dplyr::mutate(
      select_gene_confidence_ddg2p = '1',
      select_gene_moi_ddg2p_ad = ifelse(grepl("monoallelic_autosomal", allelic_requirement) & select_gene_confidence_ddg2p == '1', '1', NA),
      select_gene_moi_ddg2p_ar = ifelse(grepl("biallelic_autosomal", allelic_requirement) & select_gene_confidence_ddg2p == '1', '1', NA),
      select_gene_moi_ddg2p_adar = ifelse(grepl("monoallelic_autosomal|biallelic_autosomal", allelic_requirement) & select_gene_confidence_ddg2p == '1','1', NA),
      select_gene_moi_ddg2p_all = ifelse(select_gene_confidence_ddg2p == '1','1', NA))

  protein_coding_genes_ddg2p <- protein_coding_genes %>%
    left_join(df2, by = 'hgnc_id')

  return(protein_coding_genes_ddg2p)
}

