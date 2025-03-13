# https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads_archive/2025_01_28/DDG2P_2025-01-28.csv.gz
# accessed on Thursday 13MAR25
# accessing the previous one (not the latest, 2025-02-28) because they haven't assigned organ to the latest
# merging together definitive and strong
# then keeping as it is moderate and limited


if (!require("tidyverse")) {
  install.packages("tidyverse")
  library(tidyverse)
}

data <- read_csv("data-raw/DDG2P_2025-01-28.csv.gz")

ddg2p_genes <- data
save(ddg2p_genes, file = "data/ddg2p_genes.RData")




#df <- data %>%
  #dplyr::select(`hgnc id`, `confidence category`, `allelic requirement`, `organ specificity list`) %>%
   # dplyr::rename(hgnc_id = `hgnc id`,
    #              confidence = `confidence category`,
     #             allelic_requirement = `allelic requirement`,
      #            organ_specificity = `organ specificity list`) %>%
  #  dplyr::mutate(hgnc_id = paste0("HGNC:", hgnc_id)) %>%
   # dplyr::filter(grepl("Brain/Cognition", organ_specificity)) %>%
    #dplyr::select(hgnc_id, confidence, allelic_requirement)



#df2 <- df %>%
 # dplyr::mutate(
  #  select_gene_confidence_ddg2p = ifelse(grepl("definitive|strong", confidence), '1', NA_real_),
   # select_gene_moi_ddg2p_ad = ifelse(grepl("monoallelic_autosomal", allelic_requirement), '1', NA_real_),
    #select_gene_moi_ddg2p_ar = ifelse(grepl("biallelic_autosomal", allelic_requirement), '1', NA_real_),
#    select_gene_moi_ddg2p_adar = ifelse(grepl("monoallelic_autosoma|biallelic_autosomal", allelic_requirement), '1', NA_real_),
 #   select_gene_moi_ddg2p_x = ifelse(grepl("monoallelic_X_het|monoallelic_X_hem", allelic_requirement), '1', NA_real_),
  #  select_gene_moi_ddg2p_all = ifelse(grepl("monoallelic_autosoma|biallelic_autosomal|monoallelic_X_het|monoallelic_X_hem", allelic_requirement), '1', NA_real_),
#  )
