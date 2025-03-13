# Install and load the readxl package if not already installed
if (!require("readxl")) {
  install.packages("readxl")
  library(readxl)
}

if (!require("tidyverse")) {
  install.packages("tidyverse")
  library(tidyverse)
}

# Define the path to your Excel file
#file_path <- "sysndd_panels.xlsx"
file_path <- "data-raw/sysndd_gene_table.xlsx"

# Read the Excel file into a data frame
df <- read_excel(file_path)

df2 <- df %>%
  dplyr::select(hgnc_id, entities_inheritance_filter, entities_category) %>%
  dplyr::filter(entities_category != 'not applicable')


# Saving file as RData
sysndd_genes <- df2
save(sysndd_genes, file = "data/sysndd_genes.RData")





#df_definitive <- df2 %>%
  #  dplyr::filter(entities_category == 'Definitive') %>%
  #  distinct()

#df3 <- df_definitive %>%
  #dplyr::mutate(
    #select_gene_confidence_sysndd = ifelse(grepl("Definitive", entities_category), '1', NA_real_),
  #select_gene_moi_sysndd_ad = ifelse(grepl("Autosomal dominant", entities_inheritance_filter), '1', NA_real_),
  #select_gene_moi_sysndd_ar = ifelse(grepl("Autosomal recessive", entities_inheritance_filter), '1', NA_real_),
  #  select_gene_moi_sysndd_adar = ifelse(grepl("Autosomal dominant|Autosomal recessive", entities_inheritance_filter), '1', NA_real_),
  #  select_gene_moi_sysndd_x = ifelse(grepl("X-linked", entities_inheritance_filter), '1', NA_real_),
#  select_gene_moi_sysndd_all = ifelse(grepl("Autosomal dominant|Autosomal recessive|X-linked", entities_inheritance_filter), '1', NA_real_),
#)

#df_moderate <- df2 %>%
 # dplyr::filter(entities_category == 'Moderate')

#df_limited <- df2 %>%
 # dplyr::filter(entities_category == 'Limited')





