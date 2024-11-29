## code to prepare `NDD_HiPhive_Prioritiser_genes` dataset goes here


# data-raw/NDD_HiPhive_Prioritiser_genes.R

load('./data/NDD_HiPhive_Prioritiser_genes.RData')

usethis::use_data(hiphive, overwrite = TRUE)
