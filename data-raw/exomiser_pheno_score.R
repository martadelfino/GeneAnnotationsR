
load('./data/NDD_HiPhive_Prioritiser_genes.RData')

exomiser_pheno_score <- df

usethis::use_data(exomiser_pheno_score, overwrite = TRUE)
