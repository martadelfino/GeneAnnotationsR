# File: data-raw/generate_exomiser_pheno_score.R

load('./data/NDD_Exomiser_Pheno_Score.RData')

usethis::use_data(exomiser_pheno_score, overwrite = TRUE)
