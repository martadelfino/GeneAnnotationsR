


load('./data/mouse_gene_expression.RData')

load('./data/human_gene_expression.RData')


usethis::use_data(mouse, human, overwrite = TRUE)
