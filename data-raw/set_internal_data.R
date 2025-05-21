###### Setting up internal data #####

# Exomiser phenotype data
load('./data/NDD_Exomiser_Pheno_Score.RData')

# Exomiser HiPhive data
load('./data/NDD_HiPhive_Prioritiser_genes.RData')

# Gene expression data
load('./data/mouse_gene_expression.RData')
load('./data/human_gene_expression.RData')

# Cell essentiality data
load('./data/cell_essentiality.rdata')

# Lethal terms
### this is a vector containing a set of mammalian phenotypes ids
### extracted from Dickinson et al.(PMID: 27626380): embryonic to
### preweaning lethality terms
load('./data/LethalTerms.RData')

# Other constraint metrics
load('./data/gene_constraint_metrics.rdata')

# SysNDD
load('./data/sysndd_genes.RData')

# DDG2P
load('./data/ddg2p_genes.RData')

# ProteomicsDB protein expression
load('./data/complete_protein_expression_data.RData')


# Creating internal data file
usethis::use_data(exomiser_pheno_score,
                  hiphive,
                  mouse, human,
                  cell_essentiality,
                  lethal_terms,
                  gene_constraint_metrics,
                  sysndd_genes,
                  ddg2p_genes,
                  all_protein_data,
                  internal = TRUE, overwrite = TRUE)


