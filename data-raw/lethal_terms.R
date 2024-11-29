
### this is a vector containing a set of mammalian phenotypes ids
### extracted from Dickinson et al.(PMID: 27626380): embryonic to
### preweaning lethality terms

load('./data/LethalTerms.RData')

usethis::use_data(lethal_terms, overwrite = TRUE)
