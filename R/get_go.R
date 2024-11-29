

if (!require("tidyverse")) install.packages("tidyverse")
library("tidyverse")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")

if (!require("AnnotationDbi", quietly = TRUE))
  BiocManager::install("AnnotationDbi")
library("AnnotationDbi")

if (!require("GO.db", quietly = TRUE))
  BiocManager::install("GO.db")
library(GO.db)

if (!require("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
library("clusterProfiler")




get_go_for_specified_genes <- function(input_genes) {
  source('./R/get_protein_coding_genes.R')

  protein_coding_genes <- get_protein_coding_genes() %>%
    dplyr::select(hgnc_id, entrez_id) %>%
    dplyr::mutate(entrez_id = as.character(entrez_id))

  input_genes <- input_genes %>%
    left_join(protein_coding_genes, by = 'hgnc_id')

  input_genes_go_annotation <- toTable(org.Hs.egGO) %>%
    dplyr::filter(gene_id %in% input_genes$entrez_id) %>%
    dplyr::select(gene_id, go_id, Ontology) %>%
    distinct()

  input_genes_keys <- unique(input_genes_go_annotation$go_id)

  input_genes_go_terms <- AnnotationDbi::select(GO.db,
                                               keys = input_genes_keys,
                                               keytype="GOID",
                                               columns=c("TERM")) %>%
    dplyr::rename(go_id = GOID, go_term = TERM)

  input_genes_go_anot_terms <- input_genes_go_annotation %>%
    left_join(input_genes_go_terms) %>%
    group_by(gene_id, Ontology) %>%
    summarise(go_ids = paste0(go_id, collapse = "|"),
              go_terms = paste0(go_term, collapse = "|")) %>%
    pivot_wider(names_from = Ontology,
                values_from = c("go_ids","go_terms")) %>%
    dplyr::select(gene_id, go_ids_BP, go_terms_BP,
                  go_ids_MF, go_terms_MF,
                  go_ids_CC, go_terms_CC)%>%
    replace(is.na(.),"-")

  input_genes_go <- input_genes %>%
    left_join(input_genes_go_anot_terms, by = c("entrez_id" = "gene_id")) %>%
    dplyr::select(-entrez_id) %>%
    replace(is.na(.),"-")

  return(input_genes_go)
}



#source('./R/get_protein_coding_genes.R')

#genes <- get_protein_coding_genes()
#genes500 <- genes %>%
 # dplyr::slice(1:500) %>%
  #dplyr::select(hgnc_id)

#test <- get_go(genes500)





get_go_for_all_genes <- function() {


  full_genes_entrez <- get_protein_coding_genes() %>%
    dplyr::select(hgnc_id, entrez_id) %>%
    mutate(entrez_id = as.character(entrez_id))

  full_go_anot <- toTable(org.Hs.egGO) %>%
    filter(gene_id %in% full_genes_entrez$entrez_id) %>%
    dplyr::select(gene_id, go_id, Ontology) %>%
    distinct()

  keys <- unique(full_go_anot$go_id)

  full_go_terms <- AnnotationDbi::select(GO.db, keys = keys, keytype="GOID", columns=c("TERM") ) %>%
    dplyr::rename(go_id = GOID,
                  go_term = TERM)

  full_go_anot_term <- full_go_anot %>%
    left_join(full_go_terms) %>%
    group_by(gene_id, Ontology) %>%
    summarise(go_ids = paste0(go_id, collapse = "|"),
              go_terms = paste0(go_term, collapse = "|")) %>%
    pivot_wider(names_from = Ontology,
                values_from = c("go_ids","go_terms")) %>%
    dplyr::select(gene_id, go_ids_BP, go_terms_BP,
                  go_ids_MF, go_terms_MF,
                  go_ids_CC, go_terms_CC)%>%
    replace(is.na(.),"-")

  full_list_go_anot <- full_genes_entrez %>%
    left_join(full_go_anot_term, by = c("entrez_id" = "gene_id")) %>%
    dplyr::select(-entrez_id) %>%
    replace(is.na(.),"-")

  return(full_go_anot)
}

source('./R/get_protein_coding_genes.R')
genes <- get_protein_coding_genes()
genes500 <- genes %>%
  dplyr::slice(1:500) %>%
  dplyr::select(hgnc_id)

test1 <- get_go_for_all_genes()






get_enriched_go <- function(input_genes) {
  source('./R/get_protein_coding_genes.R')
  genes_universe <- get_protein_coding_genes() %>%
    dplyr::select(hgnc_id, symbol, entrez_id)

  genes <- input_genes %>%
    inner_join(genes_universe, by = 'hgnc_id') %>%
    pull(entrez_id) %>%
    as.character(.)

  universe <- as.character(genes_universe$entrez_id)

  # enrichment
  genes_go <- enrichGO(gene          = genes,
                       universe      = universe,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)

  genes_go_df <- genes_go@result %>%
    dplyr::filter(Count > 5) %>%
    dplyr::filter(p.adjust < 0.01)

  return(genes_go_df)
}



source('./R/get_protein_coding_genes.R')
genes <- get_protein_coding_genes()
genes500 <- genes %>%
 dplyr::slice(1:500) %>%
 dplyr::select(hgnc_id)

test2 <- get_enriched_go(genes500)









full_gene_go_exploded <- get_go_for_all_genes()

source('./R/get_protein_coding_genes.R')
genes <- get_protein_coding_genes()
genes500 <- genes %>%
 dplyr::slice(1:50) %>%
 dplyr::select(hgnc_id)

ad_go <- get_enriched_go(genes500)

#source('./R/get_protein_coding_genes.R')
genes_universe <- get_protein_coding_genes() %>%
  dplyr::select(hgnc_id, symbol, entrez_id) %>%
  mutate(entrez_id = as.character(entrez_id))


# GO SCORE

# now remove all pathways that are not AD NDD

# Create a new DataFrame containing only the pathway IDs from NDDs
only_universe_with_ad_ndd_go = full_gene_go_exploded %>%
  dplyr::filter(go_id %in% ad_go$ID)

universe_with_ad_ndd_go <- left_join(genes_universe, only_universe_with_ad_ndd_go,
                                     by = join_by('entrez_id' == 'gene_id'))

universe_with_ad_ndd_go <- left_join(universe_with_ad_ndd_go, ad_go,
                                     by = join_by('go_id' == 'ID'))

length(unique(universe_with_ad_ndd_go$hgnc_id))
# 19261

length(unique(universe_with_ad_ndd_go$go_id))
# 131




###### SUM OF THE COUNTS FOR EACH NDD PATHWAY FOR EACH GENE

# Creating a new dataframe with sums per gene
sum_universe_with_ad_ndd_go <- universe_with_ad_ndd_go %>%
  group_by(hgnc_id) %>%
  summarize(go_ndd_score = sum(Count))

length(unique(sum_universe_with_ad_ndd_go$hgnc_id))
# 19247

sum(is.na(sum_universe_with_ad_ndd_go$go_ndd_score))
# 13291

# switching NAs for 0s
sum_universe_with_ad_ndd_go <- sum_universe_with_ad_ndd_go %>%
  replace(is.na(.), 0)



