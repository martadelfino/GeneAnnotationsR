
#' Helper Function 1 - get latest DepMap release links
#'
#' This function performs a specific task. Getting the latest DepMap release links.
#'
#' @return List of links to the latest DepMap release of Model.csv and CRISPRGeneEffect.csv
#' @export
get_latest_depmap_release <- function() {
  url <- "https://depmap.org/portal/api/download/files"
  files <- read_csv(url, show_col_types = FALSE)
  latest_model <- files %>%
    dplyr::filter(filename == 'Model.csv') %>%
    slice_head(n = 1) %>%
    pull(url)
  latest_CRISPRGeneEffect <- files %>%
    dplyr::filter(filename == 'CRISPRGeneEffect.csv') %>%
    slice_head(n = 1) %>%
    pull(url)

  return(list(model = latest_model, CRISPRGeneEffect = latest_CRISPRGeneEffect))
}

#test <- get_latest_depmap_release()


#' Helper Function 2 - get latest DepMap Model.csv as a df.
#'
#' This function performs a specific task. Getting the latest DepMap Model.csv as a df.
#'
#' @param file_url Url of Model.csv, the ouput of get_latest_depmap_release() $model
#' @param string A string of the model of interest.
#' @return Df of the latest DepMap release of Model.csv.
#' @export
get_lineage <- function(file_url, string) {
  file_url <- as.character(file_url)
  lineage_file <- readr::read_csv(file_url, col_names = TRUE, col_types = cols(.default = "c")) #, show_col_types = FALSE)
  lineage <- lineage_file %>%
    dplyr::select(`ModelID`, `OncotreeLineage`) %>%
    dplyr::filter(`OncotreeLineage` %in% string) %>%
    dplyr::rename(cell_line_id = `ModelID`)
  return(lineage)
}

#test1 <- get_lineage(test$model, "CNS/Brain")



#' Helper Function 3 - get model ids of the models of interest.
#'
#' This function performs a specific task. Getting the model ids of the lineage decided by the input.
#'
#' @param vector_of_OncotreeLineage  A vector of one or more model lineages. These are the available lineages:
#' Ovary/Fallopian Tube, Myeloid, Bowel, Skin, Bladder/Urinary Tract, Lung, Kidney
#' Breast, Lymphoid, Pancreas, CNS/Brain, Soft Tissue, Bone, Fibroblast
#' Esophagus/Stomach, Thyroid, Peripheral Nervous System, Pleura, Prostate
#' Biliary Tract, Head and Neck, Uterus, Ampulla of Vater, Liver, Cervix, Eye
#' Vulva/Vagina, Adrenal Gland, Testis, Other, Normal, NA, Muscle, Embryonal, Hair
#' @return Df of the model ids.
#' @export
get_cancer_cell_line_models <- function(vector_of_OncotreeLineage, model_url) {
  # Access DepMap Model.csv

  if (length(vector_of_OncotreeLineage) == 1) {
    return(get_lineage(model_url, vector_of_OncotreeLineage))
  }

  results <- data.frame()

  for (string in vector_of_OncotreeLineage) {
    result <- get_lineage(model_url, vector_of_OncotreeLineage)
    results <- bind_rows(results, result)
  }

  return(results)
}

#vector1 <- c("Myeloid")
#vector2 <- c("Myeloid", "CNS/Brain")
#test2 <- get_cancer_cell_line_models(vector1, test$model)
#test3 <- get_cancer_cell_line_models(vector2, test$model)

#model_test <- read_csv(test$model)




#' Main Function - get the cancer cell line essentiality of the models of interest.
#'
#' This is the main function. Getting cell line essentiality scores of the models of interest.
#'
#' @param genes A df with hgnc_ids column
#' @param models A df of the model ids of interest, obtained from get_cancer_cell_line_models()
#' @param CRISPRGeneEffects The url to the CRISPRGeneEffects file, obtained from get_latest_depmap_release()
#' @param get_mean The option to get the means of the cell lines instead of the raw essentiality scores.
#' Default is FALSE.
#' @return Df of the essentiality scores. Either the raw scores or the mean, depending on get_mean.
#' @export
get_cancer_cell_line_essentiality <- function(genes, models, CRISPRGeneEffects, get_mean = FALSE) {

  Sys.setenv("VROOM_CONNECTION_SIZE" = "2097152")

  #protein_coding_genes <- get_protein_coding_genes()

  #head(CRISPRGeneEffects)
  #print('CRISPRGeneEffects_ids')
  CRISPRGeneEffect_ids <- readr::read_csv(CRISPRGeneEffects, col_names = TRUE) %>%
    dplyr::rename(cell_line_id = '...1')

  #print('temp_CGE')
  temp_CGE <- CRISPRGeneEffect_ids %>%
    inner_join(models, by = 'cell_line_id') %>%
    dplyr::select(!OncotreeLineage)

  # Removing extra characters from the gene names
  colnames(temp_CGE) <- str_split(names(temp_CGE), "\\ ", simplify=T)[,1]

  # Omitting the first column
  temp_CGE = temp_CGE[,-1]

  # Transposing the df
  t_depmap <- t(temp_CGE)

  # Adding row names as column 1
  depmap <- tibble::rownames_to_column(data.frame(t_depmap))
  colnames(depmap)[1] <- 'symbol'

  #print('depmap_hgnc')
  # Using hgnc_checker
  depmap_hgnc <- hgnc_checker(depmap$symbol) #, protein_coding_genes)

  #print('depmap')
  # Adding HGNC IDs to depmap df
  depmap <- left_join(depmap, depmap_hgnc, join_by('symbol' == 'gene_symbol'))
  depmap <- depmap %>% dplyr::relocate(hgnc_id)

  # Remove last column (with the symbol type)
  depmap <- dplyr::select(depmap, -last_col())

  print('depmap 2')
  # Get depmap annotations for only genes of interest
  depmap <- genes %>%
    left_join(depmap, by = join_by(hgnc_id))

  #depmap[, 3:ncol(depmap)] <- lapply(depmap[, 3:ncol(depmap)], as.numeric)

  depmap[, 3:ncol(depmap)] <- lapply(
    depmap[, 3:ncol(depmap)],
    function(x) {
      # Attempt conversion
      y <- as.numeric(x)

      # Identify positions where `y` is NA but `x` was not originally NA
      idx <- which(is.na(y) & !is.na(x))

      # If any problematic rows exist, print them
      if (length(idx) > 0) {
        cat("These values are being coerced to NA:\n")
        print(data.frame(row = idx, original_value = x[idx]))
      }

      # Return the coerced vector
      y
    }
  )

  #print(depmap)

  if (get_mean) {
    #print('depmap 3.1')
    depman_mean <- depmap %>%
      dplyr::mutate(cell_lines_mean = rowMeans(dplyr::select(., 3:ncol(depmap)), na.rm = TRUE)) %>%
      dplyr::select(hgnc_id, cell_lines_mean)
    return(depman_mean)

  } else {
    return(depmap)
  }

}

#vector1 <- c("Myeloid")
#test <- get_latest_depmap_release()

#models1 <- get_cancer_cell_line_models(vector1, test$model)
#source('./R/get_protein_coding_genes.R')
#genes <- get_protein_coding_genes()
#genes50 <- head(genes, 50) %>%
 # dplyr::select(hgnc_id)


#test4 <- get_cancer_cell_line_essentiality(genes50, models1, test$CRISPRGeneEffect, get_mean = FALSE)
#test5 <- get_cancer_cell_line_essentiality(genes50, models1, test$CRISPRGeneEffect, get_mean = TRUE)
