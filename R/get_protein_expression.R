
#' Helper Function 1.1 - download protein expression
#'
#'
#'
#' @return A dataframe with protein expression from protein atlas
#' @export
download_normal_prot_expr <- function() {
  # Download and unzip the file, then read the TSV directly
  url <- "https://www.proteinatlas.org/download/tsv/normal_ihc_data.tsv.zip"
  temp <- tempfile()

  # Download the zip file
  download.file(url, temp)

  # Unzip and read the TSV file into a data frame
  df <- read_tsv(unz(temp, "normal_ihc_data.tsv"))

  return(df)
}


#' Helper Function 1.2 - protein expression, filter and remove columns
#'
#'
#'
#' @return A dataframe with protein expression mapped to hgnc_id, filtered columns
#' @export
filter_and_remove_columns <- function(data, cols_to_remove) {
  # Filter rows based on Reliability column
  filtered_data <- data[data$Reliability %in% c("Enhanced", "Approved", "Supported"), ]

  # Check if the specified columns exist in the data
  cols_to_remove <- intersect(cols_to_remove, colnames(filtered_data))

  # Remove the specified columns
  filtered_data <- filtered_data[, !colnames(filtered_data) %in% cols_to_remove]

  return(filtered_data)
}


#' Helper Function 1.3 - protein expression, select brain tissues
#'
#'
#'
#' @return A dataframe with protein expression mapped to hgnc_id, brain tissues only
#' @export
select_brain_tissues <- function(data) {
  # Select rows where the IHC tissue name contains the word 'brain'
  brain_data <- data %>%
    dplyr::filter(Tissue %in% c("Caudate", "Cerebellum", "Cerebral cortex",
                                "Choroid plexus", "Hippocampus", "Hypothalamus",
                                "Substantia nigra"))
  return(brain_data)
}


#' Helper Function 1.4 - protein expression, creating binary features
#'
#'
#'
#' @return A dataframe with protein expression mapped to hgnc_id, brain tissues only
#' @export
create_category_features <- function(data) {
  # Define the target categories
  target_categories <- c("Caudate", "Cerebellum", "Cerebral cortex",
                         "Choroid plexus", "Hippocampus", "Hypothalamus",
                         "Substantia nigra")

  # Helper function: convert rating text to numeric value
  rating_to_numeric <- function(rating) {
    if (is.na(rating)) {
      return(0)
    } else {
      # Normalize to lower case for consistency
      rating_lower <- tolower(rating)
      if (rating_lower == "low") {
        return(1)
      } else if (rating_lower == "medium") {
        return(2)
      } else if (rating_lower == "high") {
        return(3)
      } else {
        # Return NA if the rating is unrecognized
        return(NA_integer_)
      }
    }
  }

  # Loop over each target category
  for (cat in target_categories) {
    # Use the category name directly for the new column, replacing spaces with underscores
    new_col_name <- gsub(" ", "_", cat)

    # Initialize the new column with NA_integer_ for all rows
    data[[new_col_name]] <- NA_integer_

    # Identify rows where the second column equals the current category.
    # Adjust this if the column indexes differ from your expectation.
    is_target <- data[[2]] == cat

    # Convert the ratings in column 4 for these rows and assign them to the new feature.
    data[[new_col_name]][is_target] <- sapply(data[[4]][is_target], rating_to_numeric)
  }

  return(data)
}


#' Helper Function 1.5 - protein expression, tissue cell type scores
#'
#'
#'
#' @return A dataframe with protein expression mapped to hgnc_id, tissue cell type scoers
#' @export
create_tissue_celltype_level_scores <- function(data) {
  # Step 1: Clean the "Cell type" column.
  # Keep only the text before the first "-" and trim any extra whitespace.
  data$cleaned_celltype <- trimws(gsub("\\s*-.*$", "", data[["Cell type"]]))

  # Step 2: Create a combined category from Tissue and cleaned Cell type.
  data$combined_cat <- paste(data[["Tissue"]], data$cleaned_celltype, sep = "_")

  # Step 3: Create a helper function to convert level strings to numeric scores (0–3).
  rating_to_numeric <- function(rating) {
    # If rating is NA, return 0
    if (is.na(rating)) return(0)

    # Convert to lower-case to standardize.
    rating_lower <- tolower(rating)
    if (rating_lower == "low") {
      return(1)
    } else if (rating_lower == "medium") {
      return(2)
    } else if (rating_lower == "high") {
      return(3)
    } else {
      # For any unexpected value, return NA (or you could choose a default)
      return(NA_integer_)
    }
  }

  # Apply the conversion to the "Level" column (adjust the column name if needed).
  data$level_numeric <- sapply(data[["Level"]], rating_to_numeric)

  # Step 4: Create numeric dummy columns for each unique combination.
  # Instead of binary 0/1, assign the numeric level if the row matches the combination, else 0.
  unique_combinations <- unique(data$combined_cat)

  # Helper function to make valid column names (replace spaces or hyphens with underscores)
  make_valid_colname <- function(x) {
    make.names(gsub("[[:space:]-]+", "_", x))
  }

  for (comb in unique_combinations) {
    new_col_name <- make_valid_colname(comb)
    # For rows matching this combined category, use the level_numeric score; otherwise use 0.
    data[[new_col_name]] <- ifelse(data$combined_cat == comb, data$level_numeric, 0)
  }

  # Optional: Remove temporary columns if you no longer need them.
  data$cleaned_celltype <- NULL
  data$combined_cat <- NULL
  data$level_numeric <- NULL

  return(data)
}


#' Helper Function 1.6 - protein expression, distinct rows
#'
#'
#'
#' @return A dataframe with protein expression mapped to hgnc_id, tissue cell type scores, distinct rows
#' @export
reduce_to_one_row <- function(data, id_col = "hgnc_id") {
  # A custom function that returns the maximum of x, or NA if all are NA.
  max_na <- function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    } else {
      return(max(x, na.rm = TRUE))
    }
  }

  reduced_data <- data %>%
    group_by(across(all_of(id_col))) %>%
    summarise(across(where(is.numeric), max_na), .groups = "drop")

  names(reduced_data) <- tolower(names(reduced_data))

  # Replace any NA values with 0
  reduced_data[is.na(reduced_data)] <- 0

  return(reduced_data)
}


#' Helper Function 1.7 - protein expression, final df
#'
#'
#'
#' @return A dataframe with protein expression mapped to hgnc_id, final df
#' @export
final_pcg_protein_expression <- function(df1_pcg, df2_filtered_removed_columns, df3_reduced_one_row) {

  df1 <- df1_pcg %>% dplyr::select(hgnc_id) %>% distinct()

  # Get unique IDs from df2
  df2_ids <- unique(df2_filtered_removed_columns$hgnc_id)

  # Perform a left join of df3 to df1
  merged_df <- df1 %>%
    left_join(df3_reduced_one_row, by = "hgnc_id") %>%
    # For every df3-originating column, if the current row’s ID is in df2
    # and its value is NA, then replace it with 0.
    mutate(across(
      .cols = all_of(setdiff(names(df3_reduced_one_row), "hgnc_id")),
      .fns = ~ if_else(hgnc_id %in% df2_ids & is.na(.x), 0, .x)
    ))

  return(merged_df)
}


#' Main Function 1 - get protein expression from HPA
#'
#'
#'
#' @return A dataframe with protein expression from HPA mapped to hgnc_id
#' @export
get_protein_expression_hpa <- function(protein_coding_genes, brain_tissues_only = TRUE) {

  protein_expression <- download_normal_prot_expr()

  protein_expression <- protein_expression %>%
    dplyr::filter(Level != "Not detected" & Level != 'Not representative') %>%
    dplyr::select(-'Gene name') %>%
    dplyr::rename('ensembl_gene_id' = 'Gene')

  hgnc_protein_expression <- protein_coding_genes %>%
    dplyr::left_join(protein_expression, by = 'ensembl_gene_id')

  hgnc_protein_expression <- hgnc_protein_expression %>%
    dplyr::select(-'ensembl_gene_id')

  if (brain_tissues_only) {
    protexp2 <- filter_and_remove_columns(hgnc_protein_expression, c('IHC tissue name',
                                                     'Gene name', 'Reliability'))
    protexp3 <- select_brain_tissues(protexp2)
    protexp4 <- create_category_features(protexp3)
    protexp5 <- create_tissue_celltype_level_scores(protexp4)
    protexp6 <- reduce_to_one_row(protexp5)
    protexp7 <- final_pcg_protein_expression(hgnc_protein_expression, protexp2, protexp6)

    return(protexp7)
  } else {
    return(hgnc_protein_expression)
  }
}


#' Helper Function 2.1 - clean protein expression raw file
#'
#'
#' @param all_protein_data description
#' @return A clean dataframe with protein expression from ProteomicsDB
#' @export
clean_protein_expression_raw_file <- function(all_protein_data) {

  #load('./data/complete_protein_expression_data.RData')

  all_protein_data_clean <- all_protein_data %>%
    dplyr::select(hgnc_id, uniprot_ids, TISSUE_ID, TISSUE_NAME,
                  UNNORMALIZED_INTENSITY, NORMALIZED_INTENSITY,
                  MIN_NORMALIZED_INTENSITY, MAX_NORMALIZED_INTENSITY,
                  SAMPLES)

  return(all_protein_data_clean)
}

#' Helper Function 2.2 - show all tissues in protein expression raw file
#'
#'
#'
#' @return A clean dataframe with all tissues in protein expression from ProteomicsDB
#' @export
protein_expression_pdb_tissues <- function(all_protein_data) {

  tissues <- all_protein_data %>%
    dplyr::select(TISSUE_ID, TISSUE_NAME) %>%
    distinct()

  return(tissues)
}


#' Main Function 2 - get protein expression from ProteomicsDB
#'
#'
#'
#' @param protein_coding_genes A dataframe with protein coding genes
#' @param specific_tissue A string with the name of the tissue to filter for
#' @param zeros A boolean indicating whether to fill missing values with 0
#' @return A dataframe with protein expression from HPA mapped to hgnc_id
#' @export
get_protein_expression_pdb <- function(protein_coding_genes, specific_tissue = NULL,
                                       zeros = TRUE, all_protein_data) {

  pcg <- protein_coding_genes %>%
    dplyr::select(hgnc_id)

  protein_expression_clean <- clean_protein_expression_raw_file(all_protein_data)

  if (is.null(specific_tissue)) {
    # Get all protein expression data
    protein_expression <- protein_expression_clean
  } else {
    # Filter for the specified tissue
    protein_expression <- protein_expression_clean %>%
      dplyr::filter(TISSUE_NAME == specific_tissue)
  }

  if (zeros) {
    # Fill missing values with 0

    pcg2 <- pcg %>%
      dplyr::left_join(protein_expression, by = 'hgnc_id',
                       suffix = c("", "_y"))

    protein_expression2 <- pcg2 %>%
      dplyr::mutate(across(c(UNNORMALIZED_INTENSITY, NORMALIZED_INTENSITY,
                             MIN_NORMALIZED_INTENSITY, MAX_NORMALIZED_INTENSITY),
                           ~ replace_na(.x, 0)))
  }
  return(protein_expression2)
}
