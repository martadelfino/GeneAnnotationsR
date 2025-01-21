
#' Helper Function 1 - download subcellular locations from protein atlas
#'
#'
#'
#' @return A dataframe with subcellular locations from protein atlas
#' @export
download_subcell_location_file <- function() {
  # Download and unzip the file, then read the TSV directly
  url <- "https://www.proteinatlas.org/download/tsv/subcellular_location.tsv.zip"
  temp <- tempfile()

  # Download the zip file
  download.file(url, temp)

  # Unzip and read the TSV file into a data frame
  df <- read_tsv(unz(temp, "subcellular_location.tsv"))

  return(df)
}

#' Helper Function 2 - keep subcell columns of high support
#'
#'
#'
#' @return A dataframe with subcellular locations from protein atlas with high support
#' @export
high_support_locations <- function(subcellular_location) {
  # Filter for high support locations
  high_support <- subcellular_location %>%
    dplyr::select(hgnc_id, `Extracellular location`, Enhanced, Supported,
                  Approved) %>%
    return(high_support)
}


#' Helper Function 3 - subcellular locations remove parenthesis
#'
#'
#'
#' @return A dataframe with subcellular locations from protein atlas with parenthesis removed
#' @export
remove_parentheses <- function(data, col_name) {
  # Check if the column exists
  if (!col_name %in% colnames(data)) {
    stop(paste("Column", col_name, "not found in dataframe."))
  }

  # Remove parentheses and their contents
  data[[col_name]] <- gsub("\\(.*?\\)", "", data[[col_name]])

  # Trim any leading or trailing whitespace
  data[[col_name]] <- trimws(data[[col_name]])

  return(data)
}


#' Helper Function 4 - create binary features for subcell
#'
#'
#'
#' @return A dataframe with subcellular locations from protein atlas in binary format
#' @export
create_binary_features_subcell <- function(data, col_names, remove_parentheses_cols = NULL) {
  # Helper function to clean column names
  clean_name <- function(name) {
    name <- tolower(name)                # Convert to lowercase
    name <- gsub("\\s+", "_", name)      # Replace spaces with underscores
    name <- gsub("[^a-z0-9_]", "_", name) # Replace non-alphanumeric characters with underscores
    name <- gsub("_+", "_", name)        # Remove consecutive underscores
    name <- gsub("^_|_$", "", name)      # Remove leading/trailing underscores
    return(name)
  }

  # Remove parentheses from specified columns
  if (!is.null(remove_parentheses_cols)) {
    for (col_name in remove_parentheses_cols) {
      if (col_name %in% colnames(data)) {
        data <- remove_parentheses(data, col_name)
      } else {
        warning(paste("Column", col_name, "not found in dataframe. Skipping."))
      }
    }
  }

  # Loop through each column name in the provided vector
  for (col_name in col_names) {
    # Check if the column exists
    if (!col_name %in% colnames(data)) {
      warning(paste("Column", col_name, "not found in dataframe. Skipping."))
      next
    }

    # Extract unique strings from the column, filtering out empty and unwanted values
    unique_values <- unique(unlist(strsplit(data[[col_name]], ";")))
    unique_values <- trimws(unique_values) # Remove extra spaces
    unique_values <- unique_values[unique_values != ""] # Remove empty strings
    unique_values <- unique_values[!is.na(unique_values)] # Remove NA values
    unique_values <- unique_values[!tolower(unique_values) %in% c("na", "n/a", "n a")] # Remove "NA", "N/A", etc.

    # Create new columns for each unique value
    for (value in unique_values) {
      cleaned_value <- clean_name(value)  # Clean the value
      col_new_name <- paste(clean_name(col_name), cleaned_value, sep = "_")
      data[[col_new_name]] <- as.integer(grepl(paste0("\\b", value, "\\b"), data[[col_name]]))
    }

    # Optionally remove the original column
    data[[col_name]] <- NULL
  }

  return(data)
}


#' Helper Function 5 - merge prefixed features in subcell
#'
#'
#'
#' @return A dataframe with subcellular locations from protein atlas merged prefixes
#' @export
merge_prefixed_features <- function(data) {
  # Get all column names
  col_names <- colnames(data)

  # Identify columns with the desired prefixes
  pattern <- "^(enhanced|supported|approved)_(.*)"
  prefixed_cols <- grep(pattern, col_names, value = TRUE)

  # Create a list to group columns by their suffixes
  suffix_groups <- split(prefixed_cols, sub(pattern, "\\2", prefixed_cols))

  # Iterate over each suffix group to merge columns
  for (suffix in names(suffix_groups)) {
    group_cols <- suffix_groups[[suffix]]

    # Create the new merged column with the original suffix name
    data[[suffix]] <- rowSums(data[group_cols]) > 0
    data[[suffix]] <- as.integer(data[[suffix]])

    # Drop the original columns
    data[group_cols] <- NULL
  }

  return(data)
}


#' Main Function - get subcellular locations
#'
#' keeping all reliability levels.
#'
#' @return A dataframe with subcellular locations mapped to hgnc_id
#' @export
get_subcellular_location <- function(protein_coding_genes,
                                     high_support_and_binary = TRUE) {

  # Subcellular location file
  subcellular <- download_subcell_location_file()

  subcellular <- subcellular %>%
    dplyr::rename(ensembl_gene_id = Gene)

  hgnc_subcellular <- protein_coding_genes %>%
    dplyr::select(hgnc_id, ensembl_gene_id) %>%
    left_join(subcellular, by = 'ensembl_gene_id') %>%
    dplyr::select(-`Gene name`, -`ensembl_gene_id`, -`Reliability`, -`GO id`)

  if (high_support_and_binary) {
    high_support_subcell <- high_support_locations(hgnc_subcellular)
    high_support_binary_subcell <- create_binary_features_subcell(high_support_subcell,
                                                       c("Extracellular location",
                                                         "Enhanced", "Supported",
                                                         "Approved"))
    cleaned_high_support_binary_subcell <- merge_prefixed_features(high_support_binary_subcell)

    return(cleaned_high_support_binary_subcell)
    } else {
    return(hgnc_subcellular)
  }
}

