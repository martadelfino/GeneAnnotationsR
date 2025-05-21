# Function to get protein expression data for a dataframe of IDs and UniProt IDs
get_protein_expression <- function(protein_df,
                                   ms_level = 1,
                                   tissue_category = "tissue;fluid",
                                   scope = 1,
                                   group_by_tissue = 1,
                                   calculation_method = 0) {

  # Check if input dataframe has required columns
  if (!all(c("hgnc_id", "uniprot_ids") %in% colnames(protein_df))) {
    stop("Input dataframe must contain 'hgnc_id' and 'uniprot_ids' columns")
  }

  # Load required packages
  if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr")
  if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

  library(httr)
  library(jsonlite)
  library(dplyr)

  # Initialize results dataframe
  all_results <- data.frame()

  # Base URL for the ProteomicsDB API
  base_url <- "https://www.proteomicsdb.org/proteomicsdb/logic/api/proteinexpression.xsodata/InputParams"

  # Loop through each row in the dataframe
  for (i in 1:nrow(protein_df)) {
    hgnc_id <- protein_df$hgnc_id[i]
    uniprot_ids <- protein_df$uniprot_ids[i]

    # Print progress
    cat(sprintf("Retrieving data for HGNC ID: %s, UniProt: %s (%d of %d)\n",
                hgnc_id, uniprot_ids, i, nrow(protein_df)))

    # Construct API parameters
    params <- list(
      PROTEINFILTER = uniprot_ids,
      MS_LEVEL = ms_level,
      TISSUE_ID_SELECTION = '',
      TISSUE_CATEGORY_SELECTION = tissue_category,
      SCOPE_SELECTION = scope,
      GROUP_BY_TISSUE = group_by_tissue,
      CALCULATION_METHOD = calculation_method,
      EXP_ID = -1
    )

    # Construct the API URL
    url_params <- paste0(
      "(PROTEINFILTER='", params$PROTEINFILTER,
      "',MS_LEVEL=", params$MS_LEVEL,
      ",TISSUE_ID_SELECTION='", params$TISSUE_ID_SELECTION,
      "',TISSUE_CATEGORY_SELECTION='", params$TISSUE_CATEGORY_SELECTION,
      "',SCOPE_SELECTION=", params$SCOPE_SELECTION,
      ",GROUP_BY_TISSUE=", params$GROUP_BY_TISSUE,
      ",CALCULATION_METHOD=", params$CALCULATION_METHOD,
      ",EXP_ID=", params$EXP_ID, ")"
    )

    api_url <- paste0(
      base_url, url_params,
      "/Results?$select=UNIQUE_IDENTIFIER,TISSUE_ID,TISSUE_NAME,TISSUE_SAP_SYNONYM,SAMPLE_ID,SAMPLE_NAME,",
      "AFFINITY_PURIFICATION,EXPERIMENT_ID,EXPERIMENT_NAME,EXPERIMENT_SCOPE,EXPERIMENT_SCOPE_NAME,",
      "PROJECT_ID,PROJECT_NAME,PROJECT_STATUS,UNNORMALIZED_INTENSITY,NORMALIZED_INTENSITY,",
      "MIN_NORMALIZED_INTENSITY,MAX_NORMALIZED_INTENSITY,SAMPLES&$format=json"
    )

    # Make the API request
    tryCatch({
      response <- GET(api_url)

      # Check if request was successful
      if (http_status(response)$category == "Success") {
        # Parse JSON response
        json_data <- fromJSON(content(response, "text", encoding = "UTF-8"))

        # Extract results
        if ("d" %in% names(json_data) && "results" %in% names(json_data$d)) {
          results <- json_data$d$results

          # Add both ID and UniProt ID as columns
          results$hgnc_id <- hgnc_id
          results$uniprot_ids <- uniprot_ids

          # Append to all results
          all_results <- bind_rows(all_results, results)
        } else {
          warning(paste("No results found for HGNC ID:", hgnc_id, "UniProt:", uniprot_ids))
        }
      } else {
        warning(paste("Failed to retrieve data for HGNC ID:", hgnc_id, "UniProt:", uniprot_ids,
                      "- Status:", http_status(response)$message))
      }
    }, error = function(e) {
      warning(paste("Error retrieving data for hgnc ID:", hgnc_id, "UniProt:", uniprot_ids,
                    "-", e$message))
    })

    # Add a small delay to avoid overloading the API
    Sys.sleep(1)
  }

  return(all_results)
}



library(tidyverse)
library(GeneAnnotationsR)



protein_coding_genes <- get_protein_coding_genes() %>%
  dplyr::select(hgnc_id, uniprot_ids) %>%
  dplyr::filter(uniprot_ids != "")


# Get protein expression data
#protein_data <- get_protein_expression(
# protein_df = protein_coding_genes,
#ms_level = 1,
#  tissue_category = "tissue;fluid",
# scope = 1,
#group_by_tissue = 1,
#  calculation_method = 0
#)



# Display the first few rows of the results
#if (nrow(protein_data) > 0) {
# head(protein_data)

# Save results to CSV
#  write.csv(protein_data, file = "data/data-raw/protein_expression_data.csv", row.names = FALSE)
# cat("Results saved to protein_expression_data.csv\n")
#} else {
# cat("No data retrieved for the specified proteins\n")
#}






###############################################################################
# Simple R script to fetch tissue data from ProteomicsDB
###############################################################################

# Load required libraries
library(httr)
library(jsonlite)

# Get tissue data from ProteomicsDB API
get_tissue_data <- function() {
  # API URL
  url <- "https://www.proteomicsdb.org/proteomicsdb/logic/api/tissuelist.xsodata/CA_AVAILABLEBIOLOGICALSOURCES_API"

  # Create API query
  query <- list(
    "$select" = "TISSUE_ID,TISSUE_NAME,TISSUE_GROUP_NAME,TISSUE_CATEGORY,SCOPE_ID,SCOPE_NAME,QUANTIFICATION_METHOD_ID,QUANTIFICATION_METHOD_NAME,MS_LEVEL,TAXCODE",
    "$format" = "json"
  )

  # Make the request
  response <- GET(url, query = query)

  # Check if successful
  if (status_code(response) == 200) {
    # Parse JSON response
    data <- fromJSON(content(response, "text", encoding = "UTF-8"))

    # Extract results
    tissues <- data$d$results

    # Print summary
    cat("Retrieved", nrow(tissues), "tissues from ProteomicsDB\n")

    return(tissues)
  } else {
    cat("Error fetching data:", http_status(response)$message, "\n")
    return(NULL)
  }
}

# Main script
cat("Fetching tissue data from ProteomicsDB...\n")

# Get the data
tissues <- get_tissue_data()

unique_tissues <- data.frame(unique(tissues$TISSUE_GROUP_NAME))





# Save to CSV if data was retrieved
if (!is.null(tissues)) {
  write.csv(tissues, "proteomicsdb_tissues.csv", row.names = FALSE)
  cat("Data saved to 'proteomicsdb_tissues.csv'\n")

  # Show first few rows
  cat("\nFirst 5 tissue entries:\n")
  print(head(tissues, 5))
}



###############################################################################
# Get protein expression data for all protein coding genes
###############################################################################


all_protein_data <- data.frame()

# Set batch size
batch_size <- 1000
total_rows <- nrow(protein_coding_genes)
num_batches <- ceiling(total_rows / batch_size)

# Process in batches
for (i in 1:num_batches) {
  # Get row indices for current batch
  start_row <- ((i-1) * batch_size) + 1
  end_row <- min(i * batch_size, total_rows)

  cat(sprintf("Processing batch %d of %d (rows %d to %d)\n", i, num_batches, start_row, end_row))

  # Extract current batch
  current_batch <- protein_coding_genes[start_row:end_row, ]

  # Get protein expression data for current batch
  batch_results <- get_protein_expression(
    protein_df = current_batch,
    ms_level = 1,
    tissue_category = "tissue;fluid",
    scope = 1,
    group_by_tissue = 1,
    calculation_method = 0
  )

  # Combine with previous results
  all_protein_data <- bind_rows(all_protein_data, batch_results)

  # Save progress
  write.csv(batch_results, file = paste0("batch_", i, "_protein_data.csv"), row.names = FALSE)
}

# Save final complete dataset
write.csv(all_protein_data, file = "complete_protein_expression_data.csv", row.names = FALSE)

save(all_protein_data, file = "complete_protein_expression_data.RData")

unique_tissues <- data.frame(unique(all_protein_data$TISSUE_NAME))



###############################################################################
# Checks for heart and brain tissues
###############################################################################

heart_tissues <- all_protein_data %>%
  filter(TISSUE_NAME == "heart" | TISSUE_NAME == "smooth muscle" | TISSUE_NAME == "lung" | TISSUE_NAME == "blood" | TISSUE_NAME == "embryonic stem cell") %>%
  select(hgnc_id, uniprot_ids, TISSUE_NAME, NORMALIZED_INTENSITY) # %>%

unique_heart_tissues <- data.frame(unique(heart_tissues$hgnc_id))

brain_tissues <- all_protein_data %>%
  filter(TISSUE_NAME == "brain" | TISSUE_NAME == "cerebral cortex" | TISSUE_NAME == "spinal cord" | TISSUE_NAME  == "cerebrospinal fluid" | TISSUE_NAME == "embryonic stem cell") %>%
  select(hgnc_id, uniprot_ids, TISSUE_NAME, NORMALIZED_INTENSITY) # %>%

unique_brain_tissues <- data.frame(unique(brain_tissues$hgnc_id))


heart_tissues2 <- all_protein_data %>%
  filter(TISSUE_NAME == "heart" | TISSUE_NAME == "smooth muscle" | TISSUE_NAME == "lung" | TISSUE_NAME == "blood") %>%
  select(hgnc_id, uniprot_ids, TISSUE_NAME, NORMALIZED_INTENSITY) # %>%

unique_heart_tissues2 <- data.frame(unique(heart_tissues2$hgnc_id))


