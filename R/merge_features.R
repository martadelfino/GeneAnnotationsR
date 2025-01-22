
#' Main Function - merging features
#'
#' Merging data/features from multiple sources.
#'
#' @param df_list List of dataframes to merge.
#' @return A dataframe with hgnc_id and all the features.
#' @export
merge_features <- function(df_list) {

  # Use Reduce to iteratively merge each data frame in the list
  merged_df <- Reduce(function(x, y) merge(x, y, by = "hgnc_id", all = TRUE), df_list)

  return(merged_df)
}
