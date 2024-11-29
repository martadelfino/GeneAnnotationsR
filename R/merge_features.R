


merge_features <- function(df_list) {

  # Use Reduce to iteratively merge each data frame in the list
  merged_df <- Reduce(function(x, y) merge(x, y, by = "hgnc_id", all = TRUE), df_list)

  return(merged_df)
}
