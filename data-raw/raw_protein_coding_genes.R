# Script to save the raw file of protein_coding_genes from EBI


protein_coding_genes <- data.table::fread(
  "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/locus_types/gene_with_protein_product.txt",
  sep = "\t",     # Define tab-separated format
  header = TRUE,  # File contains headers
  select = NULL   # You can specify column names to load, if needed
)

# Specify your folder and file name
folder_path <- "data-raw/"
file_name <- "protein_coding_genes.csv.gz"

# Construct the full file path
file_path <- file.path(folder_path, file_name)

# Open a gzfile connection to that path
con <- gzfile(file_path, "w")

# Write the data frame 'df' to the gzipped file without row names
write.csv(protein_coding_genes, con, row.names = FALSE)

# Close the connection
close(con)
