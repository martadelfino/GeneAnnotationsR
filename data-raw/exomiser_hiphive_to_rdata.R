
if (!require("tidyverse")) {
  install.packages("tidyverse")
  library(tidyverse)
}

# Define the path to your  file
file_path <- "data-raw/NDD_HiPhive_Prioritiser_mar25.genes.tsv"

# Read the file into a data frame
exomiser_pheno_score <- read_delim(file_path, delim = "\t")


save(exomiser_pheno_score, file = "data/NDD_Exomiser_Pheno_Score.RData")

