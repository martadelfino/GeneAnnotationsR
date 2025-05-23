% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_alphamissense.R
\name{map_alphamissense_ids}
\alias{map_alphamissense_ids}
\title{Helper Function 2 - Maps alphamissense ENST to ENSG and HGNC IDs}
\usage{
map_alphamissense_ids(alphamissense)
}
\arguments{
\item{alphamissense}{A df that is output of read_alphamissense_gene_hg38().
It is a df with two columns 1. ENST and 2. mean_am_pathogenicity.}
}
\value{
Dataframe of alphamissense gene hg38 scores with their mapped ENSG,
HGNC IDs.
}
\description{
This function performs a specific task. Uses biomart to map alphamissense
ENST ids to ENSG and HGNC IDs.
}
