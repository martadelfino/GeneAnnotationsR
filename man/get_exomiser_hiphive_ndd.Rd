% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_exomiser_hiphive_ndd.R
\name{get_exomiser_hiphive_ndd}
\alias{get_exomiser_hiphive_ndd}
\title{Main function - get Exomiser hiphive NDD scores}
\usage{
get_exomiser_hiphive_ndd(protein_coding_genes)
}
\arguments{
\item{protein_coding_genes}{A df of entrez_id and hgnc_id.}
}
\value{
A df with hgnc_id and the NDD hiphive scores.
}
\description{
This is a main function. Read the RData file in the data directory. This
function only works for NDD genes.
}
