% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mouse_viability.R
\name{get_mgi_viability}
\alias{get_mgi_viability}
\title{Helper Function 3 - get MGI viability.}
\usage{
get_mgi_viability(mouse_proteincoding_genes, lethal_terms)
}
\arguments{
\item{mouse_proteincoding_genes}{Output of get_mouse_protein_coding_genes()}

\item{lethal_terms}{lethal terms.}
}
\value{
Dataframe of MGI viability.
\describe{
\item{mgi_id}{MDI IF for gene}
\item{viability_mgi}{MGI viability}
}
}
\description{
This function performs a specific task. Getting MGI viability.
}
