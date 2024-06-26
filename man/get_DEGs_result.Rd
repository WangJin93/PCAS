% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_DEGs_result.R
\name{get_DEGs_result}
\alias{get_DEGs_result}
\title{Get DEGs results}
\usage{
get_DEGs_result(dataset = "LUAD_CPTAC_protein", method = "t.test")
}
\arguments{
\item{dataset}{Use dataset$Abbre to get all tables.}

\item{method}{"limma" or "t.test".}
}
\description{
Get the results of different expression analysis between tumor and normal samples in CPTAC datasets.
}
\examples{
\dontrun{
results <- get_DEGs_result(dataset = "LUAD_CPTAC_protein", method = "limma")
}
}
