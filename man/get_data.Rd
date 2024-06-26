% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_data.R
\name{get_data}
\alias{get_data}
\title{Get CPTAC data}
\usage{
get_data(
  table = "LUAD_Academia_protein",
  action = "expression",
  genes = c("GAPDH", "TNS1")
)
}
\arguments{
\item{table}{For action = expression, use dataset$Abbre to get all tables; For action = clinic, remove _protein/_mRNA/_Phospho from dataset$Abbre.}

\item{action}{"expression", "degs" or "clinic".}

\item{gene}{Gene symbols, you can input one or multiple symbols.}
}
\description{
Get the CPTAC data by using the api. All results saved in MySQL database.
}
\examples{
\dontrun{
results <- get_data(table = "LUAD_Academia_protein", action = "expression", genes = c("GAPDH","TNS1"))
}
}
