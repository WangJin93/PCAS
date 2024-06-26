% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_clinic_data.R
\name{merge_clinic_data}
\alias{merge_clinic_data}
\title{Merge clinic data}
\usage{
merge_clinic_data(cohort = "LUAD_APOLLO", data_input)
}
\arguments{
\item{cohort}{Data cohort, for example, "LUAD_APOLLO", "LUAD_CPTAC".}

\item{data_input}{Expression data obtained from get_expr_data() function.}
}
\description{
Get clinic data and merge it with expression data.
}
\examples{
\dontrun{
data_input <- get_expr_data("LUAD_APOLLO_mRNA","TP53")
results <- merge_clinic_data("LUAD_APOLLO",data_input)
}
}
