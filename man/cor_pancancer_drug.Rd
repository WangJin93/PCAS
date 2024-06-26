% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cor_pancancer_drug.R
\name{cor_pancancer_drug}
\alias{cor_pancancer_drug}
\title{Correlation analysis of drug sensitivity}
\usage{
cor_pancancer_drug(
  df,
  cor_method = "pearson",
  Target.pathway = c("Cell cycle")
)
}
\arguments{
\item{df}{The expression data of the target gene in multiple datasets, obtained by the get_expr_data() function.}

\item{cor_method}{Method for correlation analysis, default "pearson".}

\item{Target.pathway}{The signaling pathways of anti-tumor drug targets, default "Cell cycle". Use "drug_info"to get the detail infomation of these drugs.}
}
\description{
Calculate the correlation between target gene expression and anti-tumor drug sensitivity in multiple datasets.
}
\examples{
\dontrun{
dataset <- c("CCRCC_CPTAC_protein","GBM_CPTAC_protein","HNSCC_CPTAC_protein","LSCC_CPTAC_protein","LUAD_CPTAC_protein","PDAC_CPTAC_protein","UCEC_CPTAC2_protein","UCEC_CPTAC1_protein")
df <- get_expr_data(genes = "TNS1",datasets = dataset)
result <- cor_pancancer_drug(df,Target.pathway = c("Cell cycle"))
}
}
