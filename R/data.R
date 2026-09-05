#' Dataset directory of the PCAS server
#'
#' Table of contents of the CPTAC datasets served by the PCAS API: molecular
#' type, study, sample counts and the abbreviation (\code{Abbre}) used to query
#' the data.
#'
#' @format A data frame with 56 rows and 8 columns:
#' \describe{
#'   \item{Data.type}{Molecular type: Proteome, Phosphoproteome or Transcriptome}
#'   \item{Dataset}{Study / dataset description}
#'   \item{Normal}{Number of normal samples (NA = the cohort has no normal samples)}
#'   \item{Tumor}{Number of tumour samples}
#'   \item{Abbre}{Dataset abbreviation used in the query functions}
#'   \item{Study.ID}{Study identifier}
#'   \item{Disease.Type}{Disease type(s)}
#'   \item{Primary.Site}{Primary tumour site}
#' }
#' @source PCAS server (\url{https://www.jingege.wang/bioinformatics/PCAS/})
#' @name dataset_info
#' @docType data
#' @usage data("dataset_info")
NULL

#' mRNA identifier map (gene symbols to probe ids)
#'
#' Maps Ensembl mRNA ids returned by the API to gene symbols.
#'
#' @format A data frame with 60774 rows and 4 columns:
#' \describe{
#'   \item{mRNAs}{mRNA/probe id}
#'   \item{row_names}{mRNA id used for the API query}
#'   \item{Symbol}{Gene symbol}
#'   \item{gene_type}{Gene biotype, e.g. protein_coding}
#' }
#' @source PCAS server
#' @name idmap_RNA
#' @docType data
#' @usage data("idmap_RNA")
NULL

#' Protein / phosphosite identifier map
#'
#' Maps protein ids and phosphorylation-site ids (e.g.
#' \code{NP_000537.3:s315}, combined ids such as \code{NP_000025.1:s218y223t227})
#' to gene symbols.
#'
#' @format A data frame with 172404 rows and 2 columns:
#' \describe{
#'   \item{row_names}{Protein accession, or phosphosite accession:site(s)}
#'   \item{Symbol}{Gene symbol}
#' }
#' @source PCAS server
#' @name idmap_protein
#' @docType data
#' @usage data("idmap_protein")
NULL

#' Drug sensitivity matrix
#'
#' Drug-sensitivity measurements (rows = CPTAC sample ids, columns = drugs
#' named \code{<Drug>_<PubChem-like id>}).
#'
#' @format A data frame with 1553 rows and 198 columns. Row names are sample
#'   ids; columns are drugs.
#' @source PCAS server / CPTAC pharmacoproteomics
#' @name drug_CPTAC
#' @docType data
#' @usage data("drug_CPTAC")
NULL

#' Drug annotation table
#'
#' @format A data frame with 198 rows and 6 columns:
#' \describe{
#'   \item{ID}{Drug id used in the column names of \code{drug_CPTAC}}
#'   \item{Name}{Drug name}
#'   \item{Synonyms}{Drug synonyms}
#'   \item{Targets}{Drug targets}
#'   \item{Target.pathway}{Signalling pathway of the targets}
#'   \item{PubCHEM}{PubChem identifier}
#' }
#' @source PCAS server
#' @name drug_info
#' @docType data
#' @usage data("drug_info")
NULL

#' Immune infiltration table
#'
#' Immune-cell infiltration scores of CPTAC samples, computed with several
#' algorithms. Column names end with the algorithm, e.g.
#' \code{Bcells_EPIC} or \code{T_cells_TIMER}.
#'
#' @format A data frame with 1636 rows (samples; column \code{ID}) and 137
#'   infiltration-score columns.
#' @source PCAS server / TIL algorithms
#' @name TIL_CPTAC
#' @docType data
#' @usage data("TIL_CPTAC")
NULL

#' Map from infiltration cell types to algorithms
#'
#' @format A data frame with 137 rows and 2 columns:
#' \describe{
#'   \item{cell_type}{Cell-type column name in \code{TIL_CPTAC}}
#'   \item{algorithm}{Algorithm that produced the score}
#' }
#' @source PCAS server
#' @name TIL_map
#' @docType data
#' @usage data("TIL_map")
NULL

#' Gene symbol to UniProt entry map
#'
#' @format A data frame with 15162 rows and 4 columns:
#' \describe{
#'   \item{From}{UniProt/RefSeq protein accession}
#'   \item{Entry}{UniProt entry name/accession}
#'   \item{Reviewed}{"reviewed" when the entry is Swiss-Prot reviewed}
#'   \item{Symbol}{Gene symbol}
#' }
#' @source UniProt
#' @name uniport_map
#' @docType data
#' @usage data("uniport_map")
NULL

#' Phosphosite ids available per dataset
#'
#' A named list: each element is a dataset abbreviation, the values are the
#' phosphorylation-site ids measured in that dataset (used to populate the
#' site dropdowns of the Shiny app).
#'
#' @format A named list.
#' @source PCAS server
#' @name phoso_id
#' @docType data
#' @usage data("phoso_id")
NULL

#' Gene/protein symbols available on the PCAS server
#'
#' A named list of identifier sets used to populate gene dropdowns.
#'
#' @format A named list.
#' @source PCAS server
#' @name ID_list_pro
#' @docType data
#' @usage data("ID_list_pro")
NULL
