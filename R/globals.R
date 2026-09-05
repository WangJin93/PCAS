# Column names used inside ggplot2::aes() (non-standard evaluation) and the
# lazily-loaded package datasets referenced from function bodies. Declared here
# so R CMD check / codetools do not report them as undefined globals.
utils::globalVariables(c(
  "type", "value", "gene", "dataset", "label",
  "logFC", "neg_log10p", "change",
  "corr", "location", "text",
  "geneA", "geneB",
  # package datasets (LazyData)
  "dataset_info", "idmap_RNA", "idmap_protein", "drug_CPTAC", "drug_info",
  "TIL_CPTAC", "TIL_map", "uniport_map", "phoso_id", "ID_list_pro"
))
