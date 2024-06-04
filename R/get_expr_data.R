#' @title Get CPTAC expression data
#' @description
#' Get the mRNA/protein expression data in CPTAC database.
#' @import tibble reshape2
#' @param datasets Dataset names, you can input one or multiple datasets. Use 'dataset$Abbre' to get all datasets.
#' @param genes Gene symbols, you can input one or multiple symbols.
#' @examples
#' \dontrun{
#' results <- get_expr_data(datasets = "LUAD_CPTAC_mRNA", genes = c("GAPDH","TNS1"))
#' results <- get_expr_data(datasets = c("LUAD_CPTAC_protein","LSCC_CPTAC_protein"), genes = "GAPDH")
#' results <- get_expr_data(datasets = c("CCRCC_CPTAC_mRNA","GBM_CPTAC_mRNA","HNSCC_CPTAC_mRNA","LSCC_CPTAC_mRNA","LUAD_CPTAC_mRNA","PDAC_CPTAC_mRNA","UCEC_CPTAC2_mRNA","UCEC_CPTAC1_mRNA"), genes = c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7","HAVCR2","LILRB4","PDCD1","BTLA"))
#' }
#' @export
#'
get_expr_data <- function(datasets=c("LUAD_CPTAC_protein","LSCC_CPTAC_protein"), genes= c("TP53","TNS1")) {
  if (length(genes)==0){
    return(NULL)
  }else{
    message("Querying data of identifier ", paste0(genes,collapse = ", "), " from datasets ", paste0(datasets,collapse = ", "))
    cptac_data<- data.frame()
    for (x in datasets) {
      if (stringr::str_detect(x,"mRNA")) {
        id <- idmap_RNA[which(idmap_RNA$Symbol %in% genes),]
        ids <- id$row_names
      }else{
        ids <- genes
      }
      data <- get_data(x,
                       "expression",
                       ids)

      if (is.null(nrow(data))){
        message(paste0(genes,collapse = ", "), " retrive no results in ", x)
        next()
      }else{
        if (stringr::str_detect(x,"mRNA")) {
          data <- merge(id[2:3],data,by="row_names")[-1]
          colnames(data)[1] <- "row_names"
        }
        row.names(data) <- NULL
        data<-data%>% tibble::column_to_rownames(var = "row_names")	%>% t() %>% as.data.frame()
        data$type <- lapply(row.names(data), function(x) strsplit(x,"_")[[1]][length(strsplit(x,"_")[[1]])])%>% as.character()
        data$ID <- stringr::str_remove(row.names(data),"_Tumor|_Normal|_Other")
        # data <- tibble::rownames_to_column(data,var = "ID")
        # data <- reshape2::melt(data, measure.vars = genes)
        # data$value <- as.numeric(data$value)
        data <- data[c("ID","type",colnames(data)[1:(ncol(data)-2)])]
        data$dataset <- x
        cptac_data <- plyr::rbind.fill(cptac_data,data)
        }
    }
    if (nrow(cptac_data) == 0){
      message("Retrive no data.")
      return(NULL)
    }
    if (length(genes)==1){
      cptac_data[,3] <- as.numeric(cptac_data[,3])
    }else{
      cptac_data[intersect(colnames(cptac_data),genes)] <- apply(cptac_data[intersect(colnames(cptac_data),genes)], 2, as.numeric)
    }
    cptac_data <- cptac_data[c("ID","type","dataset",intersect(colnames(cptac_data),genes))]
    return(cptac_data)
  }
}

