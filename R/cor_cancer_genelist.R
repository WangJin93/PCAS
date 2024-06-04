#' @title Perform correlation analysis
#' @description
#' Perform correlation analysis of the mRNA/protein expression data in CPTAC database.
#' @import tibble reshape2
#' @param dataset1 Dataset name. Use 'dataset$Abbre' to get all datasets.
#' @param id1 Gene symbol, you can input one gene symbols.
#' @param dataset2 Dataset name. Use 'dataset$Abbre' to get all datasets.
#' @param id2 Gene symbols, you can input one or multiple symbols.
#' @param sample_type Sample type used for correlation analysis, default all types: c("Tumor", "Normal").
#' @param cor_method cor_method for correlation analysis, default "pearson".
#' @examples
#' \dontrun{
#' results <- correlation_analysis(dataset1 = "LUAD_CPTAC_protein",
#'                                 id1 = "STAT3",
#'                                 dataset2 = "LUAD_CPTAC_Phospho",
#'                                 id2 = c("NP_000499.1:s549","NP_000499.1:s551","TP503"),
#'                                 sample_type = c("Tumor","Normal"),
#'                                 cor_method = "pearson")
#' }
#' @export
#'
cor_cancer_genelist <- function(dataset1 = "LUAD_CPTAC_protein",
                              id1 = "STAT3",
                              dataset2 = "LUAD_CPTAC_mRNA",
                              id2 = c("TNS1","TP53"),
                              sample_type = c("Tumor","Normal"),
                              cor_method = "pearson"){
  data1 <- get_expr_data(dataset1, id1)
  data1$dataset<- data1$dataset %>% stringr::str_remove(.,"_protein|_Phospho|_mRNA")

  if (is.null(data1)){
    return(NULL)
  }

  data2 <- get_expr_data(dataset2, id2)
  data2$dataset<- data2$dataset %>% stringr::str_remove(.,"_protein|_Phospho|_mRNA")
  data_merge <- merge(data1[intersect(rownames(data1),rownames(data2)),],
                      data2[intersect(rownames(data1),rownames(data2)),],
                      by = c("ID","type","dataset"))

  data_merge %>%
    dplyr::filter(.data$type %in% sample_type) -> data_merge
  data_merge[4:(ncol(data_merge))] <- apply(data_merge[4:(ncol(data_merge))],2,as.numeric)
  result<-corr.test(data_merge[,4],data_merge[,5:(ncol(data_merge))],method = cor_method)
  result<-cbind(t(as.data.frame(result[1])),t(as.data.frame(result[4])))
  result<-data.frame("Symbol"= stringr::str_remove(row.names(result),"r."),"Correlation efficience"=result[,1],"P value"=result[,2])
  if (length(id2)==1 ){
    result[1,1] <- id2
  }
  row.names(result)<-NULL

  p_df<-list("cor_result" = result,"cor_data"=data_merge)
}

