#' @title Correlation analysis of immune cells infiltration
#' @description
#' Calculate the correlation between target gene expression and immune cells infiltration in multiple datasets.
#' @import dplyr psych
#' @param df The expression data of the target gene in multiple datasets, obtained by the get_expr_data() function.
#' @param cor_method Method for correlation analysis, default "pearson".
#' @param TIL_type Algorithm for calculating immune cell infiltration, default "TIMER".
#' @examples
#' \dontrun{
#' dataset <- c("CCRCC_CPTAC_protein","GBM_CPTAC_protein","HNSCC_CPTAC_protein","LSCC_CPTAC_protein","LUAD_CPTAC_protein","PDAC_CPTAC_protein","UCEC_CPTAC2_protein","UCEC_CPTAC1_protein")
#' df <- get_expr_data(genes = "TNS1",datasets = dataset)
#' result <- cor_pancancer_TIL(df,Target.pathway = c("Cell cycle"))
#' }
#' @export
#'
cor_pancancer_TIL <- function(df,
                             cor_method = "spearman",
                             TIL_type = c("TIMER")
) {

  if (is.null(df)){
    return(NULL)
  }
  df$dataset<- df$dataset %>% stringr::str_remove(.,"_protein|_Phospho|_mRNA")
  df <- df %>% dplyr::filter(.,type == "Tumor") %>%
    inner_join(.,TIL_CPTAC,by="ID")
  sss <- split(df, df$dataset)
  dataset <- names(sss)
  nrow<-length(dataset)
  sig <- TIL_map[which(TIL_map$algorithm %in% TIL_type),]$cell_type
  ncol<-length(sig) #自己定矩阵的大小
  rvalue<-matrix(nrow=nrow,ncol=ncol)
  rownames(rvalue)<-dataset
  colnames(rvalue)<-sig
  pvalue<-matrix(nrow=nrow,ncol=ncol)
  rownames(pvalue)<-dataset
  colnames(pvalue)<-sig
  for(i in 1:length(dataset)){
    sss_can <- sss[[i]]
    if (nrow(sss_can)<4) {
      next
    }else{
      rvalue[i,] = corr.test(x =sss_can[,colnames(df)[4]], y = sss_can[,c(sig)],method = cor_method)[["r"]]
      pvalue[i,] = corr.test(x =sss_can[,colnames(df)[4]], y = sss_can[,sig],method = cor_method)[["p"]]

    }
  }

  rvalue_T<-t(na.omit(rvalue))
  pvalue_T<-t(na.omit(pvalue))
  plist = list()
  plist[[1]] = rvalue_T
  plist[[2]] = pvalue_T
  plist[[3]] = sss
  names(plist) = c("r","p","sss")
  # print(plist)
  return(plist)
}
