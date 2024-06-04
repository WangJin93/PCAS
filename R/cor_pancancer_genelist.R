#' @title Perform correlation analysis
#' @description
#' Perform correlation analysis of the mRNA/protein expression data in CPTAC database.
#' @import dplyr psych
#' @param df The expression data of the target gene in multiple datasets, obtained by the get_expr_data() function.
#' @param geneset_data The expression data of a genelist in multiple datasets, obtained by the get_expr_data() function.
#' @param sample_type Sample type used for correlation analysis, default all types: c("Tumor", "Normal").
#' @param cor_method Method for correlation analysis, default "pearson".
#' @examples
#' \dontrun{
#' genelist <- c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7","HAVCR2","LILRB4","PDCD1","BTLA")
#' dataset <- c("CCRCC_CPTAC_protein","GBM_CPTAC_protein","HNSCC_CPTAC_protein","LSCC_CPTAC_protein","LUAD_CPTAC_protein","PDAC_CPTAC_protein","UCEC_CPTAC2_protein","UCEC_CPTAC1_protein")
#' df <- get_expr_data(genes = "TNS1",datasets = dataset)
#' geneset_data <- get_expr_data(genes = genelist ,datasets = dataset)
#' result <- cor_pancancer_genelist(df,geneset_data,sample_type = c("Tumor"))
#' }
#' @export
#'
cor_pancancer_genelist <- function(df,
                                  geneset_data,
                                  sample_type = c("Tumor","Normal"),
                                  cor_method = "pearson") {

  if (is.null(df)){
    return(NULL)
  }
  df$dataset<- df$dataset %>% stringr::str_remove(.,"_protein|_Phospho|_mRNA")
  geneset_data$dataset <- geneset_data$dataset %>% stringr::str_remove(.,"_protein|_Phospho|_mRNA")
  df <- merge(df, geneset_data,  by = c("ID","type","dataset"))

  df <- df %>% dplyr::filter(.,type %in% sample_type)

  sss <- split(df, df$dataset)
  dataset <- names(sss)
  nrow<-length(dataset)
  sig <- colnames(df)[5:(ncol(df))]
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
      rvalue[i,] = psych::corr.test(x =sss_can[,colnames(df)[4]], y = sss_can[,sig],method = cor_method)[["r"]]
      pvalue[i,] = psych::corr.test(x =sss_can[,colnames(df)[4]], y = sss_can[,sig],method = cor_method)[["p"]]

    }
  }

  rvalue_T<-t(rvalue)
  pvalue_T<-t(pvalue)
  plist = list()
  plist[[1]] = rvalue_T
  plist[[2]] = pvalue_T
  plist[[3]] = sss
  names(plist) = c("r","p","sss")
  return(plist)
}
