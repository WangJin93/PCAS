#' @title Correlation analysis of drug sensitivity
#' @description
#' Calculate the correlation between target gene expression and anti-tumor drug sensitivity in multiple datasets.
#' @import dplyr psych
#' @param df The expression data of the target gene in multiple datasets, obtained by the get_expr_data() function.
#' @param cor_method Method for correlation analysis, default "pearson".
#' @param Target.pathway The signaling pathways of anti-tumor drug targets, default "Cell cycle". Use "drug_info"to get the detail infomation of these drugs.
#' @examples
#' \dontrun{
#' dataset <- c("CCRCC_CPTAC_protein","GBM_CPTAC_protein","HNSCC_CPTAC_protein","LSCC_CPTAC_protein","LUAD_CPTAC_protein","PDAC_CPTAC_protein","UCEC_CPTAC2_protein","UCEC_CPTAC1_protein")
#' df <- get_expr_data(genes = "TNS1",datasets = dataset)
#' result <- cor_pancancer_drug(df,Target.pathway = c("Cell cycle"))
#' }
#' @export
#'
cor_pancancer_drug <- function(df,
                              cor_method = "pearson",
                              Target.pathway = c("Cell cycle")) {

  if (is.null(df)){
    return(NULL)
  }
  df$dataset<- df$dataset %>% stringr::str_remove(.,"_protein|_Phospho|_mRNA")
  df <- merge(df, drug_CPTAC %>% tibble::rownames_to_column("ID"),  by = c("ID"))

  sss <- split(df, df$dataset)
  dataset <- names(sss)
  nrow<-length(dataset)
  sig <- drug_info[which(drug_info$Target.pathway %in% Target.pathway),][c("ID","Name")]
  sig$drug <- paste(sig$Name,sig$ID,sep="_")
  sig <- sig$drug
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
      rvalue[i,] = psych::corr.test(x =sss_can[,colnames(df)[4]], y = sss_can[,c(sig)],method = cor_method)[["r"]]
      pvalue[i,] = psych::corr.test(x =sss_can[,colnames(df)[4]], y = sss_can[,sig],method = cor_method)[["p"]]

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
