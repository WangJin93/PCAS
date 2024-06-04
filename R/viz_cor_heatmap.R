#' @title Visualization of correlation results
#' @description
#' Presenting correlation analysis results using heat maps based on ggplot2.
#' @import dplyr psych ggplot2
#' @param r The correlation coefficient matrix r of the correlation analysis results obtained from the functions cor_pancancer_genelist(), cor_pancancer_TIL(), and cor_pancancer_drug().
#' @param p The P-value matrix p of the correlation analysis results obtained from the functions cor_pancancer_genelist(), cor_pancancer_TIL(), and cor_pancancer_drug().
#' @examples
#' \dontrun{
#' genelist <- c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7","HAVCR2","LILRB4","PDCD1","BTLA")
#' dataset <- c("CCRCC_CPTAC_mRNA","GBM_CPTAC_mRNA","HNSCC_CPTAC_mRNA","LSCC_CPTAC_mRNA","LUAD_CPTAC_mRNA","PDAC_CPTAC_mRNA","UCEC_CPTAC2_mRNA","UCEC_CPTAC1_mRNA")
#' df <- get_expr_data(genes = "TNS1",datasets = dataset)
#' geneset_data <- get_expr_data(genes = genelist ,datasets = dataset)
#' result <- cor_pancancer_genelist(df,geneset_data,sample_type = c("Tumor"))
#' viz_cor_heatmap(result$r,result$p)
#' }
#' @export
#'
viz_cor_heatmap <- function(r,p){
  r[is.na(r)] <- 0
  if (nrow(r)>1){
    gg <- hclust(dist(r))    #对行聚类
    r <- r[gg$order,]     #行，按照聚类结果排序ggplot(melt_data, aes(tumor, transcript))
  }

  melt_corr <- reshape2::melt(as.matrix(r))
  melt_p <- reshape2::melt(as.matrix(p))
  melt_data <- merge(melt_corr,melt_p,by=c("Var1","Var2"),sort = F)
  colnames(melt_data)=c("type","cell_type","corr","PValue")
  melt_data$text <- case_when(
    melt_data$PValue < 0.001 ~  "***",
    melt_data$PValue < 0.01 ~  "**",
    melt_data$PValue < 0.05 ~  "*",
    melt_data$PValue > 0.05 ~  "")
  melt_data$cell_type <- factor(melt_data$cell_type, levels = unique(melt_data$cell_type)) # 重新排序因子，决定坐标轴出图顺序

  melt_data$cell_type <- stringr::str_remove(melt_data$cell_type, "_XCELL")
  #######plot
  heat<-ggplot(melt_data, aes( cell_type,type)) +
    geom_tile(aes(fill = corr), colour = "grey", linewidth = 1)+
    scale_fill_gradient2(low = "#5C5DAF",mid = "white",high = "#EA2E2D") + # 这里可以用 windowns 小工具 takecolor 取色，看中哪个文章就吸哪个文章
    # 比如这篇 https://www.nature.com/articles/nmeth.1902
    theme_minimal() + # 不要背景
    theme(axis.title.x=element_blank(), # 去掉 title
          axis.ticks.x=element_blank(), # 去掉x 轴
          axis.title.y=element_blank(), # 去掉 y 轴
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.0,size = 14), # 调整x轴文字，字体加粗
          axis.text.y = element_text( size = 14),
          plot.margin = margin(l = 30,  unit = "pt")) + #调整y轴文字
    geom_text(aes(label=text),col ="black",size = 7) +
    labs(fill =paste0("\n","\n","\n","\n","\n","  *   p < 0.05","\n"," **  p < 0.01","\n",
                      "*** p < 0.001","\n\n","Correlation")) +   # 修改 legend 内容
    scale_x_discrete(position = "bottom")+ # 将 X 轴放置在最上面
    scale_y_discrete(position = c("right"))
  if (nrow(r)>1){
    h <- ggtree(gg,layout = "rectangular",branch.length = "none") # 绘制列聚类树
    heat<-heat %>% insert_left(h,width = 0.1)
  }
  heat
}
