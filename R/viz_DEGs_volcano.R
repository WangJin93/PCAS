#' @title Volcano plot for DEGs
#' @description
#' Plotting volcano plot for DEGs between tumor and normal samples in CPTAC datasets.
#' @import dplyr ggplot2
#' @param cohort Data cohort, for example, "LUAD_APOLLO", "LUAD_CPTAC".
#' @param data_input Expression data obtained from get_expr_data() function.
#' @examples
#' \dontrun{
#' results <- get_DEGs_result(dataset = "LUAD_CPTAC_protein", method = "limma")
#' results <- get_DEGs_result(dataset = "LUAD_CPTAC_mRNA", method = "t.test")
#' viz_DEGs_volcano(results)
#' }
#' @export
#'
viz_DEGs_volcano <- function(df,
                    p.cut=0.05,
                    logFC.cut=1,
                    show.top=FALSE,
                    show.labels=NULL){


  df$change<-ifelse(df$adj.P.Val<p.cut,ifelse(df$logFC < -logFC.cut,"Down",ifelse(df$logFC > logFC.cut,"Up","No")),"No")
  df<-arrange(df,logFC)
  if (show.top == TRUE){
    top5 <-c(df[c(1:5),1], df[c((nrow(df)-4):nrow(df)),1])
    df$label <- ifelse(df$Symbol %in% top5,df$Symbol,"")
  }
  if(!is.null(show.labels)){
    df$label <- ifelse(df$Symbol %in% show.labels,df$Symbol,"")
  }
  p <-  ggplot(data = df) +
    aes(x=logFC %>% as.numeric(),y=-log10(adj.P.Val %>% as.numeric())) +
    # geom_point(alpha = input$alphaInput, size = input$pointSize, shape = 16) +
    geom_point(size = 2, shape = 16,alpha=0.8) +


    # This needs to go here (before annotations)
    theme_light(base_size = 20)

  p <- p + aes(color=change) + scale_color_manual(values=c("blue","grey20","red"))

  p <- p + geom_hline(yintercept = c(-log10(p.cut)), linetype="dashed", color="grey30",size=0.8)+
    geom_vline(xintercept = c(logFC.cut,-logFC.cut), linetype="dashed", color="grey30",size=0.8)

  p <- p+ theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())+
    labs(x=bquote(''*Log[2]*' (Fold change)'), y=bquote(''*-Log[10]*' (P.adj value)'))

  ########## User defined labeling
  if (show.top|!is.null(show.labels)){
  p <-  p + geom_text_repel(
    data = df,
    aes(label = label),
    size = 5,
    color="black",
    nudge_x = 0.2,
    nudge_y=0.2,max.overlaps=100000,
    box.padding = unit(0.9, "lines"),
    point.padding = unit(.3+5*0.1, "lines"),show.legend=F
  )

  }


  p
}
