#' @title Visualizing genes expression
#' @description
#' Visualizing the different expression of mRNA/protein expression data between Tumor and Normal tissues in CPTAC database.
#' @import dplyr ggplot2 reshape2 ggpubr
#' @param df Gene expression data obtained from get_expr_data().
#' @param df_type The type of gene expression data, one value of "single","multi_gene", and "multi_set".
#' @param Show.P.value Whether to display the results of differential analysis, default TRUE.
#' @param Show.P.label Whether to display significance markers for differential analysis, default TRUE.
#' @param Method Methods of differential analysis, "t.test" or "limma", default "t.test".
#' @param values Color palette for normal and tumor groups. Default c("#00AFBB", "#FC4E07").
#' @param Show.n Display sample size.
#' @param Show.n.location Y-axis position displayed for sample size.
#' @examples
#' \dontrun{
#' df_single <- get_expr_data(datasets = "LUAD_CPTAC_mRNA",genes = c("TP53"))
#' df_multi_gene <- get_expr_data(datasets = "LUAD_CPTAC_protein",genes = c("TP53","TNS1"))
#' df_multi_set <- get_expr_data(datasets=c("CCRCC_CPTAC_protein","GBM_CPTAC_protein","HNSCC_CPTAC_protein","LSCC_CPTAC_protein","LUAD_CPTAC_protein" ,
#'                                          "PDAC_CPTAC_protein" , "UCEC_CPTAC2_protein", "UCEC_CPTAC1_protein"),genes = c("TNS1"))
#' viz_TvsN(df_single,df_type = "single")
#' viz_TvsN(df_multi_gene,df_type = "multi_gene")
#' viz_TvsN(df_multi_set,df_type = "multi_set")
#' }
#' @export
#'
viz_TvsN <- function(df,df_type = c("single","multi_gene","multi_set"),
                        Show.P.value = TRUE,
                        Show.P.label = TRUE,
                        Method = "t.test",
                        values = c("#00AFBB", "#FC4E07"),
                        Show.n = TRUE,
                        Show.n.location = "default") {
  if (is.null(df)) {
    return(NULL)
  }

  df <- df %>% dplyr::filter(type %in% c("Tumor","Normal"))

  ##计算每种肿瘤正常和肿瘤组织的样本量和差异分析
  if (df_type == "single"){
    colnames(df)[4] <- "value"
    df <- df  %>% na.omit()
    df$value <- as.numeric(df$value)
    if (Show.P.value == TRUE) {
      pv <- df %>%
        ggpubr::compare_means(value ~ type, data = ., method = Method, symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1), symbols = c( "***", "**", "*", "ns")),p.adjust.methods="BH")
      pv <- pv %>% dplyr::select(c("p", "p.signif", "p.adj"))
    }
    count_N<-df %>% group_by(type) %>% tally
    count_N$n <- paste("n =",count_N$n)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = type, y =value, fill = type)) +
      ggplot2::geom_violin(trim = FALSE,show.legend = F) +
      ggplot2::geom_boxplot(width = 0.2, fill = "white",show.legend = F) +
      ggplot2::xlab(NULL) +ggplot2::ylab("Expression") +
      ggplot2::scale_fill_manual(values = values)
    if (Show.n == TRUE) {
      if (Show.n.location == "default") Show.n.location <- min(df$value)-(max(df$value)-min(df$value))*0.2
      p <- p +ggplot2::geom_text(data=count_N, ggplot2::aes(label=n ,x = type,y = Show.n.location,color=values),
                                 position=ggplot2::position_dodge2(0.9),size = 6, hjust = 0.5,show.legend = F)+
        ggplot2::scale_colour_manual(values=values)
    }
    if (Show.P.value == TRUE & Show.P.label == TRUE) {
      p <- p + ggplot2::geom_text(ggplot2::aes(
        x = 1.5,
        y = max(df$value)+(max(df$value)-min(df$value))*0.1,
        label = .data$p.signif
      ),size=8,
      data = pv,
      inherit.aes = FALSE
      )
    }
    if (Show.P.value == TRUE & Show.P.label == FALSE) {
      p <- p + ggplot2::geom_text(ggplot2::aes(
        x = 1.5,
        y = max(df$value)+(max(df$value)-min(df$value))*0.1,
        label = ifelse(signif(.data$p, 3) < 0.001," P < 0.001",paste0(P =paste0("P = ",as.character(signif(.data$p, 3)))))
      ),size=6,
      data = pv,
      inherit.aes = FALSE
      )
    }
  }
  if (df_type == "multi_gene"){
    df <- reshape2::melt(df, measure.vars = colnames(df)[4:(ncol(df))])

    df$value <- as.numeric(df$value)
    df <- df  %>% na.omit()
    if (Show.P.value == TRUE) {
      pv <- df %>% as.data.frame() %>%
        ggpubr::compare_means(value ~ type, data = ., method = Method, group.by = "variable", symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1), symbols = c( "***", "**", "*", "ns")),p.adjust.methods="BH")
      pv <- pv %>% dplyr::select(c("variable", "p", "p.signif", "p.adj"))
      message("Counting P value finished")
    }
    count_N<-df %>% group_by(variable, type) %>% tally
    count_N$n <- paste("n =",count_N$n)
    p <- ggplot2::ggplot(df, aes(x = variable, y = value, fill = type)) +
      ggplot2::geom_boxplot() +
      ggplot2::xlab(NULL) +ggplot2::ylab("Expression") +
      ggplot2::scale_fill_manual(values = values)
    if (Show.n == TRUE) {
      if (Show.n.location == "default") Show.n.location <- min(df$value)-(max(df$value)-min(df$value))*0.2
      p <- p +geom_text(data=count_N, aes(label=n, y=Show.n.location,color=type), position=position_dodge2(0.9),size = 6,angle=90, hjust = 0,show.legend = F)+
        scale_colour_manual(values=values)
    }
    if (Show.P.value == TRUE & Show.P.label == TRUE) {
      p <- p + ggplot2::geom_text(aes(
        x = .data$variable,
        y = max(df$value)+(max(df$value)-min(df$value))*0.1,
        label = .data$p.signif
      ),size=6,
      data = pv,
      inherit.aes = FALSE
      )
    }
    if (Show.P.value == TRUE & Show.P.label == FALSE) {
      p <- p + ggplot2::geom_text(aes(
        x = .data$variable,
        y = max(df$value)+(max(df$value)-min(df$value))*0.1,
        label = as.character(signif(.data$p, 2))
      ),
      data = pv,
      inherit.aes = FALSE
      )
    }
  }
  if (df_type == "multi_set"){
    colnames(df)[4] <- "value"
    df$value <- as.numeric(df$value)
    df <- df  %>% na.omit()
    if (Show.P.value == TRUE) {
      pv <- df %>%
        ggpubr::compare_means(value ~ type, data = ., method = Method, group.by = "dataset", symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1), symbols = c( "***", "**", "*", "ns")),p.adjust.methods="BH")
      pv <- pv %>% dplyr::select(c("dataset", "p", "p.signif", "p.adj"))
    }
    count_N<-df %>% group_by(dataset, type) %>% tally
    count_N$n <- paste("n =",count_N$n)
    p <- ggplot(df, aes(x = dataset, y = value, fill = type)) +
      ggplot2::geom_boxplot() +
      ggplot2::xlab(NULL) +ggplot2::ylab("Expression") +
      ggplot2::scale_fill_manual(values = values)

    if (Show.n == TRUE) {
      if (Show.n.location == "default") Show.n.location <- min(df$value)-(max(df$value)-min(df$value))*0.2
      p <- p +geom_text(data=count_N, aes(label=n, y=Show.n.location,color=type), position=position_dodge2(0.9),size = 5,angle=90, hjust = 0,show.legend = F)+
        scale_colour_manual(values=values)
    }
    if (Show.P.value == TRUE & Show.P.label == TRUE) {
      p <- p + geom_text(aes(
        x = .data$dataset,
        y = max(na.omit(df$value)) * 1.1,
        label = .data$p.signif
      ),
      data = pv,
      inherit.aes = FALSE
      )
    }
    if (Show.P.value == TRUE & Show.P.label == FALSE) {
      p <- p + geom_text(aes(
        x = .data$dataset,
        y = max(na.omit(df$value)) * 1.1,
        label = as.character(signif(.data$p, 2))
      ),
      data = pv,
      inherit.aes = FALSE
      )
    }
  }

   p <- p + theme_bw() + theme(panel.grid=element_blank())+
     ggplot2::theme(
       axis.text.x = element_text(size = 16),
       axis.text.y = element_text(size = 16),
       axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16)
     )
   if (df_type == "multi_set") p <- p+theme(axis.text.x = element_text(angle = 45, hjust = 1.0),plot.margin =unit(c(0.2,0.2,0.2,2), 'cm'))
   return(p)
}
