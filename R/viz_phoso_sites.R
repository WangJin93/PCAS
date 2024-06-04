#' @title Phosphorylation site information
#' @description
#' Query phosphorylation site information of target proteins based on CPTAC database phosphorylation proteomics data or UniProt database.
#' @import drawProteins ggplot2
#' @param gene Gene/protein symbol.
#' @param phoso_infoDB Database for extracting phosphorylation site information. Default "CPTAC".
#' @examples
#' \dontrun{
#' viz_phoso_sites("TNS1")
#' viz_phoso_sites("YTHDC2",phoso_infoDB= "UniProt")
#' }
#' @export
#'
viz_phoso_sites <-  function(gene = "YTHDC2",phoso_infoDB="CPTAC"){
  if (!phoso_infoDB %in% c("UniProt","CPTAC")){
    message("The phoso_infoDB value only supports 'UniProt' and 'CPTAC'.")
    return(NULL)
  }
  uni_id <- uniport_map[which(uniport_map$Symbol == gene & uniport_map$Reviewed == "reviewed"),]$Entry %>% unique()
  drawProteins::get_features(uni_id)-> rel_json
  drawProteins::feature_to_dataframe(rel_json)-> rel_data
  draw_canvas(rel_data) -> p
  p <- draw_chains(p, rel_data,label_size = 5,labels = gene)
  if (phoso_infoDB == "CPTAC"){
    phoso_data <- idmap_protein[which(idmap_protein$Symbol == gene),]
    sites <- c()
    for (i in 1:nrow(phoso_data)) {
      id <- phoso_data$row_names[i]
      aa <- str_locate_all(id,"t|s|y")[[1]][,1]
      for (i in 1:length(aa)) {
        if (i== length(aa)){
          bb <- substr(id,aa[i],nchar(id))
        }else{
          bb<- substr(id,aa[i],aa[i+1]-1)
        }
      }
      sites <- c(sites,bb)

    }
    sites <- unique(sites)
    p_data <- data.frame(aa = substr(sites,1,1) %>% toupper(),phoso_site =paste0(uniport_map[which(uniport_map$Symbol == gene),]$From[1],":",  sites ),location = stringr::str_remove(sites,"s|y|t") %>% as.numeric(),order=1)
    p <- p +    geom_segment(data = p_data,aes(x = location,
                                               y = order+0.2,
                                               xend = location,
                                               yend = order + 0.3),
                             color = "grey50", linewidth=0.5, linetype = 1) +
      ggplot2::geom_point(data = p_data,
                          ggplot2::aes(x = location, y = order + 0.3), shape = 21,
                          colour = "black", fill = "red", size = 4, show.legend = F)
  }else{
    p_data <- phospho_site_info(rel_data)
    p <- p +    geom_segment(data = p_data,aes(x = begin,
                                               y = order+0.2,
                                               xend = begin,
                                               yend = order + 0.3),
                             color = "grey50", linewidth=0.5, linetype = 1) +
      ggplot2::geom_point(data = p_data,
                          ggplot2::aes(x = begin, y = order + 0.3), shape = 21,
                          colour = "black", fill = "red", size = 4, show.legend = F)
  }
  p <- draw_domains(p, rel_data,label_domains = F)
  p <- draw_regions(p, rel_data)
  p <- draw_motif(p, rel_data)
  # p <- draw_phospho(p, rel_data, size = 8)
  # white backgnd & change text size
  p <- p + theme_bw(base_size = 20) +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank()) +
    theme(axis.ticks =element_blank(),
          axis.text.y =element_blank()) +
    theme(panel.border =element_blank())


  p+theme(legend.position = "bottom")
  # plotly::ggplotly(p) %>%
  #   plotly::style(hovertext=paste0("site: ", p_data[,"phoso_site"]))
}
