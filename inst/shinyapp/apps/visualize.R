
themes_list <- list(
  "Base" = ggthemes::theme_few(),
  "cowplot" = cowplot::theme_cowplot(),
  "Light" = theme_light(),
  "Minimal" = theme_minimal(),
  "Classic" = theme_classic(),
  "Gray" = theme_gray(),
  "half_open" = cowplot::theme_half_open(),
  "minimal_grid" = cowplot::theme_minimal_grid()
)


p_survplot <- function(data, palette = "aaas", pval = TRUE, risk.table=T,conf.int = T,size=1,surv.median.line = c("none", "hv", "h", "v"),
                       legend.title="group",legend.labs=c("High","Low"),legend.pos= c('top', 'bottom', 'left', 'right', 'none'),
                       xlab = 'Time') {
  fit <- survival::survfit(survival::Surv(OS.time, OS.status) ~ group, data = data)
  p <- survminer::ggsurvplot(fit,
                             data = data, pval = pval, pval.method = TRUE,
                             palette = palette,  conf.int = conf.int,
                             size = 1, # change line size
                             font.legend = c(16, "black"),
                             font.x = c(16, "plain", "black"),
                             font.y = c(16, "plain", "black"),
                             font.tickslab = c(16, "plain", "black"),
                             xlab = xlab,risk.table.fontsize = 5,
                             risk.table = risk.table,pval.size = 6,
                             risk.table.y.text = FALSE, # show bars instead of names in text annotations
                             risk.table.col = "strata", # Risk table color by groups
                             ncensor.plot = F, # plot the number of censored subjects at time t
                             surv.plot.height = 0.7,
                             risk.table.height = 0.3,
                             legend.title = legend.title,
                             legend.labs = legend.labs,
                             legend = legend.pos,
                             surv.median.line = surv.median.line,
                             ggtheme = ggplot2::theme_classic()
  )
  p$plot  <- p$plot  + ggplot2::guides(color = ggplot2::guide_legend(ncol = 3))

  p$table <- p$table +  theme(axis.title =  element_text(size = 16,colour = "black"),axis.text = element_text(size = 16,colour = "black"))
  p
}

stage_plot<-function(clinic_data,feature="Tumor_Stage_simplify",is.cont = F,cont.cut){
  clinic <- clinic_data %>% dplyr::select(c("ID",feature,colnames(clinic_data)[4]))
  clinic[clinic == "Not Reported"] = NA
  clinic[clinic == "Unknown"] = NA
  clinic <- na.omit(clinic)
  if (is.cont){
    clinic[,feature] <- ifelse(clinic[,feature] > cont.cut, paste0("> ",cont.cut),paste0("≤ ",cont.cut))
  }

  clinic[,feature] <- factor(clinic[,feature])
  colnames(clinic)[2] <- "Group"
  group.n <- unique(clinic$Group)
  p<-ggplot(clinic,
            aes(x=Group,y=get(colnames(clinic_data)[4]),fill=Group))+
    geom_bar(stat="summary",fun=mean,position="dodge",color="black")+
    stat_summary(fun.data = 'mean_sd', geom = "errorbar", linewidth = 0.4,size=0.8,width=0.3,position = position_dodge(0.9),color="black")+
    labs(x=NULL,y='Relative protein expression')
  # if  ("None" %in% my.compair){
  #   p <- p
  # }else{
  #   if (group.n>2){
  #     all_compare=c()
  #     groups <- unique(clinic$Group)
  #     for (i in 1:(length(groups)-1)) {
  #       for (j in (i+1):length(groups)) {
  #         all_compare <- c(all_compare, paste(groups[i],"vs",groups[j]))
  #       }
  #     }
  #     my.compair <- strsplit(all_compare," vs ")
  #     p<-p+
  #       stat_compare_means(method = 't.test',label = "p.signif", comparisons = my.compair,size=5)
  #
  #   }else if (group.n==2){
  #     p<-p+ stat_compare_means(method = 't.test',label.y = NULL,label = "p.signif",comparisons = "all",size=5)
  #   }
  #
  # }

  df_p<-list("df"=clinic,"p"=p)
  return(df_p)
}



