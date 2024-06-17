ui.modules_pancan_corr <- function(id) {
  ns <- NS(id)
  fluidPage(
        fluidRow(
          column(
            3,
            wellPanel(
              shinyWidgets::prettyRadioButtons(ns("data_type"), "Select data type:",
                                               choiceNames = c("Proteome","Transcriptome"),#"Phosphoproteome",
                                               choiceValues = c("Proteome","Transcriptome"),#"Phosphoproteome",
                                               selected = "Proteome",
                                               shape = "round",inline = T,
                                               status = "success",
                                               animation = "tada"
              ),
              selectInput(
                inputId = ns("datasets_text"),
                label = "Select dataset::",
                choices = "",
                multiple = T
              ),
              hr(),
              selectizeInput(
                inputId = ns("Pancan_search"),
                label = "Input a protein symbol",
                choices = "TNS1",
                width = "100%",
                options = list(
                  create = TRUE,
                  maxOptions = 5,
                  plugins = list("restore_on_backspace")
                )
              ),
              # shinyjs::hidden(
              #   selectizeInput(
              #     inputId = ns("phoso_site"), # molecule identifier
              #     label = "Select Phosphorylation Site:",
              #     choices=c("NP_000537.3:s315","NP_000537.3:s314","NP_000537.3:s392","NP_000537.3:y163s166t170","NP_000537.3:t150")
              #   )
              # ),
              shinyWidgets::prettyRadioButtons(ns("data_type2"), "Select data type of Variable B:",
                                               choiceNames = c("Proteome","Transcriptome"),#"Phosphoproteome",
                                               choiceValues = c("Proteome","Transcriptome"),#"Phosphoproteome",
                                               selected = "Proteome",
                                               shape = "round",inline = T,
                                               status = "success",
                                               animation = "tada"
              ),

              selectizeInput(
                inputId = ns("geneset"),
                label = "Gene sets:",
                selected = c("Immune checkpoint"),
                choices = c("Immune checkpoint","m6A methylation","Cell senescence","EMT","Cell senescence","Ferroptosis","Cuproptosis","User defined"),
                multiple = F
              ),
              textAreaInput(
                inputId = ns("genelist"), # molecule identifier
                label = "Input gene list:",
                value = "",
                placeholder = "Paste your genelist!\neg.\nTNS1\nTP53\nTXNIP",
                height = "200"
              )
    ),
      ),
      column(
        3,
        h3("Parameter settings"),
        wellPanel(
          selectInput(
            inputId = ns("sample_type"),
            label = "Only use normal or tumor",
            choices = c("Normal","Tumor"),
            selected = "Tumor", multiple = TRUE
          ),

          selectInput(
            inputId = ns("cor_method"),
            label = "Select Correlation method",
            choices = c("spearman", "pearson"),
            selected = "pearson"
          ),
          shinyWidgets::actionBttn(
            inputId = ns("search_bttn"),
            label = "Go!",
            style = "gradient",
            icon = icon("search"),
            color = "primary",
            block = TRUE,
            size = "sm"
          )
        ),
        hr(),
        h3("Download"),
        wellPanel(
          numericInput(inputId = ns("height_scatter"), label = "Height", value = 500,max = 600),
          numericInput(inputId = ns("width_scatter"), label = "Width", value = 800),
          downloadBttn(
            outputId = ns("download"),
            style = "gradient",
            color = "default",
            block = TRUE,
            size = "sm"
          )
        )
      ),
      column(6,
             shinycssloaders::withSpinner(plotOutput(ns("hm_gene_immune_cor"), height = "auto")),
        hr(),
        # h5("NOTEs:"),
        # p("1. ", tags$a("The immune checkpoint molecules (ICMs) were collected from the study conducted by Charoentong et al (43) who reported 24 immunoinhibitory genes and 45 immunostimulatory genes. ")),
        # p("2. ", tags$a(href = "https://pubmed.ncbi.nlm.nih.gov/28052254/", "Charoentong P, Finotello F, Angelova M, et al. Pan-cancer Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and Predictors of Response to Checkpoint Blockade. Cell Rep. 2017;18(1):248-262. doi:10.1016/j.celrep.2016.12.019")),
        DT::DTOutput(outputId = ns("tbl")),
      )
    )
  )
}

server.modules_pancan_corr <- function(input, output, session) {
  ns <- session$ns

  observeEvent({input$data_type},{
    ddd <- dataset_info %>% dplyr::filter(Data.type %in% input$data_type)%>%
      dplyr::filter(!is.na(Normal))%>% .[,"Abbre"]
    updateSelectInput(session,"datasets_text",
                      label = "Select dataset:",
                      choices = ddd,
                      selected = ddd[which(stringr::str_detect(ddd,"CPTAC"))]
    )

  })
  observe({

    updateSelectizeInput(
      session,
      "Pancan_search",
      choices = ID_list_pro %>% unlist() %>% unique(),
      selected = "TNS1",
      server = TRUE
    )
  })
  observe({
    if (input$data_type == "Phosphoproteome") {
      shinyjs::show(id = "phoso_site")
    }else{
      shinyjs::hide(id = "phoso_site")

    }
  })
  observe({
    genelist <- case_when(
      input$geneset == "Immune checkpoint" ~ c("SIRPA\nCTLA4\nTIGIT\nLAG3\nVSIR\nLILRB2\nSIGLEC7\nHAVCR2\nLILRB4\nPDCD1\nBTLA"),
      input$geneset == "m6A methylation" ~ c("METTL3\nMETTL14\nWTAP\nKIAA1429\nZC3H13\nRBM15\nRBM15B\nFTO\nALKBH5\nYTHDC1\nYTHDC2\nYTHDF1\nYTHDF2\nYTHDF3\nHNRNPC\nHNRNPA2B1\nIGF2BP1\nIGF2BP2\nIGF2BP3"),
      input$geneset == "Cell senescence" ~ c("TP53\nCDKN1A\nCDKN2A\nLMNB1\nMKI67"),
      input$geneset == "EMT" ~ c("CDH1\nCDH2\nSNAIL\nVIM"),
      input$geneset == "Ferroptosis" ~ c("ACSL4\nLPCAT3\nTFR1\nSLC7A11\nGPX4\nFTH1"),
      input$geneset == "Cuproptosis" ~ c("FDX1\nLIAS\nLIPT1\nDLD\nDLAT\nPDHA1\nPDHB\nMTF1\nGLS\nCDKN2A\nATP7B\nSLC31A1"),
      input$geneset == "User defined" ~ c("")
    )
    updateTextAreaInput(session, "genelist",
                        label = "Input gene list:",
                        value = genelist)
  })
  datasets2 <- reactive({
    ddd <- dataset_info %>% dplyr::filter(Data.type == input$data_type2)
    ddd[which(stringr::str_remove(ddd$Abbre,"_protein|_Phospho|_mRNA") %in% stringr::str_remove(input$datasets_text,"_protein|_mRNA")) ,]$Abbre
  })

  genelist <- reactive({strsplit(input$genelist,"\n")[[1]]})
  geneset_data <- eventReactive(input$search_bttn, {
    ga_ids <-  genelist()
    # print(ga_ids)
    dd <- get_expr_data(genes  = ga_ids,datasets = datasets2())
    if (is.null(dd)) {
      showModal(modalDialog(
        title = "Message",   easyClose = TRUE,
        "No results were returned, because the current dataset does not contain the gene symbols you entered."
      ))
      return(NULL)

    }else{
      if (length(colnames(dd))-3 < length(ga_ids)){
        showModal(modalDialog(
          title = "Message",   easyClose = TRUE,
          paste0("The selected datasets only contain the following gene symbols: ", paste(colnames(dd)[4:(ncol(dd))],collapse = " ,") )
        ))}
    }
    return(dd)
  })

  # Show waiter for plot
  w <- waiter::Waiter$new(id = ns("hm_gene_immune_cor"), html = waiter::spin_hexdots(), color = "white")
  plot_func <- eventReactive(input$search_bttn, {
    if (nchar(input$Pancan_search) >= 1) {
      Pancan_search <- input$Pancan_search
      df <- get_expr_data(genes = Pancan_search,datasets = input$datasets_text)
      if (is.null(df)){
        showModal(modalDialog(
          title = "Message",   easyClose = TRUE,
          "No results were returned, please comfirm your inputed gene symbol."
        ))
        return(NULL)
      }

      p <- cor_pancancer_genelist(df,geneset_data(), sample_type = input$sample_type,
                            cor_method = input$cor_method)
    }
    return(p)
  })

#
  plot_data <- eventReactive(input$search_bttn, {
      p <- plot_func()

    if (is.null(p)) return(NULL)

    rvalue_T<-as.data.frame(p$r) %>%  rownames_to_column("Gene")
    pvalue_T<-as.data.frame(p$p) %>%  rownames_to_column("Gene")
    pdata<-gather(as.data.frame(rvalue_T),Dataset,corr,-Gene)
    pvalue_T<-gather(as.data.frame(pvalue_T),Dataset,PValue,-Gene)
    pdata <- merge(pdata,pvalue_T,by=c("Dataset","Gene"))
    pdata$pstar <- ifelse(pdata$PValue < 0.05,
                          ifelse(pdata$PValue < 0.001, "***", ifelse(pdata$PValue < 0.01, "**", "*")),
                          ""
    )
    return(pdata)
  })
  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })

  output$hm_gene_immune_cor <- renderPlot(width = width_scatter,
                                          height = height_scatter,{
    w$show() # Waiter add-ins
                                            viz_cor_heatmap(plot_func()$r,plot_func()$p)
                                          })

  corr_func <- eventReactive(input$tbl_rows_selected, {
    s <- input$tbl_rows_selected
    if (length(s)) {
      selected <- plot_data()[s,]
      #data<- df[which(df$tissue == type),]
      checkpoint <- selected$Gene
      cancer <- selected$Dataset
      df <- plot_func()$sss[[cancer]]
      df <- df[,c("ID","type","dataset",colnames(df)[4],checkpoint)]
      return(df)
    }
    return(df)
  })


#
  observeEvent(eventExpr = input$tbl_rows_selected,{
    s <- input$tbl_rows_selected

    if (length(s)) {
      showModal(
        modalDialog(
          title = paste("Correlation scatter plot"),
          size = "l",
          fluidPage(
            column(
              12,
              plotOutput(ns("scatter_plot"),width = "100%",height = "auto",),
              DT::DTOutput(outputId = ns("scatter_corr")),

            )
          )
        )
      )

      output$scatter_plot <- renderPlot(width = 450,
                                        height = 300,{
                                          df<-corr_func() %>% na.omit()
                                          viz_corplot(df,colnames(df)[4],colnames(df)[5],method=input$cor_method,x_lab= " exppression",y_lab=" exppression")

                                        })
      output$scatter_corr <- DT::renderDataTable(server = TRUE,{
        DT::datatable(
          corr_func(),
          rownames = T,
          extensions = c("Buttons"),
          options = list(
            pageLength = 5,
            dom = "Bfrtip",
            buttons = list(
              list(
                extend = "csv", text = "Download table", filename =  paste0(names(corr_func())[4],"_",names(corr_func())[5],"_",corr_func()[1,3]),
                exportOptions = list(
                  modifier = list(page = "all")
                )
              )
            )
          )
        )

      }
      )
      ## downloadTable
    }
  })

#
#   observeEvent(input$search_bttn, {
#     if (nchar(input$Pancan_search) >= 1) {
#       shinyjs::show(id = "save_csv")
#     } else {
#       shinyjs::hide(id = "save_csv")
#     }
#   })
#
#
  output$tbl <- DT::renderDataTable(server = FALSE, {
    if (!is.null( plot_data())){
      DT::datatable(
        plot_data(),
        rownames = T,
        selection =  "single",
        extensions = c("Buttons"),
        options = list(
          pageLength = 5,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table", filename =  paste0(input$Pancan_search,  "_pancan_TIL"),
              exportOptions = list(
                modifier = list(page = "all")
              )
            )
          )
        )
      )
    }
  })
#
#
#   ## downloadTable
#   output$downloadTable <- downloadHandler(
#     filename = function() {
#       paste0(input$Pancan_search, "_", input$profile, "_pancan_ICG.csv")
#     },
#     content = function(file) {
#       write.csv(plot_data(), file, row.names = FALSE)
#     }
#   )
#
  output$download <- downloadHandler(
    filename = function() {
      paste0(input$Pancan_search, "_pancan_genelist_cor.pdf")
    },
    content = function(file) {
        pdf(file,  width =  input$width_scatter/70 ,height = input$height_scatter/70)
      print(viz_cor_heatmap(plot_func()$r,plot_func()$p))
      dev.off()
    }
  )
}
