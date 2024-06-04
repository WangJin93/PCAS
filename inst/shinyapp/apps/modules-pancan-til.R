ui.modules_pancan_til <- function(id) {
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
                                               shape = "round",
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
              shinyjs::hidden(
                selectizeInput(
                  inputId = ns("phoso_site"), # molecule identifier
                  label = "Select Phosphorylation Site:",
                  choices=c("NP_000537.3:s315","NP_000537.3:s314","NP_000537.3:s392","NP_000537.3:y163s166t170","NP_000537.3:t150")
                )
              ),
              selectizeInput(
                inputId = ns("immune_type"),
                label = "Immune infiltration algorithm:",
                selected = c("TIMER"),
                choices = TIL_map$algorithm %>% unique(),
                multiple = F
              ),
              selectInput(
          inputId = ns("cor_method"),
          label = "Select Correlation method",
          choices = c("spearman", "pearson"),
          selected = "pearson"
        ),
        tags$hr(style = "border:none; border-top:2px solid #5E81AC;"),
        shinyWidgets::actionBttn(
          inputId = ns("search_bttn"),
          label = "Go!",
          style = "gradient",
          icon = icon("search"),
          color = "primary",
          block = TRUE,
          size = "sm"
        ),
       ),
        wellPanel(
          numericInput(inputId = ns("height_scatter"), label = "Height", value = 500,max = 600),
          numericInput(inputId = ns("width_scatter"), label = "Width", value = 800),
        downloadBttn(
          outputId = ns("download"),
          style = "gradient",
          color = "default",
          block = TRUE,
          size = "sm"
        ))
      ),
      column(9,
             shinycssloaders::withSpinner(plotOutput(ns("hm_gene_immune_cor"), height = "auto")),
        hr(),
        # h5("NOTEs:"),
        # p("1. ", tags$a(href = "http://timer.cistrome.org/", "TIL data source")),
        # p("2. ", tags$a(href = "https://pancanatlas.xenahubs.net/", "Genomic profile data source")),
        DT::DTOutput(outputId = ns("tbl")),
      )
    )
  )
}

server.modules_pancan_til <- function(input, output, session) {
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

  # observe({
  #   updateSelectInput(session,"immune_sig",
  #                     label = "Cell types :",
  #                     choices = ddd,
  #                     selected = TIL_map[which(TIL_map$algorithm %in% input$immune_type),]$cell_type
  #   )
  #
  # })
  observeEvent({input$Pancan_search
    input$datasets_text},{
      ids <- idmap_protein[which(idmap_protein$Symbol == input$Pancan_search),]$row_names
      # ids <- ids[which(lapply(ids, count_id) > 1)]
      updateSelectInput(session, "phoso_site",
                        label = "Select Phosphorylation Site:",
                        choices =ids)
    }

  )

  # Show waiter for plot
  w <- waiter::Waiter$new(id = ns("hm_gene_immune_cor"), html = waiter::spin_hexdots(), color = "white")

  plot_func <- eventReactive(input$search_bttn, {
    if (nchar(input$Pancan_search) >= 1) {
      Pancan_search <- input$Pancan_search
      if (input$data_type == "Phosphoproteome") {
        Pancan_search <- input$phoso_site
      }
      df <- get_expr_data(genes = Pancan_search,datasets = input$datasets_text)
      if (is.null(df)){
        showModal(modalDialog(
          title = "Message",   easyClose = TRUE,
          "No results were returned, please comfirm your inputed gene symbol."
        ))
        return(NULL)
      }

      p <- cor_pancancer_TIL(df,
                              cor_method = input$cor_method,
                              TIL_type = input$immune_type
      )
    }
    return(p)
  })
  plot_data <- eventReactive(input$search_bttn, {
      p <- plot_func()
    if (is.null(p)) return(NULL)
    rvalue_T<-as.data.frame(p$r) %>%  rownames_to_column("immune_cells")
    pvalue_T<-as.data.frame(p$p) %>%  rownames_to_column("immune_cells")
    pdata<-gather(as.data.frame(rvalue_T),Dataset,corr,-immune_cells)
    pvalue_T<-gather(as.data.frame(pvalue_T),Dataset,PValue,-immune_cells)
    pdata <- merge(pdata,pvalue_T,by=c("Dataset","immune_cells"))
    pdata$pstar <- ifelse(pdata$PValue < 0.05,
                          ifelse(pdata$PValue < 0.001, "***", ifelse(pdata$PValue < 0.01, "**", "*")),
                          ""
    )
    return(pdata)
  })
  TIL_corr_func <- eventReactive(input$tbl_rows_selected, {
    s <- input$tbl_rows_selected
    if (length(s)) {
      selected <- plot_data()[s,]
      #data<- df[which(df$tissue == type),]
      checkpoint <- selected$immune_cells
      Dataset <- selected$Dataset
      df <- plot_func()$sss[[Dataset]]
      df <- df[,c("ID","type","dataset",colnames(df)[4],checkpoint)]
    }
    rownames(df)<-NULL
    return(df)
  })
  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })
  output$hm_gene_immune_cor <- renderPlot(width = width_scatter,
                                          height = height_scatter,{
    w$show() # Waiter add-ins
                                            if (!is.null(plot_func())){
                                              viz_cor_heatmap(plot_func()$r,plot_func()$p)

                                            }

  })

  observeEvent(eventExpr = input$tbl_rows_selected,{
    s <- input$tbl_rows_selected

    if (length(s)) {
      showModal(
        modalDialog(
          title = paste("Gene-Immunue infitration scatter plot", input$Pancan_search),
          size = "l",
          fluidPage(
            column(
              12,align = "center",
              plotOutput(ns("scatter_plot"),width = "100%",height = "auto",),
              DT::DTOutput(outputId = ns("scatter_corr")),
            )
          )
        )
      )

      output$scatter_plot <- renderPlot(width = 450,
                                        height = 300,{
                                          df<-TIL_corr_func() %>% na.omit()
                                          viz_corplot(df,colnames(df)[4],colnames(df)[5],method=input$cor_method,x_lab= " expression",y_lab= "")

                                        })
      # output$scatter_corr <- DT::renderDataTable(
      #   TIL_corr_func(), server = TRUE,selection = 'single'
      # )
      output$scatter_corr <- DT::renderDataTable(server = FALSE, {
        DT::datatable(
          TIL_corr_func(),
          rownames = T,
          extensions = c("Buttons"),
          options = list(
            pageLength = 5,
            dom = "Bfrtip",
            buttons = list(
              list(
                extend = "csv", text = "Download table", filename =  paste0(names(TIL_corr_func())[4],"_",names(TIL_corr_func())[5],"_",TIL_corr_func()[1,3]),
                exportOptions = list(
                  modifier = list(page = "all")
                )
              )
            )
          )
        )
      })

      ## downloadTable
    }
  })



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


  ## downloadTable
  # output$downloadTable <- downloadHandler(
  #   filename = function() {
  #     paste0(input$Pancan_search, "_", input$profile, "_pancan_TIL.csv")
  #   },
  #   content = function(file) {
  #     write.csv( plot_data(), file, row.names = FALSE)
  #   }
  # )

  output$download <- downloadHandler(
    filename = function() {
      paste0(input$Pancan_search, "_pancan_TIL_cor.pdf")
    },
    content = function(file) {
        pdf(file,  width =  input$width_scatter/70 ,height = input$height_scatter/70)
      viz_cor_heatmap(plot_func()$r,plot_func()$p)
      dev.off()
    }
  )
}
