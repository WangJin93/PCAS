ui.modules_Cancer_corr <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "",
    icon = icon("database"),
    useShinyjs(),  # Set up shinyjs

    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("Variable A:"),
        shinyWidgets::prettyRadioButtons(ns("data_type"), "Select data type:",
                                          choiceNames = c("Proteome","Phosphoproteome","Transcriptome"),
                                          choiceValues = c("Proteome","Phosphoproteome","Transcriptome"),
                                          selected = "Proteome",
                                          shape = "round",
                                          status = "success",
                                          animation = "tada"
        ),
        selectInput(
          inputId = ns("cohorts_text"),
          label = "Cancer type:",
          choices = "",
          multiple = T
        ),
        selectInput(
          inputId = ns("datasets_text"),
          label = "Select dataset::",
          choices = "",
          multiple = F
        ),
        selectizeInput(
          inputId = ns("ga_id"), # molecule identifier
          label = "Input gene A symbol:",
          choices = c("TP53","TNS1"),
          options = list(
            create = TRUE,
            maxOptions = 5,
            placeholder = "e.g. TP53",
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
        hr(),
        h4("Variable B:"),
        selectInput(
          inputId = ns("datasets_text2"),
          label = "Select dataset:",
          choices = "",
          multiple = FALSE
        ),
        selectizeInput(
          inputId = ns("ga_cor_id"), # molecule identifier
          label = "Input gene B symbol:",
          choices = c("TP53","TNS1"),
          multiple = TRUE,
          options = list(
            create = TRUE,
            placeholder = "You can enter multiple genes",
            plugins = list("restore_on_backspace")
          )
        ),
        shinyjs::hidden(
          selectizeInput(
          inputId = ns("phoso_site2"), # molecule identifier
          label = "Select Phosphorylation Site:",
          multiple = TRUE,

          choices=c("")
        )
        )

      ),
      mainPanel(
        fluidPage(
          fluidRow(
            column(
              4,
              wellPanel(
                h4("Parameters:"),

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
                actionBttn(
                  inputId = ns("ga_cor_submit"),
                  label = "Submit",
                  style = "gradient",
                  icon = icon("check"),
                  color = "default",
                  block = TRUE,
                  size = "sm"
                )
              ),
              tags$br(),
              wellPanel(
                h4("Download figure:"),
                numericInput(inputId = ns("cor_height"), label = "Height", value = 5),
                numericInput(inputId = ns("cor_width"), label = "Width", value = 8),
                downloadBttn(
                  outputId = ns("cor_download"),
                  style = "gradient",
                  color = "default",
                  block = TRUE,
                  size = "sm"
                )
              )
            ),
            column(
              8,

              shinycssloaders::withSpinner(DT::dataTableOutput(ns("ga_cor_data"))),
              h4("Scatter plot of selected row:"),
              column(
                11,
                shinycssloaders::withSpinner(plotOutput(ns("ga_cor_output"),width = "600",height = "400")),
              ),
              column(
                11,
                shinycssloaders::withSpinner(DT::dataTableOutput(ns("cor_scatter_data")))
              )
            ),


          ))
      )
    )
  )
}

server.modules_Cancer_corr <- function(input, output, session) {
  ns <- session$ns
  # df <- dataset


  observeEvent(input$data_type,{
    ddd <- unique(dataset_info %>% dplyr::filter(Data.type %in% input$data_type) %>% .["Primary.Site"])
    updateSelectInput(session,"cohorts_text",
                        label = "Cancer type:",
                        choices = ddd,
                      selected = "Lung"
                      )

  })
  observeEvent({input$data_type
    input$cohorts_text},{
    ddd <- dataset_info %>% dplyr::filter(Data.type %in% input$data_type)%>%
      dplyr::filter(Primary.Site %in% input$cohorts_text)%>% .[,"Abbre"]
      updateSelectInput(session,"datasets_text",
                      label = "Select dataset:",
                      choices = ddd,
                      selected = ddd[1]
    )

  })
  observe({

    updateSelectizeInput(
      session,
      "ga_cor_id",
      choices = ID_list_pro %>% unlist() %>% unique(),
      selected = "TNS1",
      server = TRUE
    )
    updateSelectizeInput(
      session,
      "ga_id",
      choices = ID_list_pro %>% unlist() %>% unique(),
      selected = "TP53",
      server = TRUE
    )
  })

  observe({
    if (!is.null(input$datasets_text)){
          if (stringr::str_detect(input$datasets_text,"Phospho")) {
      shinyjs::show(id = "phoso_site")
    }else{
      shinyjs::hide(id = "phoso_site")

    }
    }
    if (!is.null(input$datasets_text2)){
      if (stringr::str_detect(input$datasets_text2,"Phospho")) {
        shinyjs::show(id = "phoso_site2")
      }else{
        shinyjs::hide(id = "phoso_site2")

      }
    }



  })


  observeEvent({input$ga_id
    input$datasets_text},
               updateSelectInput(session, "phoso_site",
                                 label = "Select Phosphorylation Site:",
                                 choices =intersect(
                                   phoso_id[[input$datasets_text]],
                                   idmap_protein[which(idmap_protein$Symbol == input$ga_id),]$row_names)
               )
  )
  observeEvent({input$ga_cor_id
    input$datasets_text2}, {
    ddd <- intersect( phoso_id[[input$datasets_text2]],
                                   idmap_protein[which(idmap_protein$Symbol %in% input$ga_cor_id),]$row_names)
               updateSelectInput(session, "phoso_site2",
                                 label = "Select Phosphorylation Site:",
                                 choices =  ddd,
                                 selected = ddd[1]
               )}
  )





  observeEvent(input$datasets_text,{
    if (input$datasets_text != ""){
      updateSelectInput(session, "datasets_text2",
                        label = "Select dataset:",
                        choices  = dataset_info[which(stringr::str_detect(dataset_info$Abbre, stringr::str_remove(input$datasets_text,"_protein|_Phospho|_mRNA"))) ,]$Abbre
      )


    }
  })




  #########correlation
  # Show waiter for plot
  w1 <- waiter::Waiter$new(id = ns("ga_cor_output"), html = waiter::spin_hexdots(), color = "white")
  plot_cor_func <- eventReactive(input$ga_cor_submit, {
    dataset1 <- input$datasets_text
    id1 <- input$ga_id
    if (stringr::str_detect(dataset1,"Phospho")) {
      id1 <- input$phoso_site
    }
    dataset2 <- input$datasets_text2
    id2 <- input$ga_cor_id
    if (stringr::str_detect(dataset2,"Phospho")) {
      id2 <- input$phoso_site2
    }

    results <- cor_cancer_genelist(dataset1,
                                    id1,
                                    dataset2,
                                    id2,
                                    input$sample_type,
                                    input$cor_method)
    return(results)
  })

corr_func <- eventReactive(input$ga_cor_data_rows_selected , {
    s <- input$ga_cor_data_rows_selected
    if (length(s)) {
      df<-plot_cor_func()$cor_data[c(1,4,s+4)]
    }
    # datasets_text1 <- input$datasets_text
    # id1 <- input$ga_id
    # if (stringr::str_detect(datasets_text1,"Phospho")) {
    #   id1 <- input$phoso_site
    # }
    # if (stringr::str_detect(datasets_text1,"mRNA")) {
    #   id1 <- idmap_RNA[which(idmap_RNA$Symbol == id1),]$row_names
    # }
    # # data_select1 <- get_expr_data(datasets_text1, id1)
    # # colnames(data_select1)[1] <- id1
    # # data_select1 <- tibble::rownames_to_column(data_select1,var = "ID")
    # #
    # datasets_text2 <- input$datasets_text2
    # id2 <- input$ga_cor_id
    # if (stringr::str_detect(datasets_text2,"Phospho")) {
    #   id2 <- input$phoso_site2
    # }
    # if (stringr::str_detect(datasets_text2,"mRNA")) {
    #   id2 <- idmap_RNA[which(idmap_RNA$Symbol %in% id2),]$row_names
    # }
    # colnames(df)[2:3] <- c(id1,id2)
    return(df)


  })




output$ga_cor_output <- renderPlot(width = 500,
                                   height = 300,{
                                      w1$show() # Waiter add-ins

                                     df <- corr_func()%>% na.omit()
                                     viz_corplot(df,colnames(df)[2],colnames(df)[3],method=input$cor_method,x_lab = " expression",y_lab = " expression")
                                    }
)
output$cor_download <- downloadHandler(
  filename = function() {
    paste0(colnames(corr_func())[2],"_",colnames(corr_func())[3],".pdf")
  },
  content = function(file) {
    df <- corr_func()%>% na.omit()

      pdf(file = file, onefile = FALSE,height = input$cor_height,width = input$cor_width)
    print(viz_corplot(df,colnames(df)[2],colnames(df)[3],method=input$cor_method,
                      x_lab = " expression",y_lab = " expression"))
      dev.off()

    # ggplot2::ggsave(
    #   filename = file, plot = print(p_surv(), newpage = F), device = input$device,
    #   units = "cm", width = input$width, height = input$height, dpi = 600
    # )
  }
)
  # observeEvent(eventExpr = input$ga_cor_data_rows_selected,{
  #   message(input$Pancan_search, " is queried by user from home search box.")
  #   s <- input$ga_cor_data_rows_selected
  #
  #   if (length(s)) {
  #     showModal(
  #       modalDialog(
  #         title = paste("Correlation scatter plot"),
  #         size = "l",
  #         fluidPage(
  #           div(style="display:flex",
  #               numericInput(inputId = ns("heigh_scatter"), label = "Height", value = 400,max = 600,width = 100),
  #               numericInput(inputId = ns("widt_scatter"), label = "Width", value = 500,width = 100),
  #               textInput(inputId = ns('x_lab'), "X-axis suffix", value = ' expression',width = 250),
  #               textInput(inputId = ns('y_lab'), "Y-axis suffix", value = ' expression',width = 250)
  #
  #           ),
  #           column(
  #             12,
  #             plotOutput(ns("scatter_plot"),width = "100%",height = "auto",),
  #             DT::DTOutput(outputId = ns("scatter_corr")),
  #             wellPanel(
  #               id = ns("save_scatter_csv"),
  #               downloadButton(ns("downloadTable2"), "Save as csv")
  #             )
  #
  #           ),
  #           column(
  #             12,
  #             h4("NOTEs:"),
  #             h5("1. The data query may take some time based on your network. Wait until a plot shows"),
  #             h5("2. The unit of gene expression is log2(tpm+0.001)"),
  #           )
  #         )
  #       )
  #     )
  #     width_scatter <- reactive ({ input$widt_scatter })
  #     height_scatter <- reactive ({ input$heigh_scatter })
  #     output$scatter_plot <- renderPlot(width = width_scatter,
  #                                       height = height_scatter,{
  #                                         df<-plot_cor_func()$data_select[c(1,2,s+3)]
  #                                         viz_corplot(df,colnames(df)[2],colnames(df)[3],method=input$cor_method,x_lab=input$x_lab,y_lab=input$y_lab)
  #
  #                                       })
  #     output$scatter_corr <- DT::renderDataTable(
  #       plot_cor_func()$data_select[c(1,2,s+3)], server = TRUE,selection = 'single'
  #     )
  #     ## downloadTable
  #     output$downloadTable2 <- downloadHandler(
  #       filename = function() {
  #         paste0(input$ga_id, "_", input$ga_cor_id, "_cor.csv")
  #       },
  #       content = function(file) {
  #         write.csv(corr_func(), file, row.names = FALSE)
  #       }
  #     )
  #   }
  # })
  # output$ga_cor_output <- renderPlot(width = 500,
  #                                    height = 300,{
  #                                      w1$show() # Waiter add-ins
  #                                      corr_func()$p
  #                                    }
  # )
  # output$cor_download <- downloadHandler(
  #   filename = function() {
  #     paste0(input$ga_id, " correlation in ",substr(input$datasets_text,6,9),".", input$cor_device)
  #   },
  #   content = function(file) {
  #     if(input$cor_device == "png"){
  #       png(filename = file, units = "cm", width = input$cor_width, height = input$cor_height, res = 600)
  #       print(corr_func()$p)
  #       dev.off()
  #     }
  #     if(input$cor_device == "pdf"){
  #       pdf(file = file, onefile = FALSE)
  #       print(corr_func()$p)
  #       dev.off()
  #     }
  #     # ggplot2::ggsave(
  #     #   filename = file, plot = print(p_surv(), newpage = F), device = input$device,
  #     #   units = "cm", width = input$width, height = input$height, dpi = 600
  #     # )
  #   }
  # )
  output$ga_cor_data <- DT::renderDataTable(server = FALSE, {
    if (input$ga_cor_submit>0){
      DT::datatable(
        plot_cor_func()$cor_result,
        rownames = FALSE,selection = 'single',
        extensions = c("Buttons"),
        options = list(
          pageLength = 5,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table", filename = paste(input$ga_id, "correlation results"),
              exportOptions = list(
                modifier = list(page = "all")
              )
            )
          )
        )
      )

    }
  })
  output$cor_scatter_data <- DT::renderDataTable(server = FALSE, {
    s <- input$ga_cor_data_rows_selected
    if (length(s)) {
      DT::datatable(
        corr_func(),
        rownames = T,
        extensions = c("Buttons"),
        options = list(
          pageLength = 5,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table", filename = paste0(colnames(corr_func())[1],"_",colnames(corr_func())[2]),
              exportOptions = list(
                modifier = list(page = "all")
              )
            )
          )
        )
      )

      }
  })


  }
