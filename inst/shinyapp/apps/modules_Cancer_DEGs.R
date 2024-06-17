ui.modules_diff_gene <- function(id) {
  ns <- NS(id)
  fluidPage(
fluidRow(column(3,
                wellPanel(
                  shinyWidgets::prettyRadioButtons(ns("data_type"), "Select data type:",
                                                   choiceNames = c("Proteome","Transcriptome"),
                                                   choiceValues = c("Proteome","Transcriptome"),
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
                  conditionalPanel(ns=ns,
                                   condition = "input.data_type == 'Transcriptome'",
                                   shinyWidgets::prettyRadioButtons(ns("Biotype"), "Select Biotype:",
                                                                    choiceNames = c("All","Protein coding"),
                                                                    choiceValues = c("All","Protein coding"),
                                                                    selected = "Protein coding",
                                                                    shape = "round",inline = T,
                                                                    status = "success",
                                                                    animation = "tada"
                                   )
                  ),

        HTML("<hr>"),
        shinyWidgets::prettyRadioButtons(ns("method"), "Select Method:",
                                         choiceNames = c("limma","t-test"),
                                         choiceValues = c("limma","t.test"),
                                         selected = "limma",
                                         shape = "round",inline = T,
                                         status = "success",
                                         animation = "tada"
        ),


        numericInput(ns("logFC_cutoff"), "|logFC| threshold:", value = 1,min = 0),
        selectInput(
          ns("p_cutoff"),
          "Adjust P value threshold:",
          c("0.05","0.01","0.001","0.0001"),
          selected = "0.05",
          multiple = FALSE,
          selectize = TRUE,
          width = NULL,
          size = NULL
        ),
        selectizeInput(
          ns("show.lable"),
          "Symbols to show:",
          choices = NULL,
          multiple = TRUE,
          options = list(
            create = TRUE,
            placeholder = "Enter a Symbol symbol, e.g. TP53",
            plugins = list("restore_on_backspace")
          )
        ),

        selectInput(
          ns("vol.theme"),
          "Theme:",
          c("Default","classic","minimal"),
          selected = "Default",
          multiple = FALSE,
          selectize = TRUE,
          width = NULL,
          size = NULL
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
        )),
        wellPanel(
          numericInput(inputId = ns("height_scatter"), label = "Height", value = 600),
          numericInput(inputId = ns("width_scatter"), label = "Width", value = 600),
          downloadBttn(
            outputId = ns("download"),
            style = "gradient",
            color = "default",
            block = TRUE,
            size = "sm"
          )

        )
      ),
      column(9,
             shinycssloaders::withSpinner(plotOutput(ns("volcano"), height = "auto")),
             br(),



             HTML("<hr>"),


             hr(),
        # h5("NOTEs:"),
        # p("1. The data query may take some time based on your network. Wait until a plot shows"),
        # p("2. The DEGs were analysez by the merged datasets of Pancancer and GTEx database."),
        # p("3. The analysis method is limma"),
        tags$br(),
        DT::DTOutput(outputId = ns("tbl")),

      )
    )
  )
}

server.modules_diff_gene <- function(input, output, session) {
  ns <- session$ns
  observeEvent(input$data_type,{
    ddd <- unique(dataset_info %>% dplyr::filter(Data.type %in% input$data_type) %>% .["Primary.Site"])
    updateSelectInput(session,"cohorts_text",
                      label = "Cancer type:",
                      choices = ddd,
                      selected = "Lung"
    )

  })
  observe({

    updateSelectizeInput(
      session,
      "show.lable",
      choices = ID_list_pro %>% unlist() %>% unique(),
      selected = "TNS1",
      server = TRUE
    )
  })

  observeEvent({input$data_type
    input$cohorts_text},{
      ddd <- dataset_info %>% dplyr::filter(Data.type %in% input$data_type)%>%
        dplyr::filter(!is.na(Normal))%>%
        dplyr::filter(Primary.Site %in% input$cohorts_text)%>% .[,"Abbre"]
      updateSelectInput(session,"datasets_text",
                        label = "Select dataset:",
                        choices = ddd,
                        selected = ddd[1]
      )

    })




  plot_vis <- eventReactive(input$search_bttn, {
    #创建数据库连接
    expr <- get_DEGs_result(input$datasets_text,method = input$method)
    if (input$data_type =="Transcriptome" & input$Biotype == "Protein coding"){
      expr <- expr[which(expr$gene_type == "protein_coding"),]

    }
    results <- viz_DEGs_volcano(expr,p.cut=as.numeric(input$p_cutoff) ,logFC.cut=input$logFC_cutoff,show.labels=input$show.lable)

    return(results)

  })
  # Show waiter for plot

  w <- waiter::Waiter$new(id = "wb_plot0", html = waiter::spin_hexdots(), color = "black")
  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })

  output$volcano <- renderPlot(width = width_scatter,
                               height = height_scatter,{
    w$show() # Waiter add-ins

      p<-plot_vis()

      df <-  plot_vis()$data

      #remove legend (if selected)
      if (input$vol.theme=="classic"){
        p<-p+theme_classic()
      } else if(input$vol.theme=="minimal"){
        p<-p+theme_minimal()
      }
      p
  })
  output$download <- downloadHandler(
    filename = function() {
      paste0(input$datasets_text,"_",input$method,"_Volcano_plot.pdf")
    },
    content = function(file) {
      p <- plot_vis()
      #remove legend (if selected)
      if (input$vol.theme=="classic"){
        p<-p+theme_classic()
      } else if(input$vol.theme=="minimal"){
        p<-p+theme_minimal()
      }
        pdf(file,width = input$width_scatter/70 ,height = input$height_scatter/70)
        print(p)
        dev.off()
    })
  output$tbl <- DT::renderDataTable(server = T, {
    if (input$search_bttn){
      DT::datatable(
        plot_vis()$data,
        rownames = FALSE,
        extensions = c("Buttons"),
        options = list(
          pageLength = 5,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table", filename = paste0(input$datasets_text,"_",input$method),
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
