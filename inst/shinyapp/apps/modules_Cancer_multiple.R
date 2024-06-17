ui.modules_multi_gene <- function(id) {
  ns <- NS(id)
  fluidPage(
fluidRow(column(3,
                wellPanel(
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
                  hr(),
                  textAreaInput(
                    inputId = ns("ga_ids"), # molecule identifier
                    label = "Input gene symbols:",
                    value = "YTHDC1\nYTHDC2\nYTHDF1\nYTHDF2\nYTHDF3\nIGF2BP1\nIGF2BP2\nIGF2BP3",
                    placeholder = "Paste your genelist!\neg.\nTNS1\nTP53\nTXNIP",
                    height = "200"
                  ) ,
                  selectInput(inputId = ns("method"), label = "Select statistic method", choices = c("wilcox.test", "t.test"), selected = "t.test"),

        materialSwitch(ns("Show.P.value"), "Show P value", inline = TRUE),
        materialSwitch(ns("Show.P.label"), "Show P label", inline = TRUE),
        materialSwitch(ns("Show.n"), "Show sample size", inline = T),
        materialSwitch(ns("x_axis45"), "X-axis 45°", inline = T),
        colourpicker::colourInput(inputId = ns("tumor_col"), "Tumor sample color", "#00AFBB"),
        colourpicker::colourInput(inputId = ns("normal_col"), "Normal sample color", "#FC4E07"),
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
             shinycssloaders::withSpinner(plotOutput(ns("gene_pancan_dist"), height = "auto")),
             hr(),
             # h5("NOTEs:"),
             # p("1. The data query may take some time based on your network. Wait until a plot shows"),
             # p("2. You have to turn on both 'Show P value' and 'Show P label' to show significant labels"),
             tags$br(),
             DT::DTOutput(outputId = ns("tbl")),
        #      shinyjs::hidden(
        #        wellPanel(
        #          id = ns("save_csv"),
        #          downloadButton(ns("downloadTable"), "Save as csv")
        #        )
        # )
      )
    )
  )
}

server.modules_multi_gene <- function(input, output, session) {
  ns <- session$ns

  observeEvent(input$data_type,{
    ddd <- unique(dataset_info %>% dplyr::filter(Data.type %in% input$data_type) %>% .["Primary.Site"])
    updateSelectInput(session,"cohorts_text",
                      label = "Cancer type:",
                      choices = ddd,
                      selected = "Lung"
    )
    if (input$data_type=="Phosphoproteome"){
      updateTextAreaInput(session,"ga_ids",
                          label = "Input Phosphorylation sites:",
                          value = "NP_073739.3:s1090\nNP_073739.3:s1201\nNP_073739.3:s1202\nNP_073739.3:s1208\nNP_073739.3:s1221\nNP_073739.3:s1267\nNP_073739.3:s1269\nNP_073739.3:s1277\nNP_073739.3:s1279\nNP_073739.3:s1281",
                          placeholder = "Paste your genelist!\neg.\nNP_073739.3:s1221\nNP_073739.3:s1267"
      )

    }else{
      updateTextAreaInput(session,"ga_ids",
                          label = "Input gene symbols:",
                          value = "YTHDC1\nYTHDC2\nYTHDF1\nYTHDF2\nYTHDF3\nIGF2BP1\nIGF2BP2\nIGF2BP3",
                          placeholder = "Paste your genelist!\neg.\nTNS1\nTP53\nTXNIP"
      )

    }

  })
  observeEvent({input$data_type
    input$cohorts_text},{
      ddd <- dataset_info %>% dplyr::filter(Data.type %in% input$data_type)%>%
        dplyr::filter(Primary.Site %in% input$cohorts_text)%>%
        dplyr::filter(!is.na(Normal)) %>% .[,"Abbre"]
      updateSelectInput(session,"datasets_text",
                        label = "Select dataset:",
                        choices = ddd,
                        selected = ddd[1]
      )

    })


  colors <- reactive({
    c(input$tumor_col, input$normal_col)
  })


  ga_ids <- reactive({strsplit(input$ga_ids,"\n")[[1]]})


  # Show waiter for plot
  w <- waiter::Waiter$new(id = ns("gene_pancan_dist"), html = waiter::spin_hexdots(), color = "black")

  plot_func <- eventReactive(input$search_bttn, {
    if (length(ga_ids()) >= 1) {
      df <- get_expr_data(datasets = input$datasets_text,genes = ga_ids())
      p <- viz_TvsN(df,df_type = "multi_gene",
                     Method =  input$method,
                     Show.P.value = input$Show.P.value,
                     Show.P.label = input$Show.P.label,
                     Show.n = input$Show.n,
                     values = colors())
      if (input$x_axis45){
        p <- p+ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      }
    }
    return(p)
  })

  output$colorvalues <- reactive({
    c(input$tumor_col, input$normal_col)
  })

  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })
  output$gene_pancan_dist <- renderPlot(width = width_scatter,
                                        height = height_scatter,{
    w$show() # Waiter add-ins
                                          plot_func()
                                        })

  output$downloadTable <- downloadHandler(
    filename = function() {
      paste0(input$ga_ids, "_", input$profile, "_pancan_CPTAC.csv")
    },
    content = function(file) {
      write.csv(plot_func()$data, file, row.names = FALSE)
    }
  )

  output$download <- downloadHandler(
    filename = function() {
      paste0(input$ga_ids,  "_pancan_CPTAC.pdf")
    },
    content = function(file) {
      p <- plot_func()

        pdf(file, width =  input$width_scatter/70 ,height = input$height_scatter/70)
        print(p)
        dev.off()
    }
  )

  # observeEvent(input$search_bttn, {
  #   if (nchar(input$ga_ids) >= 1) {
  #     shinyjs::show(id = "save_csv")
  #   } else {
  #     shinyjs::hide(id = "save_csv")
  #   }
  # })

  output$tbl <- DT::renderDataTable(server = FALSE, {
    DT::datatable(
      plot_func()$data,
      rownames = T,
      extensions = c("Buttons"),
      options = list(
        pageLength = 5,
        dom = "Bfrtip",
        buttons = list(
          list(
            extend = "csv", text = "Download table", filename =  paste0(input$ga_ids,  "_pancan_CPTAC"),
            exportOptions = list(
              modifier = list(page = "all")
            )
          )
        )
      )
    )
  })

}
