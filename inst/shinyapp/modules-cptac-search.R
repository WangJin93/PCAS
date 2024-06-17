ui.modules_pancptac_dist <- function(id) {
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
        selectInput(inputId = ns("method"), label = "Select statistic method", choices = c("wilcox.test", "t.test"), selected = "t.test"),
        materialSwitch(ns("pdist_show_p_value"), "Show P value", inline = TRUE),
        materialSwitch(ns("pdist_show_p_label"), "Show P label", inline = TRUE),
        materialSwitch(ns("show_n"), "Show sample size", inline = FALSE),
        numericInput(inputId = ns("show_n_position"), label = "Show n position", value = -3),
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

             tags$br(),
             DT::DTOutput(outputId = ns("tbl")),
      )
    )
  )
}

server.modules_pancptac_dist <- function(input, output, session) {
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
  observeEvent(eventExpr = input$data_type,{
    if (input$data_type =="Phosphoproteome") {
      shinyjs::show(id = "phoso_site")
      }else{
        shinyjs::hide(id = "phoso_site")

      }
  })
  count_id <- function(id){
  }
  observeEvent({input$Pancan_search
    input$datasets_text},{
      ids <- idmap_protein[which(idmap_protein$Symbol == input$Pancan_search),]$row_names
      ids <- ids[which(lapply(ids, function(x){
        sum(sapply(phoso_id[input$datasets_text], `%in%`, x=x))
      }) > 1)]

      updateSelectInput(session, "phoso_site",
                        label = "Select Phosphorylation Site:",
                        choices =ids)
    }

  )

  colors <- reactive({
    c(input$tumor_col, input$normal_col)
  })


  # Show waiter for plot
  w <- waiter::Waiter$new(id = ns("gene_pancan_dist"), html = waiter::spin_hexdots(), color = "black")

  plot_func <- eventReactive(input$search_bttn, {
    if (nchar(input$Pancan_search) >= 1) {
      Pancan_search <- input$Pancan_search
      if (input$data_type == "Phosphoproteome") {
        Pancan_search <- input$phoso_site
      }
      df <- get_expr_data(datasets=input$datasets_text,
                           genes = Pancan_search)
      p <- viz_TvsN(df,df_type = "multi_set",
        Method =  input$method,
        Show.P.value = input$pdist_show_p_value,
        Show.P.label = input$pdist_show_p_label,
        Show.n = input$show_n,
        Show.n.location =input$show_n_position,
        values = colors()
      )
        p <- p+ ylab(ifelse(input$data_type == "Phosphoproteome", Pancan_search,paste0(input$Pancan_search, " expression")))

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
    plot_func() + ggplot2::theme(
        plot.margin = margin(t = 0,
                             r = 0,
                             b = 50,
                             l = 50)
      )
  })


  output$download <- downloadHandler(
    filename = function() {
      paste0(input$Pancan_search,  "_pancan_CPTAC.pdf")
    },
    content = function(file) {
      p<- plot_func() + ggplot2::theme(
        plot.margin = margin(t = 0,
                             r = 0,
                             b = 50,
                             l = 50)
      )
      pdf(file, width =  input$width_scatter/70 ,height = input$height_scatter/70)
        print(p)
        dev.off()
    }
  )

  observeEvent(input$search_bttn, {
    if (nchar(input$Pancan_search) >= 1) {
      shinyjs::show(id = "save_csv")
    } else {
      shinyjs::hide(id = "save_csv")
    }
  })

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
            extend = "csv", text = "Download table", filename =  paste0(input$Pancan_search,  "_pancan_CPTAC"),
            exportOptions = list(
              modifier = list(page = "all")
            )
          )
        )
      )
    )
  })

}
