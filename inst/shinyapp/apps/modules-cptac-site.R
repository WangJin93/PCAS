ui.modules_cptac_site <- function(id) {
  ns <- NS(id)
  fluidPage(
fluidRow(column(3,
                wellPanel(
                  hr(),
                  h3("Visualization of protein structure and phosphorylation sites"),
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
        selectInput(inputId = ns("method"), label = "Source of phosphorylation site information", choices = c("UniProt", "CPTAC"), selected = "CPTAC"),
             tags$hr(style = "border:none; border-top:2px solid #5E81AC;"),
        hr(),
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
          numericInput(inputId = ns("height_scatter"), label = "Height", value = 400),
          numericInput(inputId = ns("width_scatter"), label = "Width", value = 1000,max = 1200),
        downloadBttn(
          outputId = ns("download"),
          style = "gradient",
          color = "default",
          block = TRUE,
          size = "sm"
        ))
      ),
      column(9,
             shinycssloaders::withSpinner(plotOutput(ns("phoso_dist"), height = "auto")),
             hr(),

             tags$br(),
             DT::DTOutput(outputId = ns("tbl")),
      )
    )
  )
}

server.modules_cptac_site <- function(input, output, session) {
  ns <- session$ns
  observe({

    updateSelectizeInput(
      session,
      "Pancan_search",
      choices = ID_list_pro %>% unlist() %>% unique(),
      selected = "TNS1",
      server = TRUE
    )
  })



  # Show waiter for plot
  w <- waiter::Waiter$new(id = ns("gene_pancan_dist"), html = waiter::spin_hexdots(), color = "black")

  plot_func <- eventReactive(input$search_bttn, {
    if (nchar(input$Pancan_search) >= 1) {
      Pancan_search <- input$Pancan_search
      p <- viz_phoso_sites(
        gene = Pancan_search,phoso_infoDB = input$method
      )
    }
    return(p)
  })


  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })
  output$phoso_dist <- renderPlot(width = width_scatter,
                                        height = height_scatter,{
    w$show() # Waiter add-ins
    plot_func()
  })


  output$download <- downloadHandler(
    filename = function() {
      paste0(input$Pancan_search,  "_phoso_site.pdf")
    },
    content = function(file) {
      p <- plot_func()
        pdf(file, width =  input$width_scatter/70 ,height = input$height_scatter/70)
        print(p)
        dev.off()
    }
  )


  output$tbl <- DT::renderDataTable(server = FALSE, {
    DT::datatable(
      plot_func()[["layers"]][[3]][["data"]],
      rownames = T,
      extensions = c("Buttons"),
      options = list(
        pageLength = 5,
        dom = "Bfrtip",
        buttons = list(
          list(
            extend = "csv", text = "Download table", filename =  paste0(input$Pancan_search,  "_phoso_site"),
            exportOptions = list(
              modifier = list(page = "all")
            )
          )
        )
      )
    )
  })

}
