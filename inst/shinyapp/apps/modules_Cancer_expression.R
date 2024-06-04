ui.modules_Cancer_expression <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "",
    icon = icon("database"),
    sidebarLayout(
      sidebarPanel(
        width = 2,
        br(),
        shinyWidgets::prettyCheckboxGroup(ns("data_type"), "Select data type:",
                                          choiceNames = c("Proteome","Phosphoproteome","Transcriptome"),
                                          choiceValues = c("Proteome","Phosphoproteome","Transcriptome"),
                                          selected = "Proteome",
                                          shape = "round",
                                          status = "success",
                                          animation = "tada"
        ),
        uiOutput(ns("cohorts_text")),
        br(),
        uiOutput(ns("datasets_text")),
        hr(),
        hr(),
        selectizeInput(
          inputId = ns("ga_id"), # molecule identifier
          label = "Input gene symbol:",
          choices = c("TP53"),
          options = list(
            create = TRUE,
            maxOptions = 5,
            placeholder = "e.g. TP53",
            plugins = list("restore_on_backspace")
          )
        ),
        br(),
        shinyjs::hidden(
          selectizeInput(
            inputId = ns("phoso_site"), # molecule identifier
            label = "Select Phosphorylation Site:",
            choices=c("NP_000537.3:s315","NP_000537.3:s314","NP_000537.3:s392","NP_000537.3:y163s166t170","NP_000537.3:t150"),
            options = list(
              create = TRUE,
              maxOptions = 5,
              placeholder = "e.g. TP53",
              plugins = list("restore_on_backspace")
            )
          )
        ),

        h4("Tips"),
        p(" selected datasets for expression/correlation/survival analysis")


      ),
      mainPanel(
        width = 9,
        tabsetPanel(type = "tabs",id=ns("tcga_single"),
                    tabPanel("Datasets",value = "Datasets",
                               DT::dataTableOutput(ns("xena_table")),
                               # use_waiter(),
                               hr(),
                               fluidRow(
                                 column(
                                   offset = 1,
                                   3,
                                   actionBttn(
                                     inputId = ns("gene_expr"),
                                     label = "Gene expression",
                                     style = "gradient",
                                     # color = "danger",
                                     icon = icon("database"),
                                     size = "sm"
                                   )
                                 ),
                                 column(
                                   3,
                                   actionBttn(
                                     inputId = ns("gene_surv"),
                                     label = "Survival analysis",
                                     style = "gradient",
                                     icon = icon("database"),
                                     size = "sm"
                                   )
                                 ),
                                 column(
                                   3,
                                   actionBttn(
                                     inputId = ns("stage"),
                                     label = "Clinic data",
                                     style = "gradient",
                                     icon = icon("database"),
                                     size = "sm"
                                   )
                                 )

                               )
                               ),
                    tabPanel("Expression analysis",value = "expression",
                             fluidPage(
                               fluidRow(
                                 column(
                                   4,
                                   wellPanel(
                                     h4("Analysis Controls"),
                                     selectInput(inputId = ns("method"), label = "Select statistic method", choices = c("wilcox.test", "t.test"), selected = "t.test"),

                                     materialSwitch(ns("pdist_show_p_value"), "Show P value", inline = FALSE),
                                     materialSwitch(ns("pdist_show_p_label"), "Show P label", inline = FALSE),
                                     materialSwitch(ns("show_n"), "Show sample size", inline = FALSE),
                                     colourpicker::colourInput(inputId = ns("tumor_col"), "Tumor sample color", "#00AFBB"),
                                     colourpicker::colourInput(inputId = ns("normal_col"), "Normal sample color", "#FC4E07"),
                                     selectInput(inputId = ns("theme"), label = "Select theme for plot", choices = names(themes_list), selected = "Base"),

                                     actionBttn(
                                       inputId = ns("ga_submit"),
                                       label = "Submit",
                                       style = "gradient",
                                       icon = icon("check"),
                                       color = "default",
                                       block = TRUE,
                                       size = "sm"
                                     )),
                                   tags$br(),
                                   wellPanel(
                                     p("Download figure:"),
                                     numericInput(inputId = ns("height_scatter"), label = "Height", value = 600),
                                     numericInput(inputId = ns("width_scatter"), label = "Width", value = 600),
                                     hr(),
                                     downloadBttn(
                                       outputId = ns("download"),
                                       style = "gradient",
                                       color = "default",
                                       block = TRUE,
                                       size = "sm"
                                     )
                                   )
                               ),
                                 column(
                                   8,
                                   shinycssloaders::withSpinner(plotOutput(ns("ga_expr_output"),width = "100%",height = "auto")),

                                   DT::dataTableOutput(ns("ga_expr_data"))
                                 )

                             ))
                             ),

                    tabPanel("survival analysis",value = "survival",
                             fluidPage(
                               fluidRow(
                                 bs4Dash::column(
                                   4,
                                   h4("Analysis Controls"),
                                   wellPanel(
                                     shinyWidgets::prettyRadioButtons(
                                       inputId = ns("cutoff_mode"),
                                       label = "Cutoff mode",
                                       choices = c("Auto", "Custom"),
                                       inline = TRUE,
                                       icon = icon("check"),
                                       animation = "jelly"
                                     ),
                                     conditionalPanel(ns=ns,
                                                      condition = "input.cutoff_mode == 'Custom'",
                                                      sliderInput(
                                                        inputId = ns("cutpoint"), label = "Cutoff (%)",
                                                        min = 25, max = 75, value = c(50)
                                                      )),
                                     selectInput(ns("color_palette"), "Color palette:",
                                                 choices = c("npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons", "rickandmorty", "custom"),
                                                 selected = "npg"
                                     ),
                                     checkboxInput(
                                       ns("pval"),
                                       "Show P value?", TRUE),
                                     checkboxInput(
                                       ns("risk.table"),
                                       "Show risk table?", TRUE),
                                     checkboxInput(
                                       ns("conf.int"),
                                       "Show confidence?", TRUE),
                                     selectInput(
                                       ns("surv.median.line"),
                                       "Show surv median line:",
                                       c("none", "hv", "h", "v"),
                                       selected = "none",
                                       multiple = FALSE,
                                       selectize = TRUE,
                                       width = NULL,
                                       size = NULL
                                     ),
                                     selectInput(
                                       ns("legend.pos"),
                                       "Legend pos:",
                                       c("left","top","bottom","right","none"),
                                       selected = "top",
                                       multiple = FALSE,
                                       selectize = TRUE,
                                       width = NULL,
                                       size = NULL
                                     )
                                     ),

                                   actionBttn(
                                     inputId = ns("surv_go"),
                                     label = "GO!",
                                     style = "gradient",
                                     icon = icon("check"),
                                     color = "default",
                                     block = TRUE,
                                     size = "sm"
                                   ),

                                   numericInput(inputId = ns("height_surv"), label = "Height", value = 450),
                                   numericInput(inputId = ns("width_surv"), label = "Width", value = 600),
                                   downloadBttn(
                                     outputId = ns("km_download"),
                                     style = "gradient",
                                     color = "default",
                                     block = TRUE,
                                     size = "sm"
                                   )

                                 ),
                                 bs4Dash::column(
                                   8,
                                   h4("Note:"),
                                   p("Only BRCA_TCGA, CCRCC_CPTAC, GBM_Pediatric, HNSCC_CPTAC, LSCC_CPTAC, LUAD_CPTAC, OV_JHU, OV_PNNL, PDAC_CPTAC, PDAC_Enrichment, PDAC_KU, UCEC_CPTAC1, UCEC_CPTAC2 cohorts have OS data!!!"),
                                   hr(),
                                   shinyjs::useShinyjs(),  # Set up shinyjs

                                   shinycssloaders::withSpinner(plotOutput(ns("KM_plot"),width = "100%",height = "auto")),
                                   br(),
                                   DT::dataTableOutput(ns("ga_output_data"))
                                   ),
                                 bs4Dash::column(
                                   8,align="center",

                                 )
                               )
                             )
                             ),
                    tabPanel("Clinic data",value = "clinic",
                             fluidPage(
                               fluidRow(
                                 column(
                                   4,
                                   wellPanel(
                                     h4("Analysis Controls"),
                                     actionBttn(
                                       inputId = ns("ga_stage_submit"),
                                       label = "Fetch clinic data!",
                                       style = "gradient",
                                       icon = icon("check"),
                                       color = "default",
                                       block = TRUE,
                                       size = "sm"
                                     ),
                                     hr(),
                                     selectInput(inputId = ns("feather"), label = "Select clinic feather to plot", choices = "", selected = ""),
                                     prettyCheckbox(ns("iscont"),"Convert continuous variables to categorical variables?",F),
                                     conditionalPanel(ns=ns,
                                                      condition = "input.iscont == true",
                                                      numericInput(ns("cont.cut"),"Cut point:",value = 0)
                                     ),
                                     hr(),
                                     selectInput(inputId = ns("theme_stage"), label = "Select theme for plot", choices = names(themes_list), selected = "Base"),
                                     hr()
                                     ),
                                   tags$br(),
                                   wellPanel(
                                     p("Download figure:"),
                                     numericInput(inputId = ns("height_stage"), label = "Height", value = 400),
                                     numericInput(inputId = ns("width_stage"), label = "Width", value = 600),
                                     hr(),
                                     downloadBttn(
                                       outputId = ns("download_stage"),
                                       style = "gradient",
                                       color = "default",
                                       block = TRUE,
                                       size = "sm"
                                     )
                                   )
                                 ),
                                 column(
                                   8,
                                   bs4Dash::tabBox(
                                     id = ns("survtab"), width = 12,

                                     tabPanel( title = "Clinic data",
                                              DT::dataTableOutput(ns("ga_stage_data"))

                                     ),
                                     tabPanel( title = "Stage plot",
                                              plotOutput(ns("ga_stage_output"),width = "100%",height = "auto")

                                     )
                                   )
                                 )

                               ))
                    )

        )

      )
    )
  )
}

server.modules_Cancer_expression <- function(input, output, session) {
  ns <- session$ns
  # df <- dataset
  OS_dataset <- c("BRCA_TCGA","CCRCC_CPTAC","GBM_Pediatric","HNSCC_CPTAC","LSCC_CPTAC","LUAD_CPTAC","OV_JHU","OV_PNNL","PDAC_CPTAC","PDAC_Enrichment","PDAC_KU","UCEC_CPTAC1","UCEC_CPTAC2")
  observe({

    updateSelectizeInput(
      session,
      "ga_id",
      choices = ID_list_pro %>% unlist() %>% unique(),
      selected = "TP53",
      server = TRUE
    )
  })
  observeEvent(input$datasets_text,{
    if (stringr::str_detect( input$datasets_text,"_Phospho")) {
      shinyjs::show(id = "phoso_site")
    }else{
      shinyjs::hide(id = "phoso_site")

    }
  })
  output$cohorts_text <- renderUI({
    selectInput(
      inputId = ns("cohorts_text"),
      label = "Cancer type:",
      choices = unique(dataset_info %>% dplyr::filter(Data.type %in% input$data_type) %>% .["Primary.Site"]),
      selected = "Lung",
      multiple = T
    )
  })
  output$datasets_text <- renderUI({
    selectInput(
      inputId = ns("datasets_text"),
      label = "Select dataset:",
      choices = dataset_info %>% dplyr::filter(Data.type %in% input$data_type)%>%
        dplyr::filter(Primary.Site %in% input$cohorts_text)%>% .["Abbre"],
      multiple = FALSE,
      selected = dataset_info %>% dplyr::filter(Data.type %in% input$data_type)%>%
        dplyr::filter(Primary.Site %in% input$cohorts_text)%>% .["Abbre"] %>%.[1]
    )
  })
  output$datasets_text2 <- renderUI({
    selectInput(
      inputId = ns("datasets_text2"),
      label = "Select dataset:",
      choices = dataset_info$Abbre,
      selected = NULL,
      multiple = FALSE
    )
  })


  output$xena_table <- DT::renderDataTable({
    dataset_info %>% dplyr::filter(Data.type %in% input$data_type)%>%
      dplyr::filter(Primary.Site %in% input$cohorts_text)%>%
      DT::datatable(
        rownames = FALSE,selection = 'single',
        options = list(
          language = list(search = "Filter with keyword:"),
          initComplete = JS(
            "function(settings, json) {",
            "$(this.api().table().header()).css({'background-color': '#4a5a6a', 'color': '#fff'});",
            "$(this.api().table().header()).css({'font-size': '95%'});",
            "$(this.api().table().body()).css({'font-size': '90%'});",
            "}"
          )
        )
      )
  })

  observeEvent(input$ga_id,
               updateSelectInput(session, "phoso_site",
                                 label = "Select Phosphorylation Site:",
                                 choices =  idmap_protein[which(idmap_protein$Symbol == input$ga_id),]$row_names
               )
  )

  observeEvent(input$gene_expr,{
                 updateTabsetPanel(session, "tcga_single",
                                   selected = "expression"
                 )})
  observeEvent(input$gene_surv,{
                    updateTabsetPanel(session, "tcga_single",
                                      selected = "survival"
                    )})
  observeEvent(input$stage,{
    updateTabsetPanel(session, "tcga_single",
                      selected = "clinic"
    )})

  # Keep selected database
  selected_database <- reactive({
    s <- input$xena_table_rows_selected
    if (length(s)) {
      return(dataset_info[s, ])
    }
  })
  observeEvent(input$xena_table_rows_selected,{
    s <- input$xena_table_rows_selected
    if (length(s)) {
      ddd <- dataset_info %>% dplyr::filter(Data.type %in% input$data_type)%>%
        dplyr::filter(Primary.Site %in% input$cohorts_text)%>% .[s,]
    }
    updateSelectInput(session, "datasets_text",
                      label = "Select dataset:",
                      selected = ddd$Abbre
    )



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
  colors <- reactive({
    c(input$tumor_col, input$normal_col)
  })


  # Show waiter for plot
  w <- waiter::Waiter$new(id = ns("ga_expr_output"), html = waiter::spin_hexdots(), color = "white")
  plot_theme <- reactive({
    themes_list[[input$theme]]
  })
  plot_func <- eventReactive(input$ga_submit, {
    datasets_text <- input$datasets_text
    id <- input$ga_id
    if (stringr::str_detect(datasets_text,"Phospho")) {
      id <- input$phoso_site
    }
    data_select <- get_expr_data(datasets_text, id)
    if (is.null(data_select)) {
      sendSweetAlert(
        session,
        title = "Error",
        text = "Error to query data and plot. Please make sure you have typed corrected gene symbol.",
        type = "error"
      )
      return(NULL)
    }
    p<- viz_TvsN(data_select,df_type = "single",
                   Method =  input$method,
                   Show.P.value = input$pdist_show_p_value,
                   Show.P.label = input$pdist_show_p_label,
                   Show.n = input$show_n,
                   values = colors())  +
      ggplot2::guides(fill = ggplot2::guide_legend(title = NULL)) +
      ylab(ifelse(stringr::str_detect(input$datasets_text,"Phospho"),input$phoso_site,input$ga_id))
    names(data_select)[1]<-id
    p_df<-list("data_select" = data_select,"p"=p)
    return(p_df)
  })
  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })
  output$ga_expr_output <- renderPlot(width = width_scatter,
                                      height = height_scatter,{
                                        w$show() # Waiter add-ins
                                          plot_func()$p+ plot_theme()+
                                            ggplot2::theme(
                                              axis.text.x = element_text(size = 16),
                                              axis.text.y = element_text(size = 16),
                                              axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16)
                                            )+
                                            ggplot2::theme(axis.text.x = element_text(hjust = .5, vjust = .5, size = 16),
                                                           axis.text.y = element_text( size = 16))
                                        }


                                      )
  output$download <- downloadHandler(
    filename = function() {
      paste0(input$ga_id, "_",input$datasets_text,".pdf")
    },
    content = function(file) {
        pdf(file = file, onefile = FALSE, width = input$width_scatter/70 ,height = input$height_scatter/70)
        print(plot_func()$p+ plot_theme()+
                ggplot2::theme(
                  axis.text.x = element_text(size = 16),
                  axis.text.y = element_text(size = 16),
                  axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16)
                )+
                ggplot2::theme(axis.text.x = element_text(hjust = .5, vjust = .5, size = 16),
                               axis.text.y = element_text( size = 16))+
                ylab(ifelse(stringr::str_detect(input$datasets_text,"Phospho"),input$phoso_site,input$ga_id))
        )
        dev.off()

      # ggplot2::ggsave(
      #   filename = file, plot = print(p_surv(), newpage = F), device = input$device,
      #   units = "cm", width = input$width, height = input$height, dpi = 600
      # )
    }
  )
  br()
  br()
  output$ga_expr_data <- DT::renderDataTable(server = FALSE, {
      DT::datatable(
        plot_func()$data_select,
        rownames = T,
        extensions = c("Buttons"),
        options = list(
          pageLength = 5,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table", filename = paste(input$ga_id, "in",input$datasets_text),
              exportOptions = list(
                modifier = list(page = "all")
              )
            )
          )
        )
      )
  })
  ####survival
  surv_data <- eventReactive(input$surv_go,{
    datasets_text <- input$datasets_text
    if (!stringr::str_remove(datasets_text,"_Phospho|_protein|_mRNA") %in% OS_dataset){
      sendSweetAlert(
        session,
        title = "Error",
        text = "Please note that only BRCA_TCGA, CCRCC_CPTAC, GBM_Pediatric, HNSCC_CPTAC, LSCC_CPTAC, LUAD_CPTAC, OV_JHU, OV_PNNL, PDAC_CPTAC, PDAC_Enrichment, PDAC_KU, UCEC_CPTAC1, UCEC_CPTAC2 cohorts have OS data!!!",
        type = "error"
      )
      return(NULL)
    }else{
      clinic <- get_data(stringr::str_remove(datasets_text,"_Phospho|_protein|_mRNA"),
                         "clinic")
      OS_data <- clinic[c("Cases_Submitter_ID","OS.time", "OS.status")]
      colnames(OS_data)[1] <- "ID"
    OS_data$OS.time <- ifelse(OS_data$OS.time < 0,0,OS_data$OS.time)
      expr_data <- get_expr_data(datasets_text, input$ga_id) %>% na.omit()

    if (is.null(expr_data)) {
      sendSweetAlert(
        session,
        title = "Error",
        text = "Error to query data and plot. Please make sure you have typed corrected gene symbol.",
        type = "error"
      )
      return(NULL)
    }else{

      expr <- expr_data  %>%
        dplyr::filter(type == "Tumor")  %>%
        inner_join(.,OS_data,by = "ID") %>%
        mutate(OS.time = as.numeric( OS.time)/365) %>%
        mutate(OS.status = as.numeric(OS.status))
      colnames(expr)[4] <- "value"
      if (input$cutoff_mode == "Auto") {
        expr2 <- expr %>%
          survminer::surv_cutpoint(
            time = "OS.time", event = "OS.status",
            variables = c("value"),
            minprop = 0.25, progressbar = TRUE
          ) %>%
          survminer::surv_categorize(labels = c("Low", "High")) %>%
          data.frame()
        expr$group <- expr2$value
      }else{
        expr <- expr %>%
          dplyr::arrange(.data$value) %>%
          dplyr::mutate(per_rank = 100 / nrow(.) * (1:nrow(.)))
        expr <- expr %>%
          dplyr::mutate(group = dplyr::case_when(
            .data$per_rank > input$cutpoint ~ "High",
            .data$per_rank <= input$cutpoint ~ "Low",
            TRUE ~ NA_character_
          ))
      }
      # expr$group <- factor(expr$group,levels = c("Low","High"))
      expr <- na.omit(expr)
      if (nrow(expr) < 10){
        sendSweetAlert(
          session,
          title = "warning ",
          text = "The valid sample size is too small (< 10), try other genes or sites.",
          type = "warning"
        )
      }
      return(expr)
    }
    }
  })
  output$ga_output_data <- DT::renderDataTable(server = FALSE, {
    DT::datatable(
      surv_data(),
      rownames = FALSE,
      extensions = c("Buttons"),
      options = list(
        pageLength = 5,
        dom = "Bfrtip",
        buttons = list(
          list(
            extend = "csv", text = "Download table", filename = paste(input$ga_id, "_surv_data_",input$datasets_text),
            exportOptions = list(
              modifier = list(page = "all")
            )
          )
        )
      )
    )
  })
  observeEvent(input$surv_go,{
    if (!is.null(surv_data())){
      shinyjs::show(id = "km_download")
    }else {
      shinyjs::hide(id = "km_download")
    }
  })

  # Show waiter for plot
  w <- waiter::Waiter$new(id = "wb_plot0", html = waiter::spin_hexdots(), color = "black")
  width_surv <- reactive ({input$width_surv })
  height_surv <- reactive ({input$height_surv })

  output$KM_plot <- renderPlot(width = width_surv,
                               height = height_surv,
                               {
    w$show() # Waiter add-ins
      if (is.null(surv_data())){
        NULL
      }else if (nrow(surv_data()) < 10){
        NULL
      }else{
        p_survplot(surv_data(),palette = input$color_palette, pval = input$pval, risk.table=input$risk.table,conf.int = input$conf.int,size=1,surv.median.line = input$surv.median.line,xlab = "Time (year)",
                   legend.title=input$ga_id,legend.pos= input$legend.pos)
      }



  })



  output$km_download <- downloadHandler(
    filename = function() {
      paste(input$ga_id, "_surv_data_",input$datasets_text,"_KMplot.pdf")
    },
    content = function(file) {
       p <-  p_survplot(surv_data(),palette = input$color_palette, pval = input$pval, risk.table=input$risk.table,conf.int = input$conf.int,size=1,surv.median.line = input$surv.median.line,xlab = "Time (year)",
                 legend.title=surv_data()[[2]],legend.pos= input$legend.pos)

      pdf(file,onefile = T, width = input$width_surv/70 ,height = input$height_surv/70)
      print(p)
      dev.off()
    })
  ####clinical
  clinic_merge <- eventReactive(input$ga_stage_submit,{
    data_select <- get_expr_data(input$datasets_text, input$ga_id)
    if (is.null(data_select)) {
      sendSweetAlert(
        session,
        title = "Error",
        text = "Error to query data and plot. Please make sure you have typed corrected gene symbol.",
        type = "error"
      )
      return(NULL)
    }
    dd <- merge_clinic_data(stringr::str_remove(input$datasets_text,"_Phospho|_protein|_mRNA"),data_select)
    return(dd)
  })
  plot_theme_stage <- reactive({
    themes_list[[input$theme_stage]]
  })

  output$ga_stage_data <- DT::renderDataTable(server = FALSE, {
    DT::datatable(
      clinic_merge(),
      rownames = FALSE,
      extensions = c("Buttons"),
      options = list(
        pageLength = 10,
        dom = "Bfrtip",
        buttons = list(
          list(
            extend = "csv", text = "Download table", filename = paste(input$ga_id, "_clinic_data_",input$datasets_text),
            exportOptions = list(
              modifier = list(page = "all")
            )
          )
        ),
        scrollX = TRUE
      )
    )
  })
  observeEvent(input$ga_stage_submit,{
    updateSelectInput(session, "feather",
                      label = "Select clinic feather to plot",
                      choices =colnames(clinic_merge())[-1:-4],selected = "Gender"
                      )
  })

  observeEvent({input$feather
    input$iscont},{
    if ( input$iscont ==T){
      vv <- clinic_merge() %>% .[,input$feather] %>% as.numeric() %>% na.omit()
      print(vv)
      updateNumericInput(session,
                         "cont.cut",
                         value = median(vv))
    }

  })
  width_stage <- reactive ({ input$width_stage })
  height_stage <- reactive ({ input$height_stage })

 output$ga_stage_output <- renderPlot(width = width_stage,
                                      height = height_stage,
                                       {
                                           w$show() # Waiter add-ins

                                           stage_plot(clinic_merge(),input$feather,is.cont = input$iscont,cont.cut = input$cont.cut)[["p"]]+ plot_theme_stage()+
                                             ggplot2::theme(
                                               axis.text.x = element_text(size = 20),
                                               axis.text.y = element_text(size = 20),
                                               axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20)
                                             )+
                                             ggplot2::theme(axis.text.x = element_text(hjust = .5, vjust = .5, size = 20),
                                                            axis.text.y = element_text( size = 20))+
                                             ylab(ifelse(stringr::str_detect(input$datasets_text,"Phospho"),input$phoso_site,input$ga_id))



                                       })
 output$download_stage <- downloadHandler(
   filename = function() {
     paste(input$ga_id, "_",input$feather,"_",input$datasets_text,".pdf")
   },
   content = function(file) {
     p <-   stage_plot(clinic_merge(),input$feather,is.cont = input$iscont,cont.cut = input$cont.cut)[["p"]]+ plot_theme_stage()+
       ggplot2::theme(
         axis.text.x = element_text(size = 20),
         axis.text.y = element_text(size = 20),
         axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20)
       )+
       ggplot2::theme(axis.text.x = element_text(hjust = .5, vjust = .5, size = 20),
                      axis.text.y = element_text( size = 20))+
       ylab(ifelse(stringr::str_detect(input$datasets_text,"Phospho"),input$phoso_site,input$ga_id))

     pdf(file,onefile = T, width = input$width_stage/70 ,height = input$height_stage/70)
     print(p)
     dev.off()
   })

  }
