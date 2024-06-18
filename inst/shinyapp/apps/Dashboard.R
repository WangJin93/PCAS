ui.modules_dash <- function(id) {
  ns <- NS(id)
  tagList(
    bs4Dash::bs4Jumbotron(
      title = "Welcome to the ProteoCancer Analysis Suite (PCAS)!",
      lead = HTML("The ProteoCancer Analysis Suite (PCAS) platform takes advantage of the CPTAC database and offers an integrated platform for proteomics, phosphoproteomics, and transcriptomics analysis in cancer research. It includes modules for analyzing individual data sets, such as gene expression, correlation, survival analysis, and clinical information. Additionally, the app allows for multi-dataset analyses across different cancer types, providing a comprehensive view of cancer biology. Advanced features like immune infiltration and drug sensitivity analysis offer insights into the tumor microenvironment. By simplifying complex data analysis, our app enables researchers to advance oncology studies and develop targeted therapies.</br></br>Install PCAS packages: <b>remotes::install_github('WangJin93/PCAS')</b>, Then execute commands: <b>PCAS::PCAS_app() </b> to run the app locally.</br><b>Citation: </b>Wang J, Song X, Wei M, Qin L, Zhu Q, Wang S, Liang T, Hu W, Zhu X, Li J. PCAS: An Integrated Tool for Multi-Dimensional Cancer Research Utilizing Clinical Proteomic Tumor Analysis Consortium Data. International Journal of Molecular Sciences. 2024; 25(12):6690. https://doi.org/10.3390/ijms25126690"),
      status = "info",
      btnName = "Github source code",
      href = "https://github.com/WangJin93/PCAS"
    ),
    fluidRow(
      bs4Dash::column(6,
                      bs4Dash::bs4UserCard(
               title = bs4UserDescription(
                 title = "Application Developer",
                 subtitle = "Jin wang",
                 image = "P1102513.jpg",
                 type = 1
               ),
               status = "info",
               width = 12,

               bs4Dash::bs4ListGroup(
                 width = 12,
                 type = "action",
                 bs4ListGroupItem(
                   h4("Email: Jinwang93@suda.edu.cn")
                 ),
                 bs4ListGroupItem(
                   h4("Affiliation: Soochow University")

                 ),
                 bs4ListGroupItem(
                   h4("Personal website: https://www.jingege.wang"),
                   href = "https://www.jingege.wang"
                 )
               )
             )
      ),
      column(6,
             tags$img(src = 'abstract.png', width = '100%')
             #
             # shinycssloaders::withSpinner(DTOutput(outputId = ns("tf_info"))),
             #   downloadButton(ns("download_tf.csv"), "Download csv table",class = "mybutton")

      )

    )
  )
}
server.modules_dash <- function(input, output, session) {
  ns <- session$ns

  # output$tf_info <- renderDT({
  #
  #   DT::datatable(
  #     tf_list,rownames = F,
  #     options = list(pageLength = 5)
  #   )
  #
  # })
  #
  # output$download_tf.csv <- downloadHandler(
  #   filename = function() {
  #     "TF_list.csv"
  #   },
  #   content = function(file) {
  #     write.csv(tf_list, file, row.names = FALSE)
  #   }
  # )

}

