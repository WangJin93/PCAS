#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(PCAS)
# 安装并加载CRAN包
if (!require("shiny")) install.packages("shiny", update = FALSE, ask = FALSE); library(shiny)
if (!require("waiter")) install.packages("waiter", update = FALSE, ask = FALSE); library(waiter)
if (!require("ggplot2")) install.packages("ggplot2", update = FALSE, ask = FALSE); library(ggplot2)
if (!require("ggpubr")) install.packages("ggpubr", update = FALSE, ask = FALSE); library(ggpubr)
if (!require("dplyr")) install.packages("dplyr", update = FALSE, ask = FALSE); library(dplyr)
if (!require("DT")) install.packages("DT", update = FALSE, ask = FALSE); library(DT)
if (!require("shinyalert")) install.packages("shinyalert", update = FALSE, ask = FALSE); library(shinyalert)
if (!require("stringr")) install.packages("stringr", update = FALSE, ask = FALSE); library(stringr)
if (!require("shinyWidgets")) install.packages("shinyWidgets", update = FALSE, ask = FALSE); library(shinyWidgets)
if (!require("psych")) install.packages("psych", update = FALSE, ask = FALSE); library(psych)
if (!require("ggthemes")) install.packages("ggthemes", update = FALSE, ask = FALSE); library(ggthemes)
if (!require("shinyjs")) install.packages("shinyjs", update = FALSE, ask = FALSE); library(shinyjs)
if (!require("bs4Dash")) install.packages("bs4Dash", update = FALSE, ask = FALSE); library(bs4Dash)
if (!require("ggrepel")) install.packages("ggrepel", update = FALSE, ask = FALSE); library(ggrepel)
if (!require("survival")) install.packages("survival", update = FALSE, ask = FALSE); library(survival)
if (!require("tibble")) install.packages("tibble", update = FALSE, ask = FALSE); library(tibble)
if (!require("survminer")) install.packages("survminer", update = FALSE, ask = FALSE); library(survminer)
if (!require("tidyr")) install.packages("tidyr", update = FALSE, ask = FALSE); library(tidyr)

# 安装并加载Bioconductor包
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", update = FALSE, ask = FALSE)

if (!require("aplot")) BiocManager::install("aplot", update = FALSE, ask = FALSE); library(aplot)
if (!require("ggtree")) BiocManager::install("ggtree", update = FALSE, ask = FALSE); library(ggtree)
if (!require("drawProteins")) BiocManager::install("drawProteins", update = FALSE, ask = FALSE); library(drawProteins)

source("apps/visualize.R")
source("apps/Dashboard.R")
source("apps/modules_Cancer_expression.R")
source("apps/modules_Cancer_DEGs.R")
source("apps/modules_Cancer_correlation.R")
source("apps/modules-cptac-search.R")
source("apps/modules_Cancer_multiple.R")
source("apps/modules-pancan-til.R")
source("apps/modules-pancan-drug.R")
source("apps/modules-cptac-site.R")
source("apps/modules-pancan-corr.R")
source("apps/mod_feedback.R")
# Define UI for application that draws a histogram
ui <- bs4Dash::dashboardPage(
  bs4Dash::dashboardHeader(
    title = tags$a(href='https://github.com/WangJin93/PCAS',
                   tags$img(src='logo.png',width="100%"),
                   style = "position: relative; top: 5px;",align="center")
  ),

  bs4Dash::dashboardSidebar(
    bs4Dash::sidebarMenu(
      bs4Dash::menuItem("Introduction",icon = icon("info-circle") , tabName = "Welcome"),
      bs4Dash::menuItem("Single Dataset Analysis",icon = icon('palette'),
               menuSubItem("Single Gene Analysis", tabName = "Single"),
               menuSubItem("Multi-Gene expression", tabName = "Multiple"),
               menuSubItem("Correlation Analysis", tabName = "correlation"),
               menuSubItem("DEPs/DEGs Analysis", tabName = "DEGs")
      ),
      bs4Dash::menuItem("Multi-Datasets Analysis",icon = icon("deezer") ,
               menuSubItem("Pancancer expression", tabName = "Pancan_cptac"),
               menuSubItem("Correlation analysis", tabName = "Genelist_corr"),
               menuSubItem("Immune infiltration", tabName = "Immune_infiltration"),
               menuSubItem("Drug sensitivity", tabName = "drug_sensitivity")
               ),
      bs4Dash::menuItem("Phoso-sites visualization",icon = icon('th') , tabName ="phoso_site"),
      bs4Dash::menuItem("Feedback",icon = icon('envelope') , tabName ="feedback")
    )
  ),

  bs4Dash::dashboardBody(
    bs4Dash::tabItems(
      tabItem("Welcome",ui.modules_dash("Welcome")),
      tabItem("Single",ui.modules_Cancer_expression("Single")),
      tabItem("correlation",ui.modules_Cancer_corr("correlation")),
      tabItem("DEGs",ui.modules_diff_gene("DEGs")),
      tabItem("Multiple",ui.modules_multi_gene("Multiple")),
      tabItem("Immune_infiltration",ui.modules_pancan_til("Immune_infiltration")),
      tabItem("Genelist_corr",ui.modules_pancan_corr("Genelist_corr")),
      tabItem("Pancan_cptac",ui.modules_pancptac_dist("module_Pancan_cptac")),
      tabItem("drug_sensitivity",ui.modules_pancan_drug("drug_sensitivity")),
      tabItem("phoso_site",ui.modules_cptac_site("phoso_site")),
      tabItem("feedback",mod_feedback_ui("feedback"))

    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  callModule(server.modules_dash, "Welcome")
  callModule(server.modules_Cancer_expression, "Single")
  callModule(server.modules_Cancer_corr, "correlation")
  callModule(server.modules_diff_gene, "DEGs")
  callModule(server.modules_multi_gene, "Multiple")
  # callModule(server.modules_cptac_gene_cor, "Pancan_correlation")
  callModule(server.modules_pancptac_dist, "module_Pancan_cptac")
  callModule(server.modules_pancan_til, "Immune_infiltration")
  callModule(server.modules_pancan_corr, "Genelist_corr")
  callModule(server.modules_pancan_drug, "drug_sensitivity")
  callModule(server.modules_cptac_site, "phoso_site")
  callModule(mod_feedback_server, "feedback")
}

# Run the application
shinyApp(ui = ui, server = server)
