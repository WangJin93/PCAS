#' feedback UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_feedback_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
             HTML("<h2 style='text-align: center;color: #666'>Feedback</h2>

<hr style='border:3 double #987cb9' width='80%' color='#987cb9' size='10'>"),
        bs4Card(
          title = "Details",
          status = "primary",
          solidHeader = FALSE,
          collapsible = FALSE,
          collapsed = FALSE,
          closable = FALSE,
          label = NULL,
          width = 12,
          tagList(
            fluidRow(
              column(6,
                textInput(
                  ns("issue_title"),
                  "Title",
                  value = "",
                  placeholder = "enter short label",
                  width = "100%"
                )
              ),
              column(6,
                     shinyWidgets::pickerInput(
                  ns("issue_labels"),
                  label = "Choose one category",
                  choices = c("bug", "enhancement"),
                  multiple = FALSE
                )
              )
            ),
            fluidRow(
              column(12,
                textAreaInput(
                  ns("issue_description"),
                  label = "Enter description",
                  value = "",
                  width = "300%",
                  cols = 80,
                  rows = 10,
                  placeholder = "markdown format permitted",
                  resize = "vertical"
                )
              )
            ),
            fluidRow(
              column(2,
                shinyWidgets::actionBttn(
                  ns("submit_issue"),
                  "Submit!",
                  icon = icon("save"),
                  style = "jelly",
                  color = "success",
                  size = "md"
                )
              )
            )
          )
        )
      )
    )
  )
}
    
#' feedback Server Function
#'
#' @noRd 
mod_feedback_server <- function(input, output, session){
  ns <- session$ns
  
  observeEvent(input$submit_issue, {
    # perform checks for mandatory inputs
    if (!shiny::isTruthy(input$issue_title)) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Oops!",
        text = "Please enter a title for your issue.",
        type = "error"
      )
      return(NULL)
    }
    
    if (!shiny::isTruthy(input$issue_labels)) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Oops!",
        text = "Please choose at least one category for your issue.",
        type = "error"
      )
      return(NULL)
    }
    
    if (!shiny::isTruthy(input$issue_description)) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Oops!",
        text = "Please enter a description for your issue.",
        type = "error"
      )
      return(NULL)
    }
    
    # ---- submit issue (credentials come from the environment, never from the
    #      repository; set PCAS_SMTP_USER / PCAS_SMTP_PASS / PCAS_SMTP_FROM /
    #      PCAS_SMTP_TO before deploying, e.g. in ~/.Renviron) -------------
    if (!requireNamespace("mailR", quietly = TRUE)) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Setup error",
        text = "The 'mailR' package is required to send feedback but is not installed.",
        type = "error"
      )
      return(NULL)
    }
    smtp_user <- Sys.getenv("PCAS_SMTP_USER", unset = "")
    smtp_pass <- Sys.getenv("PCAS_SMTP_PASS", unset = "")
    smtp_from <- Sys.getenv("PCAS_SMTP_FROM", unset = smtp_user)
    smtp_to   <- Sys.getenv("PCAS_SMTP_TO",
                            unset = "jin.wang93@outlook.com")
    if (!nzchar(smtp_user) || !nzchar(smtp_pass) || !nzchar(smtp_from)) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Setup error",
        text = "Feedback e-mail is not configured on this server (missing SMTP credentials). Please contact the maintainer directly.",
        type = "error"
      )
      return(NULL)
    }

    ok <- tryCatch({
      mailR::send.mail(
        from    = smtp_from,
        to      = c(smtp_to),
        subject = paste0(input$issue_title, " & ", input$issue_labels),
        body    = input$issue_description,
        smtp    = list(host.name = Sys.getenv("PCAS_SMTP_HOST",
                                              unset = "smtp.qq.com"),
                       port     = as.integer(Sys.getenv("PCAS_SMTP_PORT",
                                                        unset = "465")),
                       user.name = smtp_user,
                       passwd    = smtp_pass,
                       ssl       = TRUE),
        authenticate = TRUE,
        send = TRUE)
      TRUE
    }, error = function(e) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Sending failed",
        text = paste0("The issue could not be sent (", conditionMessage(e),
                      "). Please contact the maintainer directly."),
        type = "error"
      )
      FALSE
    })
    if (!ok) return(NULL)

    # show confirmation
    shinyWidgets::sendSweetAlert(
      session = session,
      title = "Issue Submitted!",
      text = "Thank you for your feedback! Your issue has been submitted to the maintainer. One of the maintainers will review your issue and contact you for additional details.",
      type = "success"
    )
    
    # reset key inputs
    updateTextInput(
      session,
      "issue_title",
      value = ""
    )
    
    updateTextAreaInput(
      session,
      "issue_description",
      value = ""
    )
    
    shinyWidgets::updatePickerInput(
      session,
      "issue_labels",
      selected = character(0)
    )
  })

 
}
    
## To be copied in the UI
# mod_feedback_ui("feedback_ui_1")
    
## To be copied in the server
# callModule(mod_feedback_server, "feedback_ui_1")
 
