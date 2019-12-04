library(shiny)

ui <- fluidPage(
  titlePanel("Single Gene Search"),
  
  sidebarLayout(
    sidebarPanel(h2("Input"),
                 textInput("gene", label = "Gene ID (AT ID)"),
                 actionButton("submit", "Submit")
                 ),
    mainPanel(h3("Results"),
              textOutput("summary"),
              textOutput("results")
              
    )
  )
)


server <- function(input, output) {
  alltarg <- read.csv("input/alltarg.txt",
                      header = TRUE,
                      row.names = 1)
  
  loc <- eventReactive(input$submit, {
    regexpr("AT[1-5]G[0-9][0-9][0-9][0-9]0", toupper(input$gene))[1]
  })
  
  gene <- eventReactive(input$submit, {
    substr(toupper(input$gene), loc(), loc()+8)
  })
  
  output$summary <- renderText({
    validate(
      need(loc() != -1, "Please enter a valid AT number")
    )
    paste("Input Gene: ", gene())
  })
  
  output$results <- renderText({
    as.character(alltarg$TF[which(alltarg$target == gene())])
  })
}

shinyApp(ui = ui, server = server)