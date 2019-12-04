library(shiny)

ui <- fluidPage(
  titlePanel("TF DEACoN"),
  
  sidebarLayout(
    sidebarPanel(h3("Input Genes Here"),
                 textAreaInput("genes", label = "Gene IDs (AT IDs)",
                               width = "100%", rows = 5, resize = "vertical",
                               placeholder = "Separators (commas, etc.) are not important. All valid AT IDs will be recognized, and other text will be ignored."),
                 actionButton("submit", "Submit")
    ),
    mainPanel(h3("Input Genes"),
              htmlOutput("summary"),
              h3("Results"),
              downloadButton("download"),
              tableOutput("results")
              
    )
  )
)


server <- function(input, output) {
  alltarg <- read.csv("input/alltarg.txt",
                      header = TRUE,
                      row.names = 1)
  
  targcounts <- read.csv("input/targcounts.txt",
                         header = TRUE,
                         row.names = 1)
  
  locs <- eventReactive(input$submit, {
    gregexpr("AT[1-5]G[0-9][0-9][0-9][0-9]0", input$genes, ignore.case = TRUE)
  })
  
  genelist <- eventReactive(input$submit, {
    toupper(regmatches(input$genes, locs())[[1]])
  })
  
  targnum <- reactive({length(genelist())})
  
  subtarg <- reactive({alltarg[which(alltarg$target %in% genelist()),]})
  
  output$summary <- renderUI({
    validate(
      need(as.vector(locs()[[1]]) != -1, "Please enter a valid AT number")
    )
    str1 <- paste(sort(genelist()), sep = " ", collapse = " ")
    str2 <- paste("Number of genes:", targnum())
    
    HTML(paste(str1, str2, sep = "<br/>"))
  })
  
  results <- reactive({
    targcounts$subtargets <- NA
    
    subtarg <- subtarg()
    
    for (i in 1:nrow(targcounts))
    {
      targcounts$subtargets[i] <- nrow(subtarg[which(subtarg$TF == targcounts$TF[i]),])
    }
    
    targcounts$subratio <- targcounts$subtargets / targnum()
    
    targcounts$logFC <- log2(targcounts$subratio / targcounts$genomeratio)
    
    targcounts$pval <- NA
    
    for(i in 1:nrow(targcounts))
    {
      targcounts$pval[i] <- as.numeric(binom.test(targcounts$subtargets[i], targnum(), targcounts$genomeratio[i],
                                                  alternative = "two.sided")$p.value)
    }
    
    #results <- targcounts[which(targcounts$logFC >= 0),]
    results <- targcounts
    
    results <- results[order(results$pval, results$logFC),]
    
  })
  
  output$results <- renderTable({  
    print(results())
  },
  bordered = TRUE,
  width = "90%",
  align = "l",
  digits = 3)
  
  output$download <- downloadHandler(
    filename = function() {
      paste("ReverseSearchOutput-", Sys.time(), ".csv", sep = "")
    },
    
    content = function(file) {
      write.csv(results(), file, row.names = FALSE)
    }
  )
}




shinyApp(ui = ui, server = server)