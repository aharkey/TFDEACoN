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
    mainPanel(h3("Results"),
              uiOutput("downloadrender"),
              DT::dataTableOutput("results"),
              htmlOutput("summary")
              )
    )
  )


server <- function(input, output) {
  alltarg <- read.csv("input/alltarg.txt",
                      header = TRUE,
                      row.names = 1)
  
  targcounts <- read.delim("input/targcounts.txt",
                         header = TRUE,
                         row.names = 1)
  
  colnames(targcounts) <- c("TF ID", "Family", "Gene name(s)", "All targets", "Background ratio")
  
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
    
    HTML(paste(paste("<h3>Input Genes</h3>", str1, sep = ""), str2, sep = "<br/>"))
  })
  
  results <- reactive({
    targcounts$"Query targets" <- NA
    
    subtarg <- subtarg()
    
    for (i in 1:nrow(targcounts))
    {
      targcounts$"Query targets"[i] <- nrow(subtarg[which(subtarg$TF == targcounts$"TF ID"[i]),])
    }
    
    targcounts$"Query ratio" <- targcounts$"Query targets" / targnum()
    
    targcounts$logFC <- log2(targcounts$"Query ratio" / targcounts$"Background ratio")
    
    targcounts$"Adj P value" <- NA
    
    for(i in 1:nrow(targcounts))
    {
      targcounts$"Adj P value"[i] <- as.numeric(binom.test(targcounts$"Query targets"[i], targnum(), targcounts$"Background ratio"[i],
                                                           alternative = "two.sided")$p.value)
    }
    
    results <- targcounts
    
    results <- results[order(results$"Adj P value", results$logFC),]
  })
  
  resultsdisplay <- reactive({
    resultsdisplay <- results()
    
    resultsdisplay[,c(5,7:9)] <- round(resultsdisplay[,c(5,7:9)], digits = 2)
    
    resultsdisplay <- resultsdisplay
  })
  
  output$results <- DT::renderDataTable({  
    DT::datatable(resultsdisplay(), rownames = FALSE)
  })
  
  output$downloadinfo <- downloadHandler(
    filename = function() {
      paste("TFDEACoN_Output_", Sys.time(), ".csv", sep = "")
    },
    
    content = function(file) {
      write.csv(results(), file, row.names = FALSE)
    })
  
  output$downloadrender <- renderUI({
    req(input$submit)
    downloadButton("downloadinfo")
  })
}




shinyApp(ui = ui, server = server)