library(shiny)

ui <- fluidPage(
  titlePanel(title = div(img(src = "logo.jpg", height = 187, width = 100),
                         "TF DEACoN")
             ),
  br(),
  
  tabsetPanel(type = "tabs",
              tabPanel("Home",
                       br(),
                       sidebarLayout(
                         sidebarPanel(h3("Input Genes"),
                                      textAreaInput("query", label = "Gene IDs (AT IDs)",
                                                    width = "100%", rows = 5, resize = "vertical",
                                                    placeholder = "Separators (commas, etc.) are not important. All valid AT IDs will be recognized, and other text will be ignored."),
                                      actionButton("submit", "Submit")
                         ),
                         mainPanel(
                           tabsetPanel(type = "tabs",
                                       tabPanel("Results",
                                                DT::dataTableOutput("resultsdisplay"),
                                                uiOutput("downloadrender")),
                                       tabPanel("Query",
                                                h3("TF DEACoN recognized the following query"),
                                                htmlOutput("confirmquery")
                                       )
                           ))
                       )),
              
              tabPanel("Help",
                       includeMarkdown("help.Rmd"))
  )
)

  
  # Read data needed:
  # alltarg = all TF to target relationships, pairwise
  alltarg <- read.csv("input/alltarg.txt",
                      header = TRUE,
                      row.names = 1)
  
  # targcounts = genome target counts for all TFs
  #   becomes basis of results dataframe
  targcounts <- read.delim("input/targcounts.txt",
                           header = TRUE,
                           row.names = 1)
  
  server <- function(input, output) {
    # Read input and find all valid gene names
    # Find location of all valid gene names
    locs <- eventReactive(input$submit, {
      gregexpr("AT[1-5]G[0-9][0-9][0-9][0-9]0", input$query, ignore.case = TRUE)
    })
    
    # Create a list of all valid gene names using loc
    # Include validation that there is at least 1 valid AT ID
    querylist <- eventReactive(input$submit, {
      validate(
        need(as.vector(locs()[[1]]) != -1, "Please enter a valid AT ID")
      )
      toupper(regmatches(input$query, locs())[[1]])
    })
    
    # Output the query list
    # This allows the user to double-check that all their gene IDs were recognized
    output$confirmquery <- renderUI({
      str1 <- paste(sort(querylist()), sep = " ", collapse = " ")
      str2 <- paste("Number of genes:", querylen())
      
      HTML(paste(str2, str1, sep = "<br/>"))
    })
    
    # Perform calculations
    # How many genes are in the query?
    querylen <- reactive({length(querylist())})
    
    # Create a smaller version of alltarg which only contains the relevant entries (with query targets)
    querytarg <- reactive({alltarg[which(alltarg$target %in% querylist()),]})
    
    # Create master "results" data frame with calculated ratio, logFC, and p-value  
    results <- reactive({
      querytarg <- querytarg()
      
      # Count the number of query targets each TF has
      results <- targcounts
      results$query.count <- NA
      for (i in 1:nrow(results))
      {results$query.count[i] <- nrow(querytarg[which(querytarg$TF == results$tf.id[i]),])}
      
      # Calculate the query ratio and logFC
      results$query.ratio <- results$query.count / querylen()
      
      results$logFC <- log2(results$query.ratio / results$genome.ratio)
      
      # Calculate the adjusted p-value
      results$adj.pval <- NA
      
      for(i in 1:nrow(results))
      {
        results$adj.pval[i] <- as.numeric(binom.test(results$query.count[i], querylen(), results$genome.ratio[i],
                                                     alternative = "two.sided")$p.value)
      }
      
      results$adj.pval <- p.adjust(results$adj.pval, method = "BH")
      
      # Arrange results from smallest to largest p-value
      results <- results[order(results$adj.pval, results$logFC),]
    })
    
    # Create version of results for displaying on screen
    output$resultsdisplay <- DT::renderDataTable({
      resultsdisplay <- results()
      
      # Round numeric results for display
      resultsdisplay[,c(5,7:9)] <- round(resultsdisplay[,c(5,7:9)], digits = 2)
      
      # Add better-looking column names
      DT::datatable(resultsdisplay, rownames = FALSE,
                    colnames = c("TF ID", "TF family", "Gene name(s)", "Genome target count", "Genome ratio", "Query target count", "Query ratio", "logFC", "Adj p-value")
      )
    })
    
    # Create the download file 
    output$downloadinfo <- downloadHandler(
      filename = function() {
        paste("TFDEACoN_Output_", Sys.time(), ".csv", sep = "")
      },
      
      content = function(file) {
        write.csv(results(), file, row.names = FALSE)
      })
    
    # Render the download button  
    output$downloadrender <- renderUI({
      req(input$submit)
      downloadButton("downloadinfo")
    })
  }
  
  shinyApp(ui = ui, server = server)