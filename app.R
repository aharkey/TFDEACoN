library(shiny)
library(shinyBS)

##### UI #####
ui <- fluidPage(
  theme = shinythemes::shinytheme("cosmo"),
  
  tags$style(HTML("
    .logo { margin-left: 5%; margin-right: 5%; }
  ")),
  
  titlePanel(title = div(img(src = "logo.jpg", height = 187, width = 100, class = "logo"),
                         "TF DEACoN"),
             windowTitle = "TF DEACoN"
  ),
  br(),
  
  navbarPage("",
             tabPanel("Home",
                      br(),
                      sidebarLayout(
                        sidebarPanel(h3("Input Genes"),
                                     textAreaInput("query", label = "Gene IDs (AT IDs)",
                                                   width = "100%", rows = 5, resize = "vertical",
                                                   placeholder = "e.g. AT5G11100, AT5G20640, AT5G20830, AT1G77380, AT5G63850, AT5G59090"),
                                     bsTooltip("query",
                                               "Separators (commas, etc.) are not important. All valid AT IDs will be recognized, and other text will be ignored."),
                                     checkboxInput("pvalfilt",
                                                   "Apply p-value filter"),
                                     bsTooltip("pvalfilt",
                                               "Set a maximum p-value cutoff"),
                                     uiOutput("pvalmenu"),
                                     checkboxInput("logFCfilt",
                                                   "Apply logFC filter"),
                                     bsTooltip("logFCfilt",
                                               "Set a minimum logFC cutoff. logFC is a measure of TF target enrichment See Help for more info."),
                                     uiOutput("logFCmenu"),
                                     actionButton("submit", "Submit"),
                                     br(), br(),
                                     actionButton("reset", "Reset")
                        ),
                        mainPanel(
                          tabsetPanel(type = "tabs",
                                      tabPanel("Results",
                                               htmlOutput("querywarn"),
                                               DT::dataTableOutput("resultsdisplay"),
                                               uiOutput("downloadrender"),
                                               br()),
                                      tabPanel("Query",
                                               htmlOutput("confirmquery")
                                      )
                          ))
                      )),
             
             tabPanel("Help",
                      includeMarkdown("help.Rmd")),
             tabPanel("Citations",
                      includeMarkdown("citations.Rmd")),
             tabPanel("Contact",
                      includeMarkdown("contact.Rmd"))
  )
)

##### Load background files #####
# Read data needed:
# alltarg = all TF to target relationships, pairwise
alltarg <- read.csv("input/alltarg.txt",
                    header = TRUE)

# targcounts = genome target counts for all TFs
#   becomes basis of results dataframe
targcounts <- read.csv("input/targcounts.csv",
                       header = TRUE)

##### Server #####
server <- function(input, output, session) {
  
  # Create logFC and pval filter options
  output$logFCmenu <- renderUI({
    if (input$logFCfilt == 0) {
      return(NULL)
    }
    
    if(input$logFCfilt == 1) {
      sliderInput("logFCvalue",
                  "logFC cutoff:",
                  min = 0,
                  max = 3,
                  value = 0)
    }
  })
  
  output$pvalmenu <- renderUI({
    if (input$pvalfilt == 0) {
      return(NULL)
    }
    
    if(input$pvalfilt == 1) {
      selectInput("pvalvalue",
                  "p-value cutoff:",
                  list(0.05, 0.01, 0.001))
    }
  })
  
  # Calculations to be performed when Submit is pressed
  observeEvent(input$submit, {
    # Find location of all valid gene names
    locs <- reactive({gregexpr("AT[1-5]G[0-9][0-9][0-9][0-9][0-9]", input$query, ignore.case = TRUE)})
    
    # Validate that there is at least 1 valid AT ID
    validate(need(as.vector(locs()[[1]]) != -1, "Please enter a valid AT ID"))
    
    # Create a list of all valid gene names using loc
    querylist <- reactive({toupper(regmatches(input$query, locs())[[1]])})
    
    # Perform calculations
    # How many genes are in the query?
    querylen <- reactive({length(querylist())})
    
    # Output the query list
    # This allows the user to double-check that all their gene IDs were recognized
    output$confirmquery <- renderUI({
      str1 <- "<h3>TF DEACoN recognized the following query</h3>"
      str2 <- paste("Number of genes:", querylen())
      str3 <- paste(sort(querylist()), sep = "</br>", collapse = "</br>")
      
      HTML(paste(str1, str2, str3, sep = "</br>"))
    })
    
    # Output warning if query list is less than 100 genes
    output$querywarn <- renderUI({
      if(querylen() <= 100)
      {
        HTML("</br><p style='color:blue;'>The query input is 100 genes or less. Note that a small sample size increases the likelihood of false negatives. See the help document for more details.</p></br>")
      }
      else
      {
        HTML("")
      }
    })
    
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
      
      # Calculate Fisher's exact test
      G <- 27655
      
      results$pval <- NA
      
      for(i in 1:nrow(results))
      {
        Q <- querylen()
        Ta <- results$genome.count[i]
        Qta <- results$query.count[i]
        
        A <- G - (Ta + Q - Qta)
        B <- Ta - Qta
        C <- Q - Qta
        D <- Qta
        
        results$pval[i] <- as.numeric(
          fisher.test(
            rbind(
              c(A, B),
              c(C, D)
            )
          )$p.value)
      }
      
      results$adj.pval <- p.adjust(results$pval, method = "BH")
      
      if(input$logFCfilt == TRUE){
        results <- results[which(results$logFC >= input$logFCvalue),]
      }
      
      if(input$pvalfilt == TRUE){
        results <- results[which(results$adj.pval <= input$pvalvalue),]
      }
      # Arrange results from smallest to largest p-value
      results <- results[order(results$adj.pval),]
    })
    
    # Create version of results for displaying on screen
    # With rounded values and display column names
    output$resultsdisplay <- DT::renderDataTable({
      resultsdisplay <- results()
      
      # Round numeric results for display
      resultsdisplay[,c(5,7:10)] <- round(resultsdisplay[,c(5,7:10)], digits = 2)
      
      # Add better-looking column names
      DT::datatable(resultsdisplay, rownames = FALSE,
                    colnames = c("TF ID", "TF family", "Gene name(s)", "Genome target count", "Genome ratio", "Query target count", "Query ratio", "logFC", "p-value", "Adj p-val")
      )
    })
    
    # Create the download file 
    output$downloadinfo <- downloadHandler(
      filename = function() {
        paste0("TFDEACoN_Output_", Sys.time(), ".csv")
      },
      
      content = function(file) {
        write.csv(results(), file, row.names = FALSE)
      })
    
    # Render the download button
    output$downloadrender <- renderUI({
      downloadButton("downloadinfo")
    })
  })
  
  # Settings to be reset when Reset genes is pressed
  observeEvent(input$reset, {
    updateTextAreaInput(session, "query", value = "")
    
    updateCheckboxInput(session, "pvalfilt", value = FALSE)
    
    updateCheckboxInput(session, "logFCfilt", value = FALSE)
    
    output$confirmquery <- renderUI({HTML("")})
    
    output$querywarn <- renderUI({HTML("")})
    
    output$resultsdisplay <- DT::renderDataTable({
      DT::datatable(data.frame(matrix("", nrow = 0, ncol = 10)), rownames = FALSE,
                    colnames = c("TF ID", "TF family", "Gene name(s)", "Genome target count", "Genome ratio", "Query target count", "Query ratio", "logFC", "p-value", "Adj p-val")
      )
    })
    
    output$downloadrender <- renderUI({HTML("")})
  })
}

shinyApp(ui = ui, server = server)