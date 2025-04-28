library(shiny)
library(ensembldb)
library(AnnotationHub)
library(Biostrings)
library(DT)
# Source our utility functions
source("exon_skipping_utils.R")

# Define UI
ui <- fluidPage(
  titlePanel("ensembldb_test"),
  
  sidebarLayout(
    sidebarPanel(
      # Input for protein selection
      textInput("proteinInput", "Enter Protein ID or Gene Symbol:", "KLHL17"),
      
      # Radio buttons to choose search type
      radioButtons("searchType", "Search By:", 
                   choices = c("Gene Symbol" = "SYMBOL", 
                               "Protein ID" = "PROTEINID"),
                   selected = "SYMBOL"),
      
      # Action button to trigger search
      actionButton("searchBtn", "Search Protein"),
      
      hr(),
      
      # Only show when exons are available
      conditionalPanel(
        condition = "output.exonsAvailable",
        selectInput("exonSelect", "Select Exon to Remove:", choices = NULL),
        actionButton("removeBtn", "Generate Modified Sequence")
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Protein Information", 
                 verbatimTextOutput("proteinInfo")),
        
        tabPanel("Exon Structure",
                 DT::dataTableOutput("exonTable")),
        
        tabPanel("Sequences",
                 h4("Original Protein Sequence"),
                 verbatimTextOutput("originalSeq"),
                 
                 hr(),
                 h4("Modified Protein Sequence (without selected exon)"),
                 verbatimTextOutput("modifiedSeq"),
                 
                 hr(),
                 h4("Analysis Details"),
                 verbatimTextOutput("analysisDetails"))
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Initialize database connection
  ensdb <- reactive({
    withProgress(message = 'Loading Ensembl database...', {
      hub <- AnnotationHub()
      hub[["AH119325"]]  # Human Ensembl database
    })
  })
  
  # Store protein data
  proteinData <- reactiveVal(NULL)
  exonData <- reactiveVal(NULL)
  transcriptData <- reactiveVal(NULL)
  analysisDetails <- reactiveVal("")
  
  # When search button is clicked
  observeEvent(input$searchBtn, {
    req(input$proteinInput)
    
    # Get protein information based on search type
    if (input$searchType == "SYMBOL") {
      # Get all transcripts for the gene symbol
      txs <- transcripts(ensdb(), filter = list(SymbolFilter(input$proteinInput)))
      transcriptData(txs)
      
      # Get proteins for these transcripts
      if (length(txs) > 0) {
        prts <- proteins(ensdb(), filter = list(TxIdFilter(txs$tx_id)))
        proteinData(prts)
        
        # Get exons for these transcripts
        exns <- exonsBy(ensdb(), by = "tx", filter = list(TxIdFilter(txs$tx_id)))
        exonData(exns)
      } else {
        proteinData(NULL)
        exonData(NULL)
        transcriptData(NULL)
        showNotification("Gene symbol not found", type = "error")
      }
    } else {
      # Direct protein ID search
      prts <- proteins(ensdb(), filter = list(ProteinIdFilter(input$proteinInput)))
      
      if (length(prts) > 0) {
        proteinData(prts)
        
        # Get the transcripts for these proteins
        txIds <- prts$tx_id
        txs <- transcripts(ensdb(), filter = list(TxIdFilter(txIds)))
        transcriptData(txs)
        
        # Get exons for these transcripts
        exns <- exonsBy(ensdb(), by = "tx", filter = list(TxIdFilter(txIds)))
        exonData(exns)
      } else {
        proteinData(NULL)
        exonData(NULL)
        transcriptData(NULL)
        showNotification("Protein ID not found", type = "error")
      }
    }
    
    # Update exon dropdown choices if exons are found
    if (!is.null(exonData()) && length(exonData()) > 0) {
      # Get exon IDs and positions for the first transcript
      firstTx <- names(exonData())[1]
      exonInfo <- exonData()[[firstTx]]
      exonChoices <- paste0("Exon ", 1:length(exonInfo), 
                           " (", start(exonInfo), "-", end(exonInfo), ")")
      
      updateSelectInput(session, "exonSelect", choices = setNames(1:length(exonInfo), exonChoices))
    }
  })
  
  # Output protein information
  output$proteinInfo <- renderPrint({
    req(proteinData())
    print(proteinData())
  })
  
  # Display if exons are available
  output$exonsAvailable <- reactive({
    return(!is.null(exonData()) && length(exonData()) > 0)
  })
  outputOptions(output, "exonsAvailable", suspendWhenHidden = FALSE)
  
  # Display exon table
  output$exonTable <- DT::renderDataTable({
    req(exonData())
    if (length(exonData()) > 0) {
      firstTx <- names(exonData())[1]
      exons <- exonData()[[firstTx]]
      
      # Create the exon table
      data.frame(
        Exon_Number = 1:length(exons),
        Exon_ID = mcols(exons)$exon_id,
        Start = start(exons),
        End = end(exons),
        Width = width(exons),
        Rank = mcols(exons)$exon_rank
      )
    }
  })
  
  # Original protein sequence
  output$originalSeq <- renderText({
    req(proteinData())
    if (length(proteinData()) > 0) {
      return(proteinData()$protein_sequence[1])
    } else {
      return("No protein sequence available")
    }
  })
  
  # Generate modified sequence when button clicked
  modifiedSequence <- reactiveVal("")
  
  observeEvent(input$removeBtn, {
    req(proteinData(), exonData(), input$exonSelect, transcriptData())
    
    # Get the transcript ID
    txId <- proteinData()$tx_id[1]
    
    # Get the selected exon
    firstTx <- names(exonData())[1]
    selectedExonIdx <- as.numeric(input$exonSelect)
    selectedExon <- exonData()[[firstTx]][selectedExonIdx]
    
    # Use our utility function to simulate exon skipping
    result <- simulateExonSkippingAtDnaLevel(ensdb(), txId, selectedExon)
    
    # Update outputs
    modifiedSequence(result$modified)
    analysisDetails(result$message)
  })
  
  # Display modified sequence
  output$modifiedSeq <- renderText({
    if (modifiedSequence() != "") {
      return(modifiedSequence())
    } else {
      return("Select an exon and click 'Generate Modified Sequence'")
    }
  })
  
  # Display analysis details
  output$analysisDetails <- renderText({
    if (analysisDetails() != "") {
      return(analysisDetails())
    } else {
      return("No analysis performed yet")
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)