library("CompHeatmaps")
library("shiny")
ui <- pageWithSidebar(
  headerPanel('Generate Heatmap with CompHeatmaps'),
  sidebarPanel(
    tags$h3("Description"),
    uiOutput("description"),
    tags$h3("Instructions"),
    uiOutput("instruction"),

    #uiOutput("instruction")
    fileInput(
      inputId = "input_dir",
      label = "Input Directory",
      placeholder="sample_data_directory"
      ),
    fileInput(
      inputId = "output_dir",
      label = "Output Directory",
      placeholder="filtered_reads_output_directory"
      ),
    actionButton(
      inputId = "run_button",
      label = "Run"
    )
  ),
  mainPanel(
    plotOutput("heatmap")
  )
)


server <- function(input, output, session) {
  output$description <- renderUI(
    HTML("<p> Welcome to the Shiny App for the CompHeatmaps R package
          (Wong, 2023). <code>CompHeatmaps</code> is an R package capable of
          visualizing abundance of raw 16S sequence reads of microbial
          communities and comparing between samples. This allows researchers
          conducting metagenomic analysis microbial communities to visually
          analyze similarities or differences of the highly conserved 16S region
          of bacteria across multiple samples. </p>"
    )
  )
  output$instruction <- renderUI(
    HTML("<p> Below, link a directory containing 16S sequence read data
         of various samples you wish to analyze. For the second input, include
         a valid output directory where the preprocessed data will go to.
         If you wish to run with the example data, you may press the
         <code>Run</code> button without linking any directories.
         ")
  )
  input_dir <- reactive({
    if (is.null(input$input_dir)) {
      return(system.file("extdata/sample_raw_16S_data", package = "CompHeatmaps"))
    }
    else {
      req(input$input_dir)
    }
  })
  output_dir <- reactive({
    if (is.null(input$output_dir)) {
      return(system.file("extdata/filtered_reads", package = "CompHeatmaps"))
    }
  })
  observeEvent(input$run_button, {
    tryCatch({
      path <- system.file("data/precomputed/result", package = "CompHeatmaps")
      result <- readRDS(path)
      #result <- CompHeatmaps::preprocess_16s_data(input_dir(), output_dir(), verbose = TRUE)
      table <- CompHeatmaps::create_abundance_table(result)
      output$heatmap <- CompHeatmaps::create_heatmap(table, proportion = 0.025)
    })
  })
}

shiny::shinyApp(ui = ui, server = server)

