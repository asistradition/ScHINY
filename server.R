library(ggplot2)
library(dplyr)

# Load the gene names that are in the data set
GENE.MAP <<- read.table('20181001_genes.tsv', col.names = c('Systemic', 'Common'), stringsAsFactors = FALSE)
ALLOWED.NAMES <<- as.vector(t(GENE.MAP))

# Condition metadata
CONDITIONS <<- read.table('conditions.csv', header=TRUE, stringsAsFactors=FALSE, sep=",")
CONDITION.LABELS <<- structure(as.character(CONDITIONS$PlotName), names=CONDITIONS$DataColumn)
RELABEL.FACETS <<- function(string.label) {return(CONDITION.LABELS[string.label])}

# Load the figure plotting functions
source('figure1_expression_summary.R')
source('figure2_umap.R')
source('figure4_expression_ridgeline.R')

plot_null <- function(meta_data, gene, input, validator = NULL) {}
sidebar_null <- function() {list()}

get_sidebar <- function(data_set) {
  if (is.null(data_set)) {sidebar_null}
  else if (data_set == "Expression Summary") {sidebar_expression_summary}
  else if (data_set == "Figure 2 - UMAP") {sidebar_umap}
  else if (data_set == "Figure 4 - Expression") {sidebar_expression_ridge}
  else {sidebar_null}
}

get_plotter <- function(data_set) {
  if (is.null(data_set)) {plot_null}
  else if (data_set == "Expression Summary") {plot_expression_summary}
  else if (data_set == "Figure 2 - UMAP") {plot_umap}
  else if (data_set == "Figure 4 - Expression") {plot_expression_ridge}
  else {plot_null}
}

get_gene_data <<- function(select.gene, select.path) {
  read.table(file.path(select.path, paste0(gsub("-", "\\.", select.gene), ".tsv")), header=FALSE)
}

default.data <- read.table(file.path(META.DATA[META.DATA$Default,]$Path, META.DATA[META.DATA$Default,]$MetaData))

server <- function(input, output) {
  
  shiny.data <- default.data
  active.data <- reactiveVal(DATA.DEFAULT)
  
  output$ExpPanel <- renderUI({get_sidebar(active.data())()})
  
  # Validate the gene input
  # Convert it to the shiny.data keys if it's in the gene map
  # Otherwise generate a matching genes error message
  validate.gene <- reactive({
    
    # Get the input string and figure out which genes it could be matching
    select.gene <- toupper(input$gene)
    could.be <- ALLOWED.NAMES[startsWith(ALLOWED.NAMES, select.gene)]
    
    # If the list of matching genes is too long, create an error message with the number of genes
    # Otherwise, number of genes and gene names
    if(length(could.be) > 50){
      could.be <- c(length(could.be), "Matching Genes")
    }
    else{
      could.be <- c(length(could.be), "Matching Genes:", could.be)
    }
    
    # If the input gene name an be mapped to the data, proceed
    # Otherwise print the error message generated above
    validate(
      need(select.gene %in% ALLOWED.NAMES, paste(could.be, collapse=" "))
    )
    
    # Make sure that there's one condition selected
    validate(
      need(!is.null(input$conditions), "At least one condition must be selected")
    )
    
    # Convert the input to the same (systemic) format as the Gene column in the data
    if(toupper(input$gene) %in% GENE.MAP$Common) {
      select.gene <- toString(GENE.MAP[GENE.MAP$Common == select.gene,]['Systemic'])
    }
    
    # Return the converted input name
    select.gene
  })
  
  validate.data <- reactive({
    if(input$dataset != active.data()) {
      shiny.data <<- read.table(file.path(META.DATA[META.DATA$Display == input$dataset,]$Path, META.DATA[META.DATA$Display == input$dataset,]$MetaData))
      active.data(input$dataset)
    }
    validate(
      need(active.data() == input$dataset, "Loading new data failed")
    )
  })
  
  # Render the plot
  output$plots <- renderPlot({
    validate.data()
    get_plotter(input$dataset)(shiny.data, input$gene, input, validator = validate.gene)
  })
}