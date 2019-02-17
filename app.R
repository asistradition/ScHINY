require(shiny)
source("settings.R")

# Define the user interface                   
ui <- fluidPage(
  fluidRow(
    column(3,
           wellPanel(
            # This is the selection for the figure aspect
             # It's outside the reactive panels so it doesn't reload all the time
            selectInput(inputId = 'dataset',
                        label = FIGURE.SELECT.LABEL,
                        choices = names(FIGURE.PLOTTER.LIST),
                        selected = DEFAULT.FIGURE),
            
            # This is the gene-specific textbox
            # It's outside the reactive panels so it doesn't reload all the time
            textInput(inputId = 'gene',
                      label = GENE.LABEL,
                      value = GENE.DEFAULT),
            
            # This is the figure aspect-specific panel
            uiOutput("ExpPanel"),
            
            # This is the data selection-specific panel
            # It's separate from the figure panel so it doesn't get reloaded when
            # Switching between figures that would keep the same data selection panel
            uiOutput("CondPanel")
           )
         ),
      column(9, plotOutput(outputId = 'plots', height = '800'))
  ),
  title = SHINY.TITLE
)

# Define the server interface
server <- function(input, output, session) {
  
  # Begin with the default data and default figure
  shiny.data <- list(meta.data = META.DATA)
  active.figure <- shiny::reactiveVal(DEFAULT.FIGURE)
  output$ExpPanel <- get.experiment.panel(DEFAULT.FIGURE)
  output$CondPanel <- get.condition.panel(DEFAULT.FIGURE)
  
  # Validate the gene input
  validate.gene <- shiny::reactive({get.validator(input$dataset)(input$gene)})
  
  # Make sure that the active.figure is correct (this is needed to make the reactive UI elements work)
  update.figure <- shiny::reactive({
    if(input$dataset != active.figure()) {
      active.figure(input$dataset)
      output$ExpPanel <- shiny::renderUI({get.experiment.panel(active.figure())()})
      output$CondPanel <- shiny::renderUI({get.condition.panel(active.figure())(selected.conditions = input$conditions, selected.genotypes = input$genotypes)})
    }
    shiny::validate(shiny::need(active.figure() == input$dataset, "Loading new figure failed"))
  })
  
  # Render the plot in response to UI changes
  output$plots <- shiny::renderPlot({
    update.figure()
    validated.genes <- validate.gene()
    shiny.data <- get.gene.data(shiny.data, validated.genes)
    get.data.plotter(active.figure())(shiny.data, validated.genes, input)
  })
}

shiny::shinyApp(ui = ui, server = server)