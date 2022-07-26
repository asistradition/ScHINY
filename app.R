require(shiny)
source("settings.R")

# Define the user interface                   
ui <- fluidPage(
  # This is a css hack for plot height I stole from stack overflow
  # https://stackoverflow.com/questions/26782041/scaling-shiny-plots-to-window-height/26785047
  tags$head(tags$style(HTML(".shiny-plot-output{height:100vh !important;}"))),
  
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
            
            checkboxInput(inputId = 'common_names',
                          label = "Common Names",
                          value = TRUE),
            
            shiny::tags$b(textOutput("Description_header")),
            
            textOutput("Description"),
            
            tags$br(),
            tags$div(shiny::a(href="https://github.com/asistradition/ScHINY", "Shiny Code (GitHub)")),
            tags$div(shiny::p(paste0("App Version: ", VERSION.STRING)))
           )
         ),
      column(9, plotOutput(outputId = 'plots'))
  ),
  title = SHINY.TITLE
)

# Define the server interface
server <- function(input, output, session) {
  
  # Begin with the default data and default figure
  shiny.data <- list()
  active.figure <- shiny::reactiveVal(DEFAULT.FIGURE)

  # Validate the gene input
  validate.gene <- shiny::reactive({get.validator(input$dataset)(input$gene, input$common_names)})
  
  # Make sure that the active.figure is correct (this is needed to make the reactive UI elements work)
  update.figure <- shiny::reactive({
    if(input$dataset != active.figure()) {
      active.figure(input$dataset)
      output$Description <- shiny::renderText({get.description.text(active.figure())()})
      output$Description_header <- shiny::renderText({"Description"})
    }
    shiny::validate(shiny::need(active.figure() == input$dataset, "Loading new figure failed"))
  })
  
  # Render the plot in response to UI changes
  output$plots <- shiny::renderPlot({
    update.figure()
    validated.genes <- validate.gene()
    get.data.plotter(active.figure())(NULL, validated.genes, input)
  })
}

shiny::shinyApp(ui = ui, server = server)