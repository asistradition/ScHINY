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
