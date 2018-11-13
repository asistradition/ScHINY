SHINY.TITLE <<- "Jackson et al Single Cell Yeast Data"

# Gene strings
GENE.LABEL <<- 'Gene Name'
GENE.DEFAULT <<- 'YEL009C'

# Data Set Loading Information
META.DATA <<- read.table('metadata.csv', header=TRUE, stringsAsFactors = FALSE, sep=",")
DATA.DEFAULT <<- META.DATA[META.DATA$Default,]$Display
DATA.LABEL <<- 'Data Visualization'

# Define the user interface                   
ui <- fluidPage(
  fluidRow(
    column(3,
            selectInput(inputId = 'dataset',
                          label = DATA.LABEL,
                          choices = META.DATA$Display,
                          selected = DATA.DEFAULT),
            textInput(inputId = 'gene',
                      label = GENE.LABEL,
                      value = GENE.DEFAULT),
            uiOutput("ExpPanel")),
      column(9, plotOutput(outputId = 'plots', height = '800'))
  ),
  title = SHINY.TITLE
)
