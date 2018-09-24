# Define the user interface                   
ui <- fluidPage(
  titlePanel(SHINY.TITLE),
  sidebarLayout(
    sidebarPanel(
      textInput(inputId = 'gene',
                label = GENE.LABEL,
                value = GENE.DEFAULT),
      uiOutput('conditions'),
      selectInput(inputId = 'yaxis',
                  label = YAXIS.SCALING.LABEL,
                  choices = c(YAXIS.UMI.COUNT, YAXIS.LIB.NORM, YAXIS.WT.NORM),
                  selected = YAXIS.SCALING.DEFAULT),
      selectInput(inputId = 'yaxislim',
                  label = YAXIS.LIMIT.LABEL,
                  choices = c(YAXIS.FIXED, YAXIS.FREE),
                  selected = YAXIS.LIMIT.DEFAULT),
      radioButtons(inputId = 'wtline',
                   label = YAXIS.HLINE.LABEL,
                   choices = c(YAXIS.SHOW, YAXIS.HIDE),
                   inline = TRUE,
                   selected = YAXIS.HLINE.DEFAULT),
      selectInput(inputId = 'dataset',
                  label = DATA.LABEL,
                  choices = META.DATA$Display,
                  selected = DATA.DEFAULT)
    ),
    mainPanel(
      plotOutput(outputId = 'plots', height = '800')
    )
  )
)
