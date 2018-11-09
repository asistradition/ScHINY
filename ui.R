SHINY.TITLE <<- "Single-Cell Yeast Gene Expression"

# Gene strings
GENE.LABEL <<- 'Gene Name'
GENE.DEFAULT <<- 'YEL009C'

# Plot visualization strings
# Units
YAXIS.SCALING.LABEL <<- 'Y-Axis Units'
YAXIS.SCALING.DEFAULT <<- 'Per-Cell %'

YAXIS.UMI.COUNT <<- 'UMI Count'
YAXIS.UMI.COUNT.LABEL <<- 'UMI Count Per Cell'

YAXIS.LIB.NORM <<- 'Per-Cell %'
YAXIS.LIB.NORM.LABEL <<- 'Per Cell Relative Transcript Abundance (%)'

YAXIS.WT.NORM <<- 'Relative To WT'
YAXIS.WT.NORM.LABEL <<- 'Per Cell Relative Transcript Abundance\n(Compared to WT)'

# Limits
YAXIS.LIMIT.LABEL <<- 'Y-Axis Scaling'
YAXIS.LIMIT.DEFAULT <<- 'Fixed'
YAXIS.FIXED <<- 'Fixed'
YAXIS.FREE <<- 'Free'

# WT Horizontal Line
YAXIS.HLINE.LABEL <<- 'WT Average'
YAXIS.HLINE.DEFAULT <<- 'Show'
YAXIS.SHOW <<- 'Show'
YAXIS.HIDE <<- 'Hide'

# Data Set Loading Information
META.DATA <<- read.table('metadata.csv', header=TRUE, stringsAsFactors = FALSE, sep=",")
DATA.PATH <<- 'data'
DATA.DEFAULT <<- META.DATA[META.DATA$Default,]$Display
DATA.LABEL <<- 'Data Set'

# Define the user interface                   
ui <- fluidPage(
  titlePanel(SHINY.TITLE),
  sidebarLayout(
    sidebarPanel(
      uiOutput('gene'),
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
      verbatimTextOutput('plots')
     # plotOutput(outputId = 'plots', height = '800')
    )
  )
)
