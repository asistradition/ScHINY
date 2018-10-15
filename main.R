library(shiny)

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

# Load the default data set
DEFAULT.FILE.PATH <<- file.path(DATA.PATH, META.DATA[META.DATA$Display == DATA.DEFAULT , 'File'])
if(META.DATA[META.DATA$Default,]$Gzipped) {
  DEFAULT.FILE.PATH <<- gzfile(DEFAULT.FILE.PATH)
}
DEFAULT.SHINY.DATA <<- read.table(DEFAULT.FILE.PATH)
DEFAULT.GENE.MAP <<- read.table(file.path(DATA.PATH, META.DATA[META.DATA$Display == DATA.DEFAULT, 'GeneFile']), col.names = c('Systemic', 'Common'), stringsAsFactors = FALSE)
DEFAULT.ALLOWED.NAMES <<- as.vector(t(DEFAULT.GENE.MAP))

source('ui.R')
source('server.R')
shinyApp(ui=ui, server=server)