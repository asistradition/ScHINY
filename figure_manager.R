## These are the functions that allow the shiny UI to interface with the figure scripts ##
## Each figure script should implement an experiment panel function and a plot function ##
## This script implements a condition panel function and functions to select figures    ##

require(shiny)

# Load the figure scripts
source(file.path(FIGURE.PATH, "network_figure.R"))
source(file.path(FIGURE.PATH, "expression_figure.R"))

# Load the validator functions
source("validators.R")

# A null plotting script for temporary use during startup
# All plot scripts MUST implement the same function signature
plot.null <- function(shiny.data, gene, input) {return(NULL)}
panel.null <- function(...) {return(list())}
describe.null <- function(...) {return("Rice Network Data")}

validate.gene.input.osi <- function(x, cn) {validate.gene.input(x, "OSI", use.common.names=cn)}
validate.gene.input.osj <- function(x, cn) {validate.gene.input(x, "OSJ", use.common.names=cn)}
validate.gene.input.og <- function(x, cn) {validate.gene.input(x, "OG", use.common.names=cn)}

validate.multigene.input.osi <- function(x, cn) {validate.multigene.input(x, "OSI", use.common.names=cn)}
validate.multigene.input.osj <- function(x, cn) {validate.multigene.input(x, "OSJ", use.common.names=cn)}
validate.multigene.input.og <- function(x, cn) {validate.multigene.input(x, "OG", use.common.names=cn)}

# List of figure plot functions with display name keys
FIGURE.PLOTTER.LIST <<- list(
  "Select Figure" = plot.null,
  "Oryza sativa indica Gene Expression" = function(...) {plot.expression.over.time("OSI", ...)},
  "Oryza sativa japonica Gene Expression" = function(...) {plot.expression.over.time("OSJ", ...)},
  "Oryza glaberrima Gene Expression" =  function(...) {plot.expression.over.time("OG", ...)},
  "Oryza sativa indica - (BBSR/EGRIN)" = function(...) {plot.figure.network(NETWORK.DATA[["OSI_BBSR_EGRIN"]], "OSI", ...)},
  "Oryza sativa japonica - (BBSR/EGRIN)" = function(...) {plot.figure.network(NETWORK.DATA[["OSJ_BBSR_EGRIN"]], "OSJ", ...)},
  "Oryza glaberrima - (BBSR/EGRIN)" =  function(...) {plot.figure.network(NETWORK.DATA[["OG_BBSR_EGRIN"]], "OG", ...)},
  "Oryza sativa indica - (BBSR/PlantRegDB)" = function(...) {plot.figure.network(NETWORK.DATA[["OSI_BBSR_PREGDB"]], "OSI", ...)},
  "Oryza sativa japonica - (BBSR/PlantRegDB)" = function(...) {plot.figure.network(NETWORK.DATA[["OSJ_BBSR_PREGDB"]], "OSJ", ...)},
  "Oryza glaberrima - (BBSR/PlantRegDB)" = function(...) {plot.figure.network(NETWORK.DATA[["OG_BBSR_PREGDB"]], "OG", ...)}
)

# List of validator scripts with display name keys
FIGURE.VALIDATOR.LIST <<- list(
  "Select Figure" = validate.null,
  "Oryza sativa indica Gene Expression" = validate.gene.input.osi,
  "Oryza sativa japonica Gene Expression" = validate.gene.input.osj,
  "Oryza glaberrima Gene Expression" =  validate.gene.input.og,
  "Oryza sativa indica - (BBSR/EGRIN)" = validate.multigene.input.osi,
  "Oryza sativa japonica - (BBSR/EGRIN)" = validate.multigene.input.osj,
  "Oryza glaberrima - (BBSR/EGRIN)" = validate.multigene.input.og,
  "Oryza sativa indica - (BBSR/PlantRegDB)" = validate.multigene.input.osi,
  "Oryza sativa japonica - (BBSR/PlantRegDB)" = validate.multigene.input.osj,
  "Oryza glaberrima - (BBSR/PlantRegDB)" = validate.multigene.input.og
)
                               

# List of figure experiment panel functions with display name keys
FIGURE.DESCRIBE.PANEL.LIST <<- list(
  "Select Figure" = describe.null,
  "Oryza sativa indica Gene Expression" = function(...) {describe.expression.over.time("OSI", ...)},
  "Oryza sativa japonica Gene Expression" = function(...) {describe.expression.over.time("OSJ", ...)},
  "Oryza glaberrima Gene Expression" =  function(...) {describe.expression.over.time("OG", ...)},
  "Oryza sativa indica - (BBSR/EGRIN)" = function(...) {describe.figure.network(NETWORK.DATA[["OSI_BBSR_EGRIN"]], "OSI", ...)},
  "Oryza sativa japonica - (BBSR/EGRIN)" = function(...) {describe.figure.network(NETWORK.DATA[["OSJ_BBSR_EGRIN"]], "OSJ", ...)},
  "Oryza glaberrima - (BBSR/EGRIN)" =  function(...) {describe.figure.network(NETWORK.DATA[["OG_BBSR_EGRIN"]], "OG", ...)},
  "Oryza sativa indica - (BBSR/PlantRegDB)" = function(...) {describe.figure.network(NETWORK.DATA[["OSI_BBSR_PREGDB"]], "OSI", ...)},
  "Oryza sativa japonica - (BBSR/PlantRegDB)" = function(...) {describe.figure.network(NETWORK.DATA[["OSJ_BBSR_PREGDB"]], "OSJ", ...)},
  "Oryza glaberrima - (BBSR/PlantRegDB)" = function(...) {describe.figure.network(NETWORK.DATA[["OG_BBSR_PREGDB"]], "OG", ...)}
)

get.data.plotter <<- function(display.name) {
  if (is.null(display.name)) {return(plot.null)}
  else if (!display.name %in% names(FIGURE.PLOTTER.LIST)) {return(plot.null)}
  else {return(FIGURE.PLOTTER.LIST[[display.name]])}
}

get.description.text <<- function(display.name) {
  if (is.null(display.name)) {return(describe.null)}
  else if (!display.name %in% names(FIGURE.PLOTTER.LIST)) {return(describe.null)}
  else {return(FIGURE.DESCRIBE.PANEL.LIST[[display.name]])}
}

# Turn a display name into a function to validate input
get.validator <<- function(display.name) {
  if (is.null(display.name)) {return(validate.always.true)}
  else if (!display.name %in% names(FIGURE.PLOTTER.LIST)) {return(validate.always.true)}
  else {return(FIGURE.VALIDATOR.LIST[[display.name]])}
}
