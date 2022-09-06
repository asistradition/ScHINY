## These are the functions that allow the shiny UI to interface with the figure scripts ##
## Each figure script should implement an experiment panel function and a plot function ##
## This script implements a condition panel function and functions to select figures    ##

require(shiny)

# Load the figure scripts
source(file.path(FIGURE.PATH, "figure_1_schematic.R"))
source(file.path(FIGURE.PATH, "figure_2_umap.R"))

# Load the validator functions
source("validators.R")

# A null plotting script for temporary use during startup
# All plot scripts MUST implement the same function signature
plot.null <- function(shiny.data, gene, input) {return(NULL)}
panel.null <- function(...) {return(list())}
describe.null <- function(...) {return("Single Cell Yeast Expression Data")}

cond.panel.standard <- function(selected.conditions=NULL, selected.genotypes=NULL) {
  if(is.null(selected.conditions)) {selected.conditions <- CONDITION.SELECT.DEFAULT}
  if(is.null(selected.genotypes)) {selected.genotypes <- GENOTYPE.LABELS}
  list(
    checkboxGroupInput(inputId = 'conditions', label = 'Conditions', choices = rev(CONDITION.LABELS), selected = selected.conditions),
    checkboxGroupInput(inputId = 'genotypes', label = 'Genotypes', choices = rev(GENOTYPE.LABELS), selected = selected.genotypes)
    )
}

# List of figure plot functions with display name keys
FIGURE.PLOTTER.LIST <<- list("Select Figure" = plot.null, 
                             "Figure 1 - Schematic" = plot.figure.1.schematic,
                             "Figure 2 - UMAP" = plot.figure.2.umap)

# List of figure experiment panel functions with display name keys
FIGURE.EXP.PANEL.LIST <<- list("Select Figure" = panel.null,
                               "Figure 1 - Schematic" = exp.panel.figure.1.schematic,
                               "Figure 2 - UMAP" = exp.panel.figure.2.umap)

# List of figure condition panel functions with display name keys
FIGURE.COND.PANEL.LIST <<- list("Select Figure" = panel.null,
                                "Figure 1 - Schematic" = panel.null,
                                "Figure 2 - UMAP" = panel.null)

# List of validator scripts with display name keys
FIGURE.VALIDATOR.LIST <<- list("Select Figure" = validate.null,
                               "Figure 1 - Schematic" = validate.null,
                               "Figure 2 - UMAP" = process.gene.input)

# List of figure experiment panel functions with display name keys
FIGURE.DESCRIBE.PANEL.LIST <<- list("Select Figure" = describe.null,
                                    "Figure 1 - Schematic" = describe.figure.1.schematic,
                                    "Figure 2 - UMAP" = describe.figure.2.umap)

# Turn a display name into a function to generate a plot
get.data.plotter <<- function(display.name) {
  if (is.null(display.name)) {return(plot.null)}
  else if (!display.name %in% names(FIGURE.PLOTTER.LIST)) {return(plot.null)}
  else {return(FIGURE.PLOTTER.LIST[[display.name]])}
}

# Turn a display name into a function to generate the experiment UI
get.experiment.panel <<- function(display.name) {
  if (is.null(display.name)) {return(panel.null)}
  else if (!display.name %in% names(FIGURE.EXP.PANEL.LIST)) {return(panel.null)}
  else {return(FIGURE.EXP.PANEL.LIST[[display.name]])}
}

# Turn a display name into a function to generate the condition UI
get.condition.panel <<- function(display.name) {
  if (is.null(display.name)) {return(panel.null)}
  else if (!display.name %in% names(FIGURE.COND.PANEL.LIST)) {return(panel.null)}
  else {return(FIGURE.COND.PANEL.LIST[[display.name]])}
}

get.description.text <<- function(display.name) {
  if (is.null(display.name)) {return(describe.null)}
  else if (!display.name %in% names(FIGURE.COND.PANEL.LIST)) {return(describe.null)}
  else {return(FIGURE.DESCRIBE.PANEL.LIST[[display.name]])}
}

# Turn a display name into a function to validate input
get.validator <<- function(display.name) {
  if (is.null(display.name)) {return(validate.always.true)}
  else if (!display.name %in% names(FIGURE.EXP.PANEL.LIST)) {return(validate.always.true)}
  else {return(FIGURE.VALIDATOR.LIST[[display.name]])}
}
