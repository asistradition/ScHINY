## These are the functions that allow the shiny UI to interface with the figure scripts ##
## Each figure script should implement an experiment panel function and a plot function ##
## This script implements a condition panel function and functions to select figures    ##

require(shiny)

# Load the figure scripts
source(file.path(FIGURE.PATH, "figure_1_schematic.R"))
source(file.path(FIGURE.PATH, "figure_2_ridgeline.R"))
source(file.path(FIGURE.PATH, "figure_2_umap.R"))
source(file.path(FIGURE.PATH, 'figure_3_condition_umap.R'))
source(file.path(FIGURE.PATH, "figure_4_violin.R"))

# Load the validator functions
source("validators.R")

# A null plotting script for temporary use during startup
# All plot scripts MUST implement the same function signature
plot.null <- function(shiny.data, gene, input) {return(NULL)}
panel.null <- function(...) {return(list())}

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
                             "Figure 2 - Ridgeline" = plot.figure.2.ridgeline,
                             "Figure 2 - UMAP" = plot.figure.2.umap,
                             "Figure 3 - Condition UMAP" = plot.figure.3.umap,
                             "Figure 4 - Violin" = plot.figure.4.violin)

# List of figure experiment panel functions with display name keys
FIGURE.EXP.PANEL.LIST <<- list("Select Figure" = panel.null,
                               "Figure 1 - Schematic" = exp.panel.figure.1.schematic,
                               "Figure 2 - Ridgeline" = exp.panel.figure.2.ridgeline,
                               "Figure 2 - UMAP" = exp.panel.figure.2.umap,
                               "Figure 3 - Condition UMAP" = exp.panel.figure.3.umap,
                               "Figure 4 - Violin" = exp.panel.figure.4.violin)

# List of figure condition panel functions with display name keys
FIGURE.COND.PANEL.LIST <<- list("Select Figure" = panel.null,
                                "Figure 1 - Schematic" = panel.null,
                                "Figure 2 - Ridgeline" = panel.null,
                                "Figure 2 - UMAP" = cond.panel.figure.2.umap,
                                "Figure 3 - Condition UMAP" = cond.panel.figure.3.umap,
                                "Figure 4 - Violin" = cond.panel.figure.4.violin)

# List of validator scripts with display name keys
FIGURE.VALIDATOR.LIST <<- list("Select Figure" = validate.null,
                               "Figure 1 - Schematic" = validate.null,
                               "Figure 2 - Ridgeline" = validate.multigene.input,
                               "Figure 2 - UMAP" = process.gene.input,
                               "Figure 3 - Condition UMAP" = process.gene.input,
                               "Figure 4 - Violin" = validate.multigene.input)

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

# Turn a display name into a function to validate input
get.validator <<- function(display.name) {
  if (is.null(display.name)) {return(validate.always.true)}
  else if (!display.name %in% names(FIGURE.EXP.PANEL.LIST)) {return(validate.always.true)}
  else {return(FIGURE.VALIDATOR.LIST[[display.name]])}
}

# Take the result of a validation function and load the gene data that it referrs to
get.gene.data <- function(data.list, genes.to.load) {
  if(is.null(genes.to.load)) {return(data.list)}
  for (gene.name in genes.to.load) {
    if(!gene.name %in% names(data.list)) {
      data.list[[gene.name]] <- read.table(gzfile(file.path(DATA.PATH, paste0(gsub("-", "\\.", gene.name), ".tsv.gz"))), header=TRUE)
    }
  }
  if (length(data.list) > MAX.CONCURRENT.LOADED) {data.list <- data.list[c("meta.data", genes.to.load)]}
  return(data.list)
}

# Take the loaded data and convert it to a dataframe suitable for plotting
process.data.list <- function(data.list, genes.to.load, data.type="counts", scale.data = FALSE) {
  processed.data.frame <- data.list[["meta.data"]]
  if (is.null(genes.to.load)) (return(processed.data.frame))
  for (gene.name in genes.to.load) {
    if (scale.data) {new.data <- scale(data.list[[gene.name]][data.type])}
    else {new.data <- data.list[[gene.name]][data.type]}
    colnames(new.data) <- gene.name
    processed.data.frame <- cbind(processed.data.frame, new.data)
  }
  return(processed.data.frame)
}

