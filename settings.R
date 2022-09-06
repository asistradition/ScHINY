VERSION.STRING <- "September 5, 2022 (RAPA v1)"

## This file should be sourced at the top of all the shiny app scripts ##
## It loads necessary packages and populates the global namespace ##

# Load shiny
library(shiny)

# Load plotting stuff
library(ggplot2)
library(ggridges)
library(RColorBrewer)
library(cowplot)
library(scales)
library(viridis)
library(network)
library(sna)
library(GGally)

# Load tidyverse stuff
library(dplyr)
library(readr)

## Function to load data from tsv.gz files ##
load.tsv.gz <- function(file.name, index.name='gene') {
  x <- readr::read_tsv(file.name)
  if('...1' %in% colnames(x)) {
    x[index.name] <- x['...1']
    x['...1'] <- NULL
  }
  return(x)
}

# Path and file names
if(!exists('DATA.PATH')) {DATA.PATH <<- 'data'}
if(!exists('FIGURE.PATH')) {FIGURE.PATH <<- 'figures'}
if(!exists('IMAGE.PATH')) {IMAGE.PATH <<- 'images'}

if(!exists('GENE.META.DATA.FILE')) {GENE.META.DATA.FILE <<- '2021_RAPA_GENE_METADATA.tsv.gz'}
if(!exists('META.DATA.FILE')) {META.DATA.FILE <<- '2021_RAPA_METADATA.tsv.gz'}

# Set labels and defaults for the UI
if(!exists('SHINY.TITLE')) SHINY.TITLE <<- "Continuous Single-Cell Response to Rapamycin"
if(!exists('GENE.LABEL')) GENE.LABEL <<- 'Gene Name'
if(!exists('GENE.DEFAULT')) GENE.DEFAULT <<- 'YEL009C'
if(!exists('FIGURE.SELECT.LABEL')) {FIGURE.SELECT.LABEL <<- "Figure"}
if(!exists('DEFAULT.FIGURE')) {DEFAULT.FIGURE <<- "Select Figure"}

# Set plot parameters
if(!exists('UMAP.ALPHA')) {UMAP.ALPHA <<- 0.5}
if(!exists('UMAP.SIZE')) {UMAP.SIZE <<- 1}
if(!exists('LABEL.FONT.SIZE')) {LABEL.FONT.SIZE <<- 14}
if(!exists('TITLE.FONT.SIZE')) {TITLE.FONT.SIZE <<- 16}

# Load the base data - genotype, condition, and UMAP coords
if(!exists('META.DATA')) {
  META.DATA <<- load.tsv.gz(file.path(DATA.PATH, META.DATA.FILE), index.name='Cell')
}
if(!exists('GENE.META.DATA')) {
  GENE.META.DATA <<- load.tsv.gz(file.path(DATA.PATH, GENE.META.DATA.FILE))
}
if(!exists('DATA.TYPES')) {DATA.TYPES <<- c('counts', 'logcounts', 'activity')}
if(!exists('MAX.CONCURRENT.LOADED')) {MAX.CONCURRENT.LOADED <<- 50}

example.data <- readRDS("data/YEL009C.rds")

# Load the gene names that are in the data set
if(!exists('ALLOWED.NAMES')) {ALLOWED.NAMES <<- as.vector(t(GENE.META.DATA[, c("gene", "CommonName")]))}

# Set levels, labels, and colors for categorical variable TIME
if(!exists('TIME.COLORS')) {
  TIME.COLORS <<- c('#edf8b1',
                    '#c7e9b4',
                    '#7fcdbb',
                    '#41b6c4',
                    '#1d91c0',
                    '#225ea8',
                    '#253494',
                    '#081d58')
}
if(!exists('TIME.LEVELS')) {
  TIME.LEVELS <<- c(1, 2, 3, 4, 5, 6, 7, 8)
}

if(!exists('TIME.LABELS')) {
  TIME.LABELS <<- c("-20/-10", "-10/0", "0/10", "10/20", "20/30", "30/40", "40/50", "50/60")
}

# Set levels, labels, and colors for categorical variable REPLICATE
if(!exists('REPLICATE.COLORS')) {
  REPLICATE.COLORS <<- c('#7570b3', '#e6ab02')
}
if(!exists('REPLICATE.LEVELS')) {
  REPLICATE.LEVELS <<- c(1, 2)
}

if(!exists('REPLICATE.LABELS')) {
  REPLICATE.LABELS <<- c("Expt. 1", "Expt. 2")
}

# Load the figure management script (which will bring in the figure scripts)
if(!exists('get.data.plotter')) {source('figure_manager.R')}

## Here are functions to load gene-specific data and to process that data into a data.frame for plotting ##

# Take the result of a validation function and load the gene data that it refers to
get.gene.data <- function(data.list, genes.to.load) {
  
  if(is.null(genes.to.load)) {return(data.list)}
  
  for (gene.name in genes.to.load) {
    
    # Check and see if it's in the user cache
    # If it's a miss, load it from file
    if(!gene.name %in% names(data.list)) {
      data.list[[gene.name]] <- readRDS(file.path(DATA.PATH, paste0(gsub("-", "\\.", gene.name), ".rds")))
    }
  }
  
  # If the cache is too big, clear it down to just what we need here
  if (length(data.list) > MAX.CONCURRENT.LOADED) {data.list <- data.list[genes.to.load]}
  
  return(data.list)
}

## Also here are several common functions for converting between human-readable labels and labels within the data ##

yeast.systemic.to.common <- function(name.vec) {
  translated.names <- match(name.vec, GENE.META.DATA$gene, nomatch=NA)
  translated.names <- GENE.META.DATA[[translated.names, "CommonName"]]
  no.translation <- is.na(translated.names)
  translated.names[no.translation] <- name.vec[no.translation]
  return(translated.names)
}
