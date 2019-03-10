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

## Function to load data from tsv.gz files ##

load.tsv.gz <- function(file.name, ...) {read.table(base::gzfile(file.name), ...)}

# Path and file names
if(!exists('DATA.PATH')) {DATA.PATH <<- 'data'}
if(!exists('FIGURE.PATH')) {FIGURE.PATH <<- 'figures'}
if(!exists('IMAGE.PATH')) {IMAGE.PATH <<- 'images'}
if(!exists('GENE.MAP.FILE')) {GENE.MAP.FILE <<- '20181001_genes.tsv'}
if(!exists('META.DATA.FILE')) {META.DATA.FILE <<- 'jackson_2019_meta_data.tsv.gz'}
if(!exists('NETWORK.FILE')) {NETWORK.FILE <<- 'jackson_2019_network.tsv.gz'}

# Set labels and defaults for the UI
if(!exists('SHINY.TITLE')) SHINY.TITLE <<- "Single Cell Saccharomycs Cerevisiae"
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
if(!exists('META.DATA')) {META.DATA <<- load.tsv.gz(file.path(DATA.PATH, META.DATA.FILE), header=TRUE, sep="\t")}
if(!exists('NETWORK.DATA')) {NETWORK.DATA <<- load.tsv.gz(file.path(DATA.PATH, NETWORK.FILE), header=TRUE, sep="\t")}
if(!exists('DATA.TYPES')) {DATA.TYPES <<- c('counts', 'logcounts', 'activity')}
if(!exists('MAX.CONCURRENT.LOADED')) {MAX.CONCURRENT.LOADED <<- 50}

# Load the gene names that are in the data set
if(!exists('GENE.MAP')) {GENE.MAP <<- read.table(file.path(DATA.PATH, GENE.MAP.FILE), col.names = c('Systemic', 'Common'), stringsAsFactors = FALSE)}
if(!exists('ALLOWED.NAMES')) {ALLOWED.NAMES <<- as.vector(t(GENE.MAP))}

# Set levels, labels, and colors for categorical variables CONDITION and GENOTYPE
if(!exists('CONDITION.COLORS')) {CONDITION.COLORS <<- c('lightseagreen', 'darkgoldenrod1', 'darkolivegreen2', 'deepskyblue2', 'darkorchid2', 'deeppink1', 'forestgreen', 'green2', 'red4', 'palevioletred1', 'slategray3')}
if(!exists('CONDITION.LEVELS')) {CONDITION.LEVELS <<- c("CStarve", "Urea", "AmmoniumSulfate", "Proline", "Glutamine", "MinimalEtOH", "MinimalGlucose", "YPEtOH", "YPDRapa", "YPDDiauxic", "YPD")}
if(!exists('CONDITION.LABELS')) {CONDITION.LABELS <<- c("CSTARVE", "NLIM-UREA", "NLIM-NH4", "NLIM-PRO", "NLIM-GLN", "MMEtOH", "MMD", "YPEtOH", "RAPA", "DIAUXY", "YPD")}
if(!exists('CONDITION.FACET.LABELER')) {CONDITION.FACET.LABELER <<- function(string.label) {return(CONDITION.LABELS[match(string.label, CONDITION.LEVELS, nomatch=NA)])}}
if(!exists('CONDITION.LABEL.TO.LEVEL')) {CONDITION.LABEL.TO.LEVEL <<- function(string.label) {return(CONDITION.LEVELS[match(string.label, CONDITION.LABELS, nomatch=NA)])}}
if(!exists('CONDITION.SELECT.DEFAULT')) {CONDITION.SELECT.DEFAULT <<- c("NLIM-NH4", "MMD", "RAPA", "YPD")}

# Set levels, labels, and colors for categorical variables GENOTYPE
if(!exists('GENOTYPE.COLORS')) {GENOTYPE.COLORS <<- brewer.pal(12, 'Set3')}
if(!exists('GENOTYPE.LEVELS')) {GENOTYPE.LEVELS <<- c("stp2", "stp1", "rtg3", "rtg1", "gcn4",  "dal82", "dal81", "dal80", "gat1", "gln3", "gzf3", "WT(ho)")}
if(!exists('GENOTYPE.LABELS')) {GENOTYPE.LABELS <<- c("Δstp2", "Δstp1", "Δrtg3", "Δrtg1", "Δgcn4",  "Δdal82", "Δdal81", "Δdal80", "Δgat1", "Δgln3", "Δgzf3", "WT")}
if(!exists('GENOTYPE.LABEL.TO.LEVEL')) {GENOTYPE.LABEL.TO.LEVEL <<- function(string.label) {return(GENOTYPE.LEVELS[match(string.label, GENOTYPE.LABELS, nomatch=NA)])}}

# Set levels, labels, and colors for categorical variables CLUSTER
if(!exists('CLUSTER.LEVELS')) {CLUSTER.LEVELS <<- c("1", "2", "3", "4", "5", "6", "7", "8", "9")}
if(!exists('CLUSTER.COLORS')) {CLUSTER.COLORS <<- brewer.pal(9, 'Set1')}

# Load the figure management script (which will bring in the figure scripts)
if(!exists('get.data.plotter')) {source('figure_manager.R')}

## Here are functions to load gene-specific data and to process that data into a data.frame for plotting ##

# Take the result of a validation function and load the gene data that it referrs to
get.gene.data <- function(data.list, genes.to.load) {
  if(is.null(genes.to.load)) {return(data.list)}
  for (gene.name in genes.to.load) {
    if(!gene.name %in% names(data.list)) {
      data.list[[gene.name]] <- load.tsv.gz(file.path(DATA.PATH, paste0(gsub("-", "\\.", gene.name), ".tsv.gz")), header = TRUE)
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

## Also here are several common functions for converting between human-readable labels and labels within the data ##

yeast.systemic.to.common <- function(name.vec) {
  translated.names <- match(name.vec, GENE.MAP$Systemic, nomatch=NA)
  translated.names <- GENE.MAP[translated.names, "Common"]
  no.translation <- is.na(translated.names)
  translated.names[no.translation] <- name.vec[no.translation]
  return(translated.names)
}

cond.trans.list <- structure(CONDITION.LEVELS, names = CONDITION.LABELS)
geno.trans.list <- structure(GENOTYPE.LEVELS, names = GENOTYPE.LABELS)
condition.label.to.level <- function(cond.vec) {return(as.vector(cond.trans.list[cond.vec]))}
genotype.label.to.level <- function(geno.vec) {return(geno.trans.list[geno.vec])}
