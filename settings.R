VERSION.STRING <- "July 21, 2022"

## This file should be sourced at the top of all the shiny app scripts ##
## It loads necessary packages and populates the global namespace ##

# Load shiny
library(shiny)

# Load plotting stuff
library(ggplot2)
library(RColorBrewer)
library(network)
library(GGally)
library(data.table)
library(stringr)

# Load tidyverse stuff
library(dplyr)

source("rice_parse_settings.R")

## Function to load data from tsv.gz files ##

load.tsv.gz <- function(file.name, ...) {read.table(base::gzfile(file.name), ...)}
cut.max.mcc <- function(net.df) {
  n <- which(net.df[, "MCC"] == max(net.df[, "MCC"], na.rm = T))
  return(net.df[1:n, ])
}

# Path and file names
if(!exists('DATA.PATH')) {DATA.PATH <<- 'data'}
if(!exists('FIGURE.PATH')) {FIGURE.PATH <<- 'figures'}
if(!exists('IMAGE.PATH')) {IMAGE.PATH <<- 'images'}

if(!exists('NETWORK.RDS.FILE')) {NETWORK.RDS.FILE <<- 'Oryza_Networks_072722.rds'}
if(!exists('EXPRESSION.RDS.FILE')) {EXPRESSION.RDS.FILE <<- 'MZ2020_expression_072722.rds'}


if(!exists('SPECIES_DISPLAY_NAMES')) {
  SPECIES.DISPLAY.NAMES <- list(OSI='Oryza sativa indica',
                                OSJ='Oryza sativa japonica',
                                OG='Oryza glaberrima')
}

if(!exists('GENE.MAP.FILE')) {
  GENE.MAP.FILE <<- list(OSI='Oryza_sativa_indica_gene_metadata.tsv',
                         OSJ='Oryza_sativa_japonica_gene_metadata.tsv',
                         OG='Oryza_glaberrima_gene_metadata.tsv')
  }

if(!exists('NETWORK.FILE')) {
  NETWORK.FILE <<- list(
    OSI_BBSR_EGRIN='Oryza_sativa_indica_network_BBSR_EGRIN.tsv.gz',
    OSJ_BBSR_EGRIN='Oryza_sativa_japonica_network_BBSR_EGRIN.tsv.gz',
    OG_BBSR_EGRIN='Oryza_glaberrima_network_BBSR_EGRIN.tsv.gz',
    OSI_BBSR_PREGDB='Oryza_sativa_indica_network_BBSR_PREGDB.tsv.gz',
    OSJ_BBSR_PREGDB='Oryza_sativa_japonica_network_BBSR_PREGDB.tsv.gz',
    OG_BBSR_PREGDB='Oryza_glaberrima_network_BBSR_PREGDB.tsv.gz'
  )
}

if(!exists('EXPRESSION.FILE')) {
  EXPRESSION.FILE <<- list(OSI='MZ2020_Oryza_indica_tpm.tsv.gz',
                            OSJ='MZ2020_Oryza_sativa_tpm.tsv.gz',
                            OG='MZ2020_Oryza_glaberrima_tpm.tsv.gz')
}

if(!exists('EXPRESSION.METADATA.FILE')) {
  EXPRESSION.METADATA.FILE <<- "metadata_MZ2021.tsv"
}

if(!exists('EXPRESSION.METADATA.TIME')) {
  EXPRESSION.METADATA.TIME <<- "Time"
}
if(!exists('EXPRESSION.METADATA.LANDRACE')) {
  EXPRESSION.METADATA.LANDRACE <<- "Genotype_Name"
}
if(!exists('EXPRESSION.METADATA.CONDITION')) {
  EXPRESSION.METADATA.CONDITION <<- "Condition"
}



# Set labels and defaults for the UI
if(!exists('SHINY.TITLE')) SHINY.TITLE <<- "Rice Networks"
if(!exists('GENE.LABEL')) GENE.LABEL <<- 'Gene Name'
if(!exists('GENE.DEFAULT')) GENE.DEFAULT <<- 'OsLG*'
if(!exists('FIGURE.SELECT.LABEL')) {FIGURE.SELECT.LABEL <<- "Figure"}
if(!exists('DEFAULT.FIGURE')) {DEFAULT.FIGURE <<- "Select Figure"}

# Set plot parameters
if(!exists('LABEL.FONT.SIZE')) {LABEL.FONT.SIZE <<- 14}
if(!exists('TITLE.FONT.SIZE')) {TITLE.FONT.SIZE <<- 16}

# Load network data
if(!exists('NETWORK.DATA')) {
  if(file.exists(file.path(DATA.PATH, NETWORK.RDS.FILE))) {
    NETWORK.DATA <<- readRDS(file.path(DATA.PATH, NETWORK.RDS.FILE))
  }
  else {
    NETWORK.DATA <<- lapply(
      X = NETWORK.FILE,
      FUN = function(x) {
        cut.max.mcc(
          load.tsv.gz(
            file.path(DATA.PATH, x),
            header=TRUE,
            sep="\t",
            stringsAsFactors = FALSE
          )
        )
      }
    )
  }
}

# Load expression data
if(!exists('EXPRESSION.DATA')) {
  if(file.exists(file.path(DATA.PATH, EXPRESSION.RDS.FILE))) {
    EXPRESSION.DATA <<- readRDS(file.path(DATA.PATH, EXPRESSION.RDS.FILE))
  }
  else {
    EXPRESSION.DATA <<- lapply(
      X = EXPRESSION.FILE,
      FUN = function(x) {
        load.tsv.gz(
          file.path(DATA.PATH, x),
          header=TRUE,
          sep="\t",
          stringsAsFactors = FALSE
        )
      }
    )
  }
}

if(!exists('EXPRESSION.METADATA')) {
  EXPRESSION.METADATA <<- lapply(
    X = EXPRESSION.DATA,
    FUN = function(x) {
      md <- read.table(
        file.path(DATA.PATH, EXPRESSION.METADATA.FILE),
        header=TRUE,
        sep="\t",
        stringsAsFactors = FALSE
      )
      md <- md[md$Sample %in% x$X,]
      rownames(md) <- md$Sample
      colnames(md) <- str_replace_all(colnames(md), "\\.{1,}", "_")
      
      md <- md[x$X, ]
      md <- cbind(md, parse.sample.mz2022(rownames(md)))
      
      return(md)
    }
  )
}

# Load the gene names that are in the data set
if(!exists('GENE.MAP')) {
  GENE.MAP <<- lapply(
    GENE.MAP.FILE,
    function (x) {
      as.data.frame(fread(file.path(DATA.PATH, x), stringsAsFactors = F))[, c("gene_id", "gene_name")]
    }
  )
}

if(!exists('ALLOWED.NAMES')) {
  ALLOWED.NAMES <<- lapply(
    GENE.MAP,
    function(x) {as.vector(t(x[, c("gene_id", "gene_name")]))}
  )
}

# Load the figure management script (which will bring in the figure scripts)
if(!exists('get.data.plotter')) {source('figure_manager.R')}

systemic.to.common <- function(name.vec, map.table) {
  
  if (length(name.vec) == 0) {
    return(name.vec)
  }
  
  translated.names <- match(name.vec, map.table$gene_id, nomatch=NA)
  translated.names <- map.table[translated.names, "gene_name"]
  no.translation <- is.na(translated.names)
  translated.names[no.translation] <- name.vec[no.translation]
  return(translated.names)
}
