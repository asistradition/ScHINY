require(ggplot2)
require(viridis)
require(scales)

UMAP.SETTINGS <<- list()

UMAP.SETTINGS[['Experimental Replicate']] <- list(
  meta.col='Experiment',
  color.title="Replicate",
  color.levels=REPLICATE.LEVELS,
  color.labels=REPLICATE.LABELS,
  color.palette=REPLICATE.COLORS,
  plot.title="Experimental Replicate"
)

UMAP.SETTINGS[['Experimental Collection Time']] <- list(
  meta.col='Pool',
  color.title="Expt. Time Pool",
  color.levels=TIME.LEVELS,
  color.labels=TIME.LABELS,
  color.palette=TIME.COLORS,
  plot.title="Experimental Collection Time [Pools]"
)

UMAP.SETTINGS[['Cell Cycle Phase']] <- list(
  meta.col='CC',
  color.title="Phase",
  color.levels=CC.LEVELS,
  color.labels=CC.LABELS,
  color.palette=CC.COLORS,
  plot.title="Estimated Cell Cycle Phase [Spellman 1998]"
)

UMAP.SETTINGS[['Rapamycin Response Time']] <- list(
  meta.col='program_0_time',
  color.title="Rapa Time",
  color.bounds=c(-10, 60),
  color.palette="mako",
  plot.title="Rapamycin Response Time [min]"
)

UMAP.SETTINGS[['Cell Cycle Time']] <- list(
  meta.col='program_1_time',
  color.title="CC Time",
  color.bounds=c(0, 88),
  color.palette="turbo",
  plot.title="Cell Cycle Time [min]"
)

UMAP.SETTINGS[['Gene Expression']] <- list(
  layer='expression',
  color.title="Counts",
  color.palette="viridis",
  plot.title=" Gene Expression [counts]",
  plot.title.paste.gene=T
)

UMAP.SETTINGS[['Denoised Expression']] <- list(
  layer='denoised',
  color.title="Counts",
  color.palette="viridis",
  plot.title=" Gene Expression [counts]",
  plot.title.paste.gene=T
)

UMAP.SETTINGS[['Velocity']] <- list(
  layer='velocity',
  color.title="Counts",
  color.palette=c(muted("red"), "white", muted("blue")),
  plot.title=" Gene Velocity [counts/min]",
  plot.title.paste.gene=T,
  gradient.palette=T
)

UMAP.COLS <<- c("UMAP_1", "UMAP_2")

plot.figure.2.umap <- function(shiny.data, gene.vec, input) {
  
  if(is.null(input$color_by)) {return(NULL)}
  
  color.by.config <- UMAP.SETTINGS[[input$color_by]]
  
  # Use the metadata col if that's what we're plotting
  if (!is.null(color.by.config$meta.col)) {
    color.by.col <- color.by.config$meta.col
    umap.data <- META.DATA[, c(color.by.col, UMAP.COLS)]
  }
  else if (!is.null(color.by.config$layer)) {
    color.by.col <- validate.gene.input(gene.vec)
    umap.data <- META.DATA
    umap.data[, color.by.col] <- shiny.data[[color.by.col]][[color.by.config$layer]]
  }
  
  # Add gene name to title if that flag is set
  if (!is.null(color.by.config$plot.title.paste.gene)) {
    plot.title <- paste0(yeast.systemic.to.common(color.by.col), color.by.config$plot.title)
  }
  else {
    plot.title <- color.by.config$plot.title
  }
  
  ## Create base UMAP data frame and shuffle it to prevent overplotting ##
  umap.data <- as.data.frame(umap.data[sample(nrow(umap.data)), ])
  
  ## Factorize if needed ##
  if('color.levels' %in% names(color.by.config)) {
    umap.data[, color.by.col] <- as.factor(umap.data[, color.by.col])
  }
  
  ## Build the core of the plot ##
  umap.plt <- ggplot(umap.data, aes(x=UMAP_1, y=UMAP_2)) + 
    geom_point(aes_string(color=color.by.col), alpha=UMAP.ALPHA, size=UMAP.SIZE) + 
    theme_classic() +
    labs(color=color.by.config$color.title) +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          text=element_text(size=LABEL.FONT.SIZE),
          plot.title = element_text(size=TITLE.FONT.SIZE, face="bold", hjust=0.5)) +
    labs(title = plot.title)
  
  ## Add color scale ##
  
  # Use manual scale
  if('color.levels' %in% names(color.by.config)) {
    umap.plt <- umap.plt + 
      scale_color_manual(labels = color.by.config$color.labels,
                         breaks = color.by.config$color.levels,
                         values = color.by.config$color.palette) +
    guides(color = guide_legend(override.aes = list(size=6, alpha=1)))}
  
  # Use the defined limits for a viridis scale
  else if('color.bounds' %in% names(color.by.config)) {
    umap.plt <- umap.plt + 
      scale_colour_viridis(
        limits=color.by.config$color.bounds,
        option=color.by.config$color.palette,
        oob=squish)
  }
  
  # Make our own 0.01 / 0.99 limits for a viridis scale
  else {
    
    lims <- quantile(umap.data[, color.by.col], c(0.01, 0.99))
    
    # If the lower limit is > 0 bound it at zero
    if (lims[1] > 0) {lims[1] <- 0}
    
    # If the lower limit is negative and the upper limit is positive
    # Bound on abs
    if (lims[1] < 0 & lims[2] > 0) {
      lims[1] <- -1 * max(abs(lims))
      lims[2] <- max(abs(lims))
    }
    
    if('color.log.scale' %in% names(color.by.config)) {
      umap.plt <- umap.plt + 
        scale_colour_viridis(
          limits=lims,
          option=color.by.config$color.palette,
          trans = scales::pseudo_log_trans(sigma = 1, base=2),
          oob=squish)      
    }
    else if ('gradient.palette' %in% names(color.by.config)) {
      umap.plt <- umap.plt + 
        scale_colour_gradient2(
          limits=lims,
          low=color.by.config$color.palette[1],
          mid=color.by.config$color.palette[2],
          high=color.by.config$color.palette[3],
          oob=squish)
    }
    else {
      umap.plt <- umap.plt + 
        scale_colour_viridis(
          limits=lims,
          option=color.by.config$color.palette,
          oob=squish)
    }
  }

  return(umap.plt)
}

exp.panel.figure.2.umap <- function() {
  list(selectInput(inputId = 'color_by',
                   label = "Color By",
                   choices = c(
                     "Experimental Collection Time",
                     "Experimental Replicate",
                     "Rapamycin Response Time",
                     "Cell Cycle Time",
                     "Cell Cycle Phase",
                     "Gene Expression",
                     "Denoised Expression",
                     "Velocity"
                   ),
                   selected = "Experimental Collection Time"))
}

describe.figure.2.umap <- function(...) {"Uniform Manifold Approximation and Projection (UMAP) projection"}