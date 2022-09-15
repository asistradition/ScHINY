require(ggplot2)
require(viridis)
require(scales)
require(patchwork)

SCATTER.SETTINGS <<- list()

SCATTER.SETTINGS[['Experimental Collection Time']] <- list(
  meta.col='Pool',
  color.title="Expt. Time Pool",
  color.levels=TIME.LEVELS,
  color.labels=TIME.LABELS,
  color.palette=TIME.COLORS,
  add.vline.at.x0=T
)

SCATTER.SETTINGS[['Cell Cycle Phase']] <- list(
  meta.col='CC',
  color.title="Phase",
  color.levels=CC.LEVELS,
  color.labels=CC.LABELS,
  color.palette=CC.COLORS
)

SCATTER.SETTINGS[['Expression / Time']] <- list(
  x='time',
  y='denoised',
  x.axis.label="Time [min]",
  y.axis.label='Denoised Counts'
)

SCATTER.SETTINGS[['Velocity / Time']] <- list(
  x='time',
  y='velocity',
  x.axis.label="Time [min]",
  y.axis.label='Velocity',
  center.x.on.zero=T
)

SCATTER.SETTINGS[['Velocity / Expression']] <- list(
  x='denoised',
  y='velocity',
  x.axis.label="Denoised Counts",
  y.axis.label='Velocity',
  center.x.on.zero=T
)

plot.figure.3.scatter <- function(shiny.data, gene.vec, input) {
  
  if(is.null(input$color_by) | is.null(input$plot_type)) {return(NULL)}
  
  color.by.config <- SCATTER.SETTINGS[[input$color_by]]
  gene.id <- validate.gene.input(gene.vec)
  
  scatter.data <- data.frame(
    time = get.gene.time(gene.id),
    color = factor(META.DATA[, color.by.config$meta.col, drop=T])
  )
  
  n.facet <- length(input$plot_type)
  n.facet.cols <- min(n.facet, 2)
  n.facet.rows <- ceiling(n.facet / n.facet.cols)
  
  # Build a full grid with patchwork
  return(patchwork::wrap_plots(
   lapply(
     input$plot_type,
     .make.facet.plot,
     scatter.data,
     gene.id,
     shiny.data[[gene.id]],
     color.by.config
   ),
   ncol=n.facet.cols,
   nrow=n.facet.rows,
   byrow=T,
   guides='collect'
  ))
}

.make.facet.plot <- function(plot.type, core.data, gene.name, gene.data, color.by.config) {
  
  facet.settings <- SCATTER.SETTINGS[[plot.type]]
  
  x <- facet.settings$x
  y <- facet.settings$y
  
  # Put the data into the dataframe
  if (!(x %in% colnames(core.data))) {
    core.data[x] <- gene.data[[x]]
  }
  
  if (!(y %in% colnames(core.data))) {
    core.data[y] <- gene.data[[y]]
  }
  
  # Get limits
  if (x == 'time') {
    x.lim <- get.gene.time.limits(gene.name)
    x.lab <- paste(
      get.gene.time.name(gene.name),
      facet.settings$x.axis.label
    )
  }
  else {
    x.lim <- find.data.limits(core.data[, x])
    x.lab <- facet.settings$x.axis.label
  }
  
  y.lim <- find.data.limits(core.data[, y])

  # Make plot
  plt <- ggplot(core.data, aes_string(x, y, color='color')) +
    geom_point(alpha=UMAP.ALPHA) +
    theme_classic() +
    labs(
      x=x.lab,
      y=facet.settings$y.axis.label,
      color=facet.settings$color.title
    ) +
    scale_x_continuous(limits = x.lim) +
    scale_y_continuous(limits = y.lim) +
    scale_color_manual(labels = color.by.config$color.labels,
                       breaks = color.by.config$color.levels,
                       values = color.by.config$color.palette) +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          text=element_text(size=LABEL.FONT.SIZE),
          plot.title = element_text(size=TITLE.FONT.SIZE, face="bold", hjust=0.5)) +
    guides(color = guide_legend(override.aes = list(size=6, alpha=1)))
    
  if ('center.x.on.zero' %in% names(facet.settings)) {
    plt <- plt + 
      geom_hline(aes(yintercept = 0), color='black') +
      theme(axis.line.x = element_blank())
  }

  return(plt)
}

exp.panel.figure.3.scatter <- function() {
  list(
    checkboxGroupInput(
      inputId="plot_type",
      label="Plot Type",
      choices=c(
        "Expression / Time",
        "Velocity / Time",
        "Velocity / Expression"
      ),
      selected=c(
        "Expression / Time",
        "Velocity / Time"
      )
    ),
    selectInput(
      inputId = 'color_by',
      label = "Color By",
      choices = c(
       "Experimental Collection Time",
       "Cell Cycle Phase"
      ),
      selected = "Experimental Collection Time"
    )
  )
}

describe.figure.3.scatter <- function(...) {"Scatter plots for specific genes"}