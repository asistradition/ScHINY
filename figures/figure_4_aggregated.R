require(ggplot2)
require(viridis)
require(scales)
require(patchwork)
require(dplyr)
require(tidyr)

AGGREGATE.SETTINGS <<- list()

AGGREGATE.SETTINGS[['Counts']] <- list(
  dataset='expression',
  y.label="Expression [counts]"
)

AGGREGATE.SETTINGS[['Denoised Counts']] <- list(
  dataset='denoised',
  y.label="Denoised Expression [counts]"
)

AGGREGATE.SETTINGS[['Velocity']] <- list(
  dataset='velocity',
  y.label="Velocity [counts/min]"
)

AGGREGATE.SETTINGS[['Half-life']] <- list(
  dataset='decay_constants',
  y.label='Half-life [min]',
  convert.to.halflife=T
)

AGGREGATE.SETTINGS[['TFA (Expression)']] <- list(
  dataset='tfa_expression',
  y.label='TF Activity'
)

AGGREGATE.SETTINGS[['TFA (Velocity)']] <- list(
  dataset='tfa_velocity',
  y.label='TF Activity'
)

AGGREGATE.SETTINGS[['TFA (Decay)']] <- list(
  dataset='tfa_decay_latent',
  y.label='TF Activity'
)

plot.figure.4.aggregate <- function(shiny.data, gene.vec, input) {
  
  if(is.null(input$plot_type)) {return(NULL)}
  
  gene.id <- validate.gene.input(gene.vec)
  x.label <- paste(
    get.gene.time.name(gene.id),
    "Time [min]"
  )
  
  if ("rapamycin" %in% names(shiny.data[[gene.id]])) {
    plot.data <- shiny.data[[gene.id]][['rapamycin']]
  }
  else if ("cell_cycle" %in% names(shiny.data[[gene.id]])) {
    plot.data <- shiny.data[[gene.id]][['cell_cycle']]
  }
  
  plot.dtypes <- lapply(
    input$plot_type,
    function(x) {AGGREGATE.SETTINGS[[x]]$dataset}
  )
  
  plot.data <- plot.data %>%
    dplyr::select(-gene) %>%
    dplyr::filter(dataset %in% unlist(plot.dtypes)) %>%
    tidyr::pivot_longer(
      !c(experiment, dataset, agg_func),
      names_to = "time",
      values_to = "value"
    ) %>%
    dplyr::mutate(time=as.numeric(time))
  
  n.facet <- length(input$plot_type)
  n.facet.cols <- min(n.facet, 2)
  n.facet.rows <- ceiling(n.facet / n.facet.cols)
  
  # Build a full grid with patchwork
  return(patchwork::wrap_plots(
    lapply(
      input$plot_type,
      .make.facet.plot.4,
      plot.data,
      gene.id,
      x.label
    ),
    ncol=n.facet.cols,
    nrow=n.facet.rows,
    byrow=T,
    guides='collect'
  ))
  
}

.make.facet.plot.4 <- function(plot.type, plot.data, gene.id, x.lab) {
  
  plot.settings <- AGGREGATE.SETTINGS[[plot.type]]
  
  plot.data <- dplyr::filter(plot.data, dataset == plot.settings$dataset)
  
  plot.data$alpha <- 0.75
  plot.data$alpha[plot.data$experiment == "ALL"] <- 1
  
  plot.data$dataset <- factor(plot.data$dataset)
  
  if ("convert.to.halflife" %in% names(plot.settings)) {
    plot.data[, 'value'] <- log(2) / plot.data[, 'value']
  }
  
  median.data <- plot.data[plot.data["agg_func"] == "median", ]
  std.data <- plot.data[plot.data["agg_func"] == "stdev", ]
  
  x.lim <- get.gene.time.limits(gene.id)
  y.lim <- find.data.limits(median.data[, 'value', drop=T], q=c(0, 1))
  
  line.plot <- ggplot(median.data,
                      aes(time, value, color=experiment, alpha=alpha)) +
    geom_point() +
    geom_line() +
    theme_classic() +
    labs(
      x=x.lab,
      y=plot.settings$y.label,
      color="Experiment"
    ) +
    scale_x_continuous(limits = x.lim) +
    scale_y_continuous(limits = y.lim) +
    scale_color_manual(labels = c("Both", REPLICATE.LABELS),
                       breaks = c("ALL", REPLICATE.LEVELS),
                       values = c("black", REPLICATE.COLORS)) +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          text=element_text(size=LABEL.FONT.SIZE),
          plot.title = element_text(size=TITLE.FONT.SIZE, face="bold", hjust=0.5))
  
  return(line.plot)
}

exp.panel.figure.4.aggregate <- function() {
  list(
    checkboxGroupInput(
      inputId="plot_type",
      label="Plot Type",
      choices=c(
        "Counts",
        "Denoised Counts",
        "Velocity",
        #"Half-life",
        "TFA (Expression)",
        "TFA (Velocity)",
        "TFA (Decay)"
      ),
      selected=c(
        "Denoised Counts",
        "Velocity",
        #"Half-life",
        "TFA (Decay)"
      )
    )
  )
}

describe.figure.4.aggregate <- function(...) {"Median aggregate plots from sliding windows for specific genes"}
