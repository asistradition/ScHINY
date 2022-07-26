require(ggplot2)
require(dplyr)

# Plot expression over time of a specific gene

plot.expression.over.time <- function(species, shiny.data, gene.vec, input) {
  # Get the network data that matches inputs
  
  expr.data <- data.frame(
    expression=EXPRESSION.DATA[[species]][, gene.vec],
    time=EXPRESSION.METADATA[[species]][, EXPRESSION.METADATA.TIME],
    landrace=EXPRESSION.METADATA[[species]][, EXPRESSION.METADATA.LANDRACE],
    cond=EXPRESSION.METADATA[[species]][, EXPRESSION.METADATA.CONDITION]
  ) %>%
    dplyr::group_by(time, landrace, cond) %>%
    dplyr::summarize(expr=median(expression), n=n())
  
  # Create a real network if there's anything to plot; otherwise create nothing
  plt <- ggplot(expr.data, aes(x=factor(time), y=expr, color=landrace, linetype=cond, group=interaction(landrace, cond))) +
    geom_point() +
    stat_summary(fun=sum, geom="line") +
    theme_classic() +
    scale_y_continuous(limits=c(0, NA)) +
    labs(x="Time [min]", y=paste0(gene.vec, " Expression [TPM]"), color="Landrace", linetype="Stress")
  
  return(plt)
}

describe.expression.over.time <- function(species, ...) {
  paste0(
    SPECIES.DISPLAY.NAMES[[species]],
    " gene expression plotted over time, separated by landrace and stress condition"
  )
}