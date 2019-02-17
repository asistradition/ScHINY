require(ggplot2)

FIGURE.1.FILE.NAME <- "Figure_1-01.png"

# Don't judge me
plot.figure.1.schematic <- function(shiny.data, gene, input) {
  fig1.img <- grid::rasterGrob(png::readPNG(file.path(IMAGE.PATH, FIGURE.1.FILE.NAME)), interpolate=TRUE)
  fig1.plt.hack <- qplot(1:10, 1:10, geom = "blank") +
    annotation_custom(fig1.img, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    theme_void()
  return(fig1.plt.hack)
}

exp.panel.figure.1.schematic <- function() {list()}
