require(ggplot2)
require(network)
require(GGally)

NETWORK.COLS <<- c("regulator", "target", "prior")

# Plot gene expression distributions as a violin plot
plot.figure.6.network <- function(shiny.data, gene.vec, input) {
  net.slice <- NETWORK.DATA[NETWORK.DATA$regulator %in% gene.vec | NETWORK.DATA$target %in% gene.vec, NETWORK.COLS]
  net.slice$regulator <- yeast.systemic.to.common(net.slice$regulator)
  net.slice$target <- yeast.systemic.to.common(net.slice$target)
  edge.colors <- net.slice$prior
  edge.colors[edge.colors == 0] <- "red"
  edge.colors[edge.colors == 1] <- "gray50"
  net <- network(net.slice, matrix.type = 'edgelist', directed = TRUE, ignore.eval=FALSE, names.eval=c("prior"))
  plt <- ggnet2(net, arrow.size = 12, arrow.gap = 0.025, edge.color = edge.colors, label=TRUE)
  
  return(plt)
}