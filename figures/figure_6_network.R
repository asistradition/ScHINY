require(ggplot2)
require(network)
require(GGally)

NETWORK.COLS <<- c("regulator", "target", "prior")

# Plot gene expression distributions as a violin plot
plot.figure.6.network <- function(shiny.data, gene.vec, input) {
  # Get the network data that matches inputs
  net.slice <- NETWORK.DATA[NETWORK.DATA$regulator %in% gene.vec | NETWORK.DATA$target %in% gene.vec, NETWORK.COLS]
  
  # Create a real network if there's anything to plot; otherwise create a adjacency matrix of 0s
  if (nrow(net.slice) > 0) {
    net.slice$regulator <- yeast.systemic.to.common(net.slice$regulator)
    net.slice$target <- yeast.systemic.to.common(net.slice$target)
    edge.colors <- net.slice$prior
    edge.colors[edge.colors == 0] <- "red"
    edge.colors[edge.colors == 1] <- "gray50"
    net <- network(net.slice, matrix.type = 'edgelist', directed = TRUE, ignore.eval=FALSE, names.eval=c("prior"))
    plt <- ggnet2(net, arrow.size = 12, arrow.gap = 0.025, edge.color = edge.colors, label=TRUE)
  }
  else {
    net.slice <- data.frame(regulator=yeast.systemic.to.common(gene.vec), target=yeast.systemic.to.common(gene.vec))
    net <- network(net.slice, matrix.type = 'edgelist', directed = FALSE)
    plt <- ggnet2(net, label=TRUE, edge.alpha = 0)
    
    # The force algorithms return NaN when there's 1 thing. Fix that.
    if (length(gene.vec) == 1) {plt$data$x = c(0); plt$data$y = c(0)}
  }
  
  
  return(plt)
}