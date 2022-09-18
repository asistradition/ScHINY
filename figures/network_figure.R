require(ggplot2)
require(network)
require(GGally)

NETWORK.COLS <<- c("regulator", "target", "prior")

# Plot networks
plot.figure.network <- function(network.df, species, shiny.data, gene.vec, input) {
  # Get the network data that matches inputs
  
  net.slice <- acquire.network.data(network.df, gene.vec, species, use.common.names = input$common_names)

  # Create a real network if there's anything to plot; otherwise create nothing
  if (nrow(net.slice) > 0) {
    edge.colors <- net.slice$prior
    edge.colors[edge.colors == 0] <- "red"
    edge.colors[edge.colors == 1] <- "gray50"
    
    net <- 
    plt <- ggnet2(
      network(
        net.slice,
        matrix.type = 'edgelist',
        directed = TRUE,
        ignore.eval=FALSE,
        loops = F
      ),
      arrow.size = 12,
      arrow.gap = 0.025,
      edge.color = edge.colors,
      label=TRUE
    )
  }
  
  else {
    if (input$common_names) {
      net.slice <- data.frame(
        regulator=systemic.to.common(gene.vec, GENE.MAP[[species]]),
        target=systemic.to.common(gene.vec, GENE.MAP[[species]]),
        stringsAsFactors=FALSE
      )
    }
    else {
      net.slice <- data.frame(
        regulator=gene.vec,
        target=gene.vec,
        stringsAsFactors=FALSE
      )
    }

    plt <- ggnet2(
      network(
        net.slice,
        matrix.type = 'edgelist',
        directed = FALSE,
        ignore.eval = FALSE,
        loops = T
      ),
      label=TRUE,
      edge.alpha = 0
    )
    
    # The force algorithms return NaN when there's 1 thing. Fix that.
    if (length(gene.vec) == 1) {plt$data$x = c(0); plt$data$y = c(0)}
  }
  
  return(plt)
}

acquire.network.data <- function(network.df, gene.vec, species, use.common.names=T) {
  net.slice <- network.df[network.df$regulator %in% gene.vec | network.df$target %in% gene.vec, NETWORK.COLS]
  net.slice <- net.slice[net.slice$regulator != net.slice$target, ]
  
  if (nrow(net.slice) > 0) {
    if (use.common.names) {
      net.slice$regulator <- systemic.to.common(net.slice$regulator, GENE.MAP[[species]])
      net.slice$target <- systemic.to.common(net.slice$target, GENE.MAP[[species]])   
    }
    net.slice <- net.slice[!duplicated(net.slice[,c("regulator", "target")]), ]
  }
  
  return(net.slice)
}

describe.figure.network <- function(net.df, species, ...) {
  paste0(
    SPECIES.DISPLAY.NAMES[[species]],
    " network graph with ",
    toString(nrow(net.df)),
    " regulatory edges; known interaction edges from the prior in gray ",
    "and new inferred interaction edges in red"
  )
}