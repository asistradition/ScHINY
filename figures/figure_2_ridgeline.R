require(ggplot2)
require(cowplot)

plot.figure.2.ridgeline <- function(shiny.data, gene.vec, input) {
  ridgeline.data <- process.data.list(shiny.data, gene.vec, data.type = "counts", scale.data = FALSE)
  ridgeline.data <- ridgeline.data[, c(gene.vec, 'Condition', 'Genotype_Group')]
  ridgeline.data <- reshape2::melt(ridgeline.data, id.vars = c('Condition', 'Genotype_Group'), variable.name = "Gene", value.name="Count")
  
  # Calculate the median expression values for each strain (TF deletion)
  ridgeline.medians <- dplyr::group_by(ridgeline.data, Gene, Condition, Genotype_Group)
  ridgeline.medians <- dplyr::summarize(ridgeline.medians, gene_median=median(Count))
  
  # Make a faceted ridgeline plot for each gene
  plt.list <- lapply(as.list(gene.vec), function(x) {make.ridgeline.plt(ridgeline.data, ridgeline.medians, x) + labs(title=yeast.systemic.to.common(x))})
  
  # Pull the first plot off and hack it into just being labels
  plt.labels <- plt.list[[1]] + 
    scale_x_continuous(limits = c(0, 0), breaks = c(0)) + 
    theme_void() +
    theme(plot.title = element_blank(), axis.text.y = element_text(size=LABEL.FONT.SIZE, vjust=0, hjust=1, face="bold", color="black"),
          axis.title.x = element_blank(), legend.position = 'none', plot.margin = unit(c(0, 0, 0, 0), "cm"))
  plt.labels$layers[[2]] <- NULL
  
  # Stack up the label plots with some whitespace on the top and bottom
  # Stack up the ridgeline facet plots left to right with the labels on the left
  plot.stack <- plot_grid(plot_grid(ggdraw(), plt.labels, ggdraw(), ncol = 1, align="hv", rel_heights = c(.025, 1, .05)), 
                          plot_grid(plotlist = plt.list, ncol = length(plt.list), align="hv"), 
                          ncol=2, rel_widths = c(.1, 1), align="hv")
  
  # Stack up all the plots with an x-axis title on the bottom
  x.axis.title <- ggdraw() + draw_label("Transcript Counts", fontface='bold')
  return(plot_grid(plot.stack, x.axis.title, ncol=1, rel_heights = c(1,0.05)))
}

exp.panel.figure.2.ridgeline <- function() {list()}

# Plot function with the common options for each plot
make.ridgeline.plt <- function(ridge.data, median.data, gene.name) {
  ridge.data <- ridge.data[ridge.data$Gene == gene.name, ]
  x.max <- max(quantile(ridge.data$Count, probs = 0.99)[[1]], 1)
  ggplot(ridge.data, aes(x=Count, y=Condition)) +
    stat_binline(aes(height=..ndensity.., fill=factor(Condition)), binwidth=1, center=0, draw_baseline = FALSE, scale = 0.9, na.rm = TRUE) + 
    geom_point(data=median.data[median.data$Gene == gene.name, ], aes(x=gene_median, y=Condition), na.rm = TRUE) +
    theme_ridges(grid = FALSE) +
    scale_y_discrete(limits=CONDITION.LEVELS, labels=CONDITION.LABELS) +
    scale_x_continuous(limits = c(-0.5, x.max), breaks = function(x) {ceiling(scales::pretty_breaks(n=3)(x))}) +
    scale_fill_manual(values=CONDITION.COLORS) +
    theme(axis.text.x = element_text(size = LABEL.FONT.SIZE, angle=90, vjust = 0.5, hjust = 0.5),
          plot.title = element_text(size = TITLE.FONT.SIZE, face = "bold.italic", hjust=0.5), axis.text.y = element_blank(), axis.title.y = element_blank(), 
          axis.title.x = element_blank(), legend.position = 'none', plot.margin = unit(c(0, 0, 0, 0), "cm"))
}