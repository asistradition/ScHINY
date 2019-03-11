require(ggplot2)
require(viridis)
require(scales)

plot.figure.2.umap <- function(shiny.data, gene.vec, input) {
  
  if(is.null(input$color_by) || is.null(input$conditions) || is.null(input$genotypes)) {return(NULL)}
  
  if(input$color_by == "Condition") {
    color_by <- "Condition"; color.title <- "Condition"; plot.title <- "Condition"
    umap.data <- process.data.list(shiny.data, NULL)
    }
  else if(input$color_by == "Genotype") {
    color_by <- "Genotype_Group"; color.title <- "Genotype"; plot.title <- "Genotype"
    umap.data <- process.data.list(shiny.data, NULL)
    }
  else if(input$color_by == "Gene Expression") {
    color_by <- validate.gene.input(gene.vec)
    color.title <- "Log2(Expr)"
    plot.title <- paste0("Log2(\"MultiBatch Normalized\") Gene Expression [", yeast.systemic.to.common(color_by), "]")
    umap.data <- process.data.list(shiny.data, color_by, data.type = "logcounts", scale.data = FALSE)
  }
  else {return(NULL)}

  if (color_by %in% c("Condition", "Genotype_Group")) {umap.data <- umap.data[, c('Condition', 'Genotype_Group', 'UMAP1', 'UMAP2')]}
  else {umap.data <- umap.data[, c(color_by, 'Condition', 'Genotype_Group', 'UMAP1', 'UMAP2')]}
  
  umap.data <- umap.data[umap.data$Condition %in% condition.label.to.level(input$conditions),]
  umap.data <- umap.data[umap.data$Genotype_Group %in% genotype.label.to.level(input$genotypes),]
  
  ## Create base UMAP as a point plot ##
  umap.data <- umap.data[sample(nrow(umap.data)),]
  umap.plt <- ggplot(umap.data, aes(x=UMAP1, y=UMAP2)) + 
    geom_point(aes_string(color=color_by), alpha=UMAP.ALPHA, size=UMAP.SIZE) + 
    theme_classic() + labs(color=color.title) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          text=element_text(size=LABEL.FONT.SIZE), plot.title = element_text(size=TITLE.FONT.SIZE, face="bold", hjust=0.5)) +
    labs(title = plot.title)
  
  ## Add color scale ##
  
  if(is.null(input$color_by) || input$color_by == "Condition") {umap.plt <- umap.plt + 
    scale_color_manual(labels = CONDITION.LABELS, breaks = CONDITION.LEVELS, values = CONDITION.COLORS) +
    guides(color = guide_legend(override.aes = list(size=6, alpha=1), reverse = TRUE))}
  else if(input$color_by == "Genotype") {umap.plt <- umap.plt + 
    scale_color_manual(labels = GENOTYPE.LABELS, breaks = GENOTYPE.LEVELS, values = GENOTYPE.COLORS) +
    guides(color = guide_legend(override.aes = list(size=6, alpha=1), reverse = TRUE))}
  else if(input$color_by == "Gene Expression") {
    scale.max <- quantile(umap.data[, color_by], probs = 0.99)[[1]]
    umap.plt <- umap.plt + scale_colour_viridis(option='B', limits = c(0,scale.max), oob = squish)
  }
  return(umap.plt)
}

exp.panel.figure.2.umap <- function() {
  list(selectInput(inputId = 'color_by',
                   label = "Color By",
                   choices = c("Condition", "Genotype", "Gene Expression"),
                   selected = "Condition"))
}

cond.panel.figure.2.umap <- function(...) {cond.panel.standard(selected.conditions = CONDITION.LABELS, selected.genotypes = GENOTYPE.LABELS)}

describe.figure.2.umap <- function(...) {"Uniform Manifold Approximation and Projection (UMAP) projection of log-transformed and batch-normalized single-cell data"}