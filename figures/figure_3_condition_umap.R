require(ggplot2)
require(viridis)
require(scales)

plot.figure.3.umap <- function(shiny.data, gene.vec, input) {
  
  if(is.null(input$color_by) || is.null(input$condition) || is.null(input$genotypes)) {return(NULL)}
  
  if(input$color_by == "Cluster") {
    color_by <- "Cluster"; color.title <- "Cluster"; plot.title <- "Condition"
    umap.data <- process.data.list(shiny.data, NULL)
  }
  else if(input$color_by == "Genotype") {
    color_by <- "Genotype_Group"; color.title <- "Genotype"; plot.title <- "Genotype"
    umap.data <- process.data.list(shiny.data, NULL)
  }
  else if(input$color_by == "Gene Expression") {
    color_by <- validate.gene.input(gene.vec)
    color.title <- "Log2(Expr)"
    plot.title <- paste0("Log2(Transcript Count) Gene Expression [", yeast.systemic.to.common(color_by), "]")
    umap.data <- process.data.list(shiny.data, color_by, data.type = "counts", scale.data = FALSE)
  }
  else {return(NULL)}
  
  umap.data <- umap.data[umap.data$Condition == condition.label.to.level(input$condition), ]
  umap.data <- umap.data[umap.data$Genotype_Group %in% genotype.label.to.level(input$genotypes),]
  
  if (color_by == "Genotype_Group") {umap.data <- umap.data[, c('Genotype_Group', 'UMAP.FIG3.1', 'UMAP.FIG3.2')]}
  else if (color_by == "Cluster") {
    umap.data <- umap.data[, c('Cluster', 'UMAP.FIG3.1', 'UMAP.FIG3.2')]
    umap.data$Cluster <- as.factor(umap.data$Cluster)
  }
  else {
    umap.data <- umap.data[, c(color_by, 'Condition', 'Genotype_Group', 'UMAP.FIG3.1', 'UMAP.FIG3.2')]
    umap.data[color_by] <- log2(umap.data[color_by])
  }
  
  ## Create base UMAP as a point plot ##
  umap.data <- umap.data[sample(nrow(umap.data)),]
  umap.plt <- ggplot(umap.data, aes(x=UMAP.FIG3.1, y=UMAP.FIG3.2)) + 
    geom_point(aes_string(color=color_by), alpha=UMAP.ALPHA, size=UMAP.SIZE) + 
    theme_classic() + labs(color=color.title) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          text=element_text(size=LABEL.FONT.SIZE), plot.title = element_text(size=TITLE.FONT.SIZE, face="bold", hjust=0.5)) +
    labs(title = plot.title, x="UMAP1", y="UMAP2")
  
  ## Add color scale ##
  
  if(is.null(input$color_by) || input$color_by == "Cluster") {umap.plt <- umap.plt + 
    scale_color_manual(labels = CLUSTER.LEVELS, breaks = CLUSTER.LEVELS, values = structure(CLUSTER.COLORS, names = CLUSTER.LEVELS)) +
    guides(color = guide_legend(override.aes = list(size=6, alpha=1)))}
  else if(input$color_by == "Genotype") {umap.plt <- umap.plt + 
    scale_color_manual(labels = GENOTYPE.LABELS, breaks = GENOTYPE.LEVELS, values = GENOTYPE.COLORS) +
    guides(color = guide_legend(override.aes = list(size=6, alpha=1), reverse = TRUE))}
  else if(input$color_by == "Gene Expression") {
    scale.max <- quantile(umap.data[, color_by], probs = 0.99)[[1]]
    umap.plt <- umap.plt + scale_colour_viridis(option='B', limits = c(0,scale.max), oob = squish)
  }
  return(umap.plt)
}

exp.panel.figure.3.umap <- function() {
  list(selectInput(inputId = 'condition',
                   label = 'Condition',
                   choices = CONDITION.LABELS,
                   selected = "YPD"),
       selectInput(inputId = 'color_by',
                   label = "Color By",
                   choices = c("Cluster", "Genotype", "Gene Expression"),
                   selected = "Cluster"))
}

cond.panel.figure.3.umap <- function(...) {list(checkboxGroupInput(inputId = 'genotypes', label = 'Genotypes', choices = rev(GENOTYPE.LABELS), selected = rev(GENOTYPE.LABELS)))}

describe.figure.3.umap <- function(...) {"Cells from each growth condition are separately normalized and transformed into 2-dimensional space by UMAP."}
