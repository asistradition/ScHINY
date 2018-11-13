library(ggplot2)
library(dplyr)
library(ggridges)

# Plot visualization strings
# Units
XAXIS.SCALING.LABEL <- 'Y-Axis Units'
XAXIS.SCALING.DEFAULT <- 'Per-Cell %'

XAXIS.UMI.COUNT <- 'UMI Count'
XAXIS.UMI.COUNT.LABEL <- 'UMI Count Per Cell'

XAXIS.LIB.NORM <- 'Per-Cell %'
XAXIS.LIB.NORM.LABEL <- 'Per Cell Relative Transcript Abundance (%)'

# Limits
XAXIS.LIMIT.LABEL <- 'X-Axis Scaling'
XAXIS.LIMIT.DEFAULT <- 'Fixed'
XAXIS.FIXED <- 'Fixed'
XAXIS.FREE <- 'Free'

fake_quantile_mean <- function(data, probs=seq(0,0.5,1)) {
  return(rep(mean(data), length(probs)))
}

sidebar_expression_ridge <- function() {
  list(
    checkboxGroupInput(inputId = 'conditions',
                       label = 'Conditions',
                       choices = CONDITIONS$DataColumn,
                       selected = CONDITIONS[CONDITIONS$Default,]$DataColumn),
    selectInput(inputId = 'yaxis',
                label = XAXIS.SCALING.LABEL,
                choices = c(XAXIS.UMI.COUNT, XAXIS.LIB.NORM),
                selected = XAXIS.SCALING.DEFAULT),
    selectInput(inputId = 'yaxislim',
                label = XAXIS.LIMIT.LABEL,
                choices = c(XAXIS.FIXED, XAXIS.FREE),
                selected = XAXIS.LIMIT.DEFAULT)
  )
}

plot_expression_ridge <- function(meta_data, gene, input, validator = NULL) {
  if (!is.null(validator)) {gene <- validator()}
  select.data <- meta_data
  select.data['Gene'] <- read.table(file.path(META.DATA[META.DATA$Display == input$dataset, 'Path'], paste0(gene, ".tsv")), header=FALSE)
  
  select.data %>%
    select(Gene, Condition, Genotype_Group, Genotype, TotalUMI) %>%
    filter(Condition %in% input$conditions) -> select.data
  
  # Set stuff for the data scaling
  if (is.null(input$yaxis)) {x.axis.scale <- XAXIS.SCALING.DEFAULT}
  else {x.axis.scale <- input$yaxis}
  
  if (is.null(input$yaxislim)) {x.axis.lim <- XAXIS.LIMIT.DEFAULT}
  else {x.axis.lim <- input$yaxislim}
  
  if(x.axis.scale == XAXIS.UMI.COUNT) {
    select.data$Gene = select.data$Gene
    y.text = XAXIS.UMI.COUNT.LABEL
  }
  else if (x.axis.scale == XAXIS.LIB.NORM) {
    select.data$Gene <- select.data$Gene / select.data$TotalUMI * 100
    y.text = XAXIS.LIB.NORM.LABEL
  }
  
  select.data %>%
    group_by(Condition, Genotype) %>%
    summarize(gene_mean=mean(Gene), Genotype_Group=first(Genotype_Group), Count=n()) %>%
    as.data.frame() -> mdata
  
  # Generate a title line with the common and systemic names
  plot.title.str <- paste("(", GENE.MAP[GENE.MAP$Systemic == gsub("\\.", "-", gene), 'Common'], ")", sep="")
  plot.title.str <- paste(gsub("\\.", "-", gene), plot.title.str, "Expression\n")
  
  # Draw plots for the data
  pl <- ggplot(select.data, aes(x=Gene, y=Genotype_Group)) +
    labs(title=plot.title.str, x=y.text, y="Genotype") +
    stat_density_ridges(aes(fill=factor(Genotype_Group)), scale = 1, rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 2, quantile_fun = fake_quantile_mean) + 
    geom_point(data=mdata, aes(x=gene_mean, y=Genotype_Group, fill=factor(Genotype_Group))) +
    theme_ridges() +
    theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14, face="bold", hjust=0.5),
          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14, face="bold", hjust=0.5),
          plot.title = element_text(size = 16, face = "bold", hjust=0.5),
          strip.text = element_text(size = 16, face = "bold"),
          legend.position = 'none') +
    scale_x_continuous(limits = c(0, NA))
    
  # Facet_wrap on condition with desired scaling
  if (x.axis.lim == XAXIS.FIXED){
    pl <- pl + facet_wrap(~Condition, ncol=2, scales='fixed', labeller = labeller(Condition = RELABEL.FACETS))
  }
  else if (x.axis.lim == XAXIS.FREE){
    pl <- pl + facet_wrap(~Condition, ncol=2, scales='free', labeller = labeller(Condition = RELABEL.FACETS))
  }
  
  # Return the plot
  pl
}