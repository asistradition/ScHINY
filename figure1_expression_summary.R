library(ggplot2)
library(dplyr)

# Plot visualization strings
# Units
YAXIS.SCALING.LABEL <- 'Y-Axis Units'
YAXIS.SCALING.DEFAULT <- 'Per-Cell %'

YAXIS.UMI.COUNT <- 'UMI Count'
YAXIS.UMI.COUNT.LABEL <- 'UMI Count Per Cell'

YAXIS.LIB.NORM <- 'Per-Cell %'
YAXIS.LIB.NORM.LABEL <- 'Per Cell Relative Transcript Abundance (%)'

YAXIS.WT.NORM <- 'Relative To WT'
YAXIS.WT.NORM.LABEL <- 'Per Cell Relative Transcript Abundance\n(Compared to WT)'

# Limits
YAXIS.LIMIT.LABEL <- 'Y-Axis Scaling'
YAXIS.LIMIT.DEFAULT <- 'Fixed'
YAXIS.FIXED <- 'Fixed'
YAXIS.FREE <- 'Free'

# WT Horizontal Line
YAXIS.HLINE.LABEL <- 'WT Average'
YAXIS.HLINE.DEFAULT <- 'Show'
YAXIS.SHOW <- 'Show'
YAXIS.HIDE <- 'Hide'

sidebar_expression_summary <- function() {
  list(
    checkboxGroupInput(inputId = 'conditions',
                       label = 'Conditions',
                       choices = CONDITIONS$DataColumn,
                       selected = CONDITIONS[CONDITIONS$Default,]$DataColumn),
    selectInput(inputId = 'yaxis',
                label = YAXIS.SCALING.LABEL,
                choices = c(YAXIS.UMI.COUNT, YAXIS.LIB.NORM, YAXIS.WT.NORM),
                selected = YAXIS.SCALING.DEFAULT),
    selectInput(inputId = 'yaxislim',
                label = YAXIS.LIMIT.LABEL,
                choices = c(YAXIS.FIXED, YAXIS.FREE),
                selected = YAXIS.LIMIT.DEFAULT),
    radioButtons(inputId = 'wtline',
                 label = YAXIS.HLINE.LABEL,
                 choices = c(YAXIS.SHOW, YAXIS.HIDE),
                 inline = TRUE,
                 selected = YAXIS.HLINE.DEFAULT)
  )
}

plot_expression_summary <- function(meta_data, gene, input, validator = NULL) {
  if (!is.null(validator)) {gene <- validator()}
  select.data <- meta_data
  select.data['Gene'] <- read.table(file.path(META.DATA[META.DATA$Display == input$dataset, 'Path'], paste0(gene, ".tsv")), header=FALSE)
  
  select.data %>%
    select(Gene, Condition, Genotype_Group, Genotype, TotalUMI, Num_Cells) %>%
    filter(Condition %in% input$conditions) -> select.data
  
  # Set stuff for the data scaling
  if (is.null(input$yaxis)) {y.axis.scale <- YAXIS.SCALING.DEFAULT}
  else {y.axis.scale <- input$yaxis}
  
  if (is.null(input$wtline)) {wt.line <- YAXIS.HLINE.DEFAULT}
  else {wt.line <- input$wtline}
  
  if (is.null(input$yaxislim)) {y.axis.lim <- YAXIS.LIMIT.DEFAULT}
  else {y.axis.lim <- input$yaxislim}
  
  if(y.axis.scale == YAXIS.UMI.COUNT) {
    select.data$Gene = select.data$Gene / select.data$Num_Cells
    y.text = YAXIS.UMI.COUNT.LABEL
  }
  else if (y.axis.scale == YAXIS.LIB.NORM | y.axis.scale ==  YAXIS.WT.NORM) {
    select.data$Gene <- select.data$Gene / select.data$TotalUMI * 100
    y.text = YAXIS.LIB.NORM.LABEL
  }
  
  # Calculate the WT mean
  select.data %>%
    filter(Genotype_Group %in% "WT(ho)") %>%
    group_by(Condition) %>%
    summarize(wt=mean(Gene)) -> wt_mean
  select.data <- merge(select.data, wt_mean, by='Condition')
  
  # Normalize data to WT if that option is selected
  if (y.axis.scale == YAXIS.WT.NORM) {
    y.text = YAXIS.WT.NORM.LABEL
    select.data$Gene <- select.data$Gene / select.data$wt
    select.data[is.na(select.data)] = NaN
    wt_mean$wt = 1
  }
  
  # Generate a title line with the common and systemic names
  plot.title.str <- paste("(", GENE.MAP[GENE.MAP$Systemic == gsub("\\.", "-", gene), 'Common'], ")", sep="")
  plot.title.str <- paste(gsub("\\.", "-", gene), plot.title.str, "Expression\n")
  
  # Draw plots for the data
  pl <- ggplot(select.data, aes(Genotype_Group, Gene)) +
    labs(title=plot.title.str, x="Genotype", y=y.text, color="Genotype") +
    geom_point(aes(color=factor(Genotype_Group)), size=3, alpha=0.75) +
    stat_summary(fun.y='mean', size=20, geom='point', shape='-') +
    theme_bw() +
    theme(axis.text.x = element_text(size = 14, angle=90), axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 14), legend.title = element_text(size = 16, face = 'bold'),
          plot.title = element_text(size = 16, face = "bold", hjust=0.5),
          strip.text = element_text(size = 14, face = "bold"),
          legend.position = "none")
    
    # Draw a dashed line onto the plots at the WT mean if that option is selected
    if (wt.line == YAXIS.SHOW){pl <- pl + geom_hline(data = wt_mean, aes(yintercept=wt), linetype='dashed', size=0.25)}
  
  # Facet_wrap on condition with desired scaling
  if (y.axis.lim == YAXIS.FIXED){pl <- pl + facet_wrap(~Condition, ncol=2, scales='fixed', labeller = labeller(Condition = RELABEL.FACETS))}
  else if (y.axis.lim == YAXIS.FREE){pl <- pl + facet_wrap(~Condition, ncol=2, scales='free', labeller = labeller(Condition = RELABEL.FACETS))}
  
  # Return the plot
  pl
}