require(ggplot2)

# Plot gene expression distributions as a violin plot
plot.figure.4.violin <- function(shiny.data, gene.vec, input) {
  violin.data <- process.data.list(shiny.data, gene.vec, data.type = "logcounts", scale.data = FALSE)
  violin.data <- violin.data[, c(gene.vec, 'Condition', 'Genotype_Group')]
  
  if(!is.null(input$conditions)) {use_conditions <- input$conditions}
  else {return(NULL)}
    
  if(!is.null(input$genotypes)) {use_genotypes <- input$genotypes}
  else {return(NULL)}
  
  violin.data <- violin.data[violin.data$Condition %in% condition.label.to.level(use_conditions),]
  violin.data <- violin.data[violin.data$Genotype_Group %in% genotype.label.to.level(use_genotypes),]
  violin.data <- reshape2::melt(violin.data, id.vars = c('Condition', 'Genotype_Group'), variable.name = "Gene", value.name="Count")
  
  plt <- ggplot(violin.data, aes(x=Genotype_Group, y=Count)) + 
    geom_violin(aes(fill=Genotype_Group), scale = 'width') +
    theme_classic() + 
    scale_fill_manual(labels = GENOTYPE.LABELS, breaks = GENOTYPE.LEVELS, values = GENOTYPE.COLORS) +
    scale_x_discrete(labels = use_genotypes, limits = GENOTYPE.LABEL.TO.LEVEL(use_genotypes)) +
    geom_boxplot(outlier.shape = NA, width=0.25) +
    theme(axis.text.x = element_text(size = 8, angle=90, vjust=0.5), axis.title.x = element_text(size = 12),
          axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 12),
          legend.text = element_text(size = 10), legend.title = element_text(size = 12, face = 'bold'),
          plot.title = element_text(size = 12, face = "bold", hjust=0.5),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          strip.background = element_blank(), strip.text = element_text(size=10, face="bold"),
          legend.position = "none") +
    labs(x="Genotype", y="Log2(Normalized Count)") +
    facet_wrap(. ~ Gene + Condition, labeller = as_labeller(safe.relabeler))
  return(plt)
}

exp.panel.figure.4.violin <- function() {list()}

cond.panel.figure.4.violin <- function(...) {cond.panel.standard(selected.conditions = CONDITION.SELECT.DEFAULT, selected.genotypes = GENOTYPE.LABELS)}

safe.relabeler <- function(label.names) {
  label.names <- yeast.systemic.to.common(label.names)
  relabel.names <- CONDITION.FACET.LABELER(label.names)
  relabel.names[is.na(relabel.names)] <- label.names[is.na(relabel.names)]
  return(relabel.names)
}