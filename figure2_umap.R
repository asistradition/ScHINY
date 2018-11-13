library(ggplot2)
library(dplyr)

# Plot visualization strings
# Boundry z-scores
BOUNDRY.LABEL <- "Restrict z-scores to [-2,2]"
BOUNDRY.DEFAULT <- "True"
BOUNDRY.YES <- "True"
BOUNDRY.NO <- "False"
BOUNDRY.MAX <- 2
BOUNDRY.MIN <- -2
UMI.MAX <- 5000

# Color Options
UMAP.COLOR.LABEL <- "Color By"
UMAP.COLOR.GENE <- "Gene"
UMAP.COLOR.COND <- "Condition"
UMAP.COLOR.GENOTYPE <- "Genotype"
UMAP.COLOR.UMI <- "UMI Count"
UMAP.COLOR.DEFAULT <- "Condition"

UMAP.COLOR.COND.MAP <- structure(c("#002EFFFF", "#FF8B00FF", "#5DFF00FF", "#FF008BFF", "#00B9FFFF", "#FFFF0808", "#FF0000FF", "#00FFB9FF", "#BBBBBBBB", "#66666666", "#E800FFFF"),
                                 names = c("AmmoniumSulfate", "CStarve", "Glutamine", "MinimalEtOH", "MinimalGlucose", "Proline", "Urea", "YPD", "YPDDiauxic", "YPDRapa", "YPEtOH"))

UMAP.COLOR.GENO.MAP <- structure(c("#002EFFFF", "#FF8B00FF", "#5DFF00FF", "#FF008BFF", "#00B9FFFF", "#FFFF0808", "#FF0000FF", "#00FFB9FF", "#BBBBBBBB", "#66666666", "#E800FFFF", "#5D00FFFF"),
                                 names = c("dal80", "dal81", "dal82", "gat1", "gcn4", "gln3", "gzf3", "rtg1", "rtg3", "stp1", "stp2", "WT(ho)"))


sidebar_umap <- function() {
  list(
    radioButtons(inputId = 'umapcolorby',
                 label = YAXIS.HLINE.LABEL,
                 choices = c(UMAP.COLOR.GENE, UMAP.COLOR.COND, UMAP.COLOR.GENOTYPE, UMAP.COLOR.UMI),
                 inline = FALSE,
                 selected = UMAP.COLOR.DEFAULT),
    selectInput(inputId = 'boundry',
                label = BOUNDRY.LABEL,
                choices = c(BOUNDRY.YES, BOUNDRY.NO),
                selected = BOUNDRY.DEFAULT),
    checkboxGroupInput(inputId = 'conditions',
                       label = 'Conditions',
                       choices = CONDITIONS$DataColumn,
                       selected = CONDITIONS$DataColumn),
    checkboxGroupInput(inputId = 'genotypes',
                       label = 'Genotypes',
                       choices = names(UMAP.COLOR.GENO.MAP),
                       selected = names(UMAP.COLOR.GENO.MAP))
  )
}

plot_umap <- function(meta_data, gene, input, validator = NULL) {
  select.data <- meta_data
  
  # Sanity check color by input
  if (is.null(input$umapcolorby)) {color.by <- UMAP.COLOR.DEFAULT}
  else {color.by <- input$umapcolorby}
  
  # Set parameters based on what we want to color by
  if (color.by == UMAP.COLOR.GENE) {
    colorer = quote(Gene)
    color.label = "z-score"
    if (!is.null(validator)) {gene <- validator()}
    gene.data <- read.table(file.path(META.DATA[META.DATA$Display == input$dataset, 'Path'], paste0(gsub("-", "\\.", gene), ".tsv")), header=FALSE)
    if (input$boundry == BOUNDRY.YES) {
      gene.data[gene.data > BOUNDRY.MAX] <- BOUNDRY.MAX
      gene.data[gene.data < BOUNDRY.MIN] <- BOUNDRY.MIN
    }
    select.data['Gene'] <- gene.data$V1
    plot.title.str <- paste("(", GENE.MAP[GENE.MAP$Systemic == gsub("\\.", "-", gene), 'Common'], ")", sep="")
    plot.title.str <- paste(gsub("\\.", "-", gene), plot.title.str, "Expression\n")
  }
  else if (color.by == UMAP.COLOR.UMI) {
    colorer = quote(umi)
    color.label = "UMI Count"
    select.data['Gene'] <- 0
    plot.title.str <- "Strain Genotype"
    if (input$boundry == BOUNDRY.YES) {
      select.data[select.data$umi > UMI.MAX, "umi"] <- UMI.MAX
    }
  }
  else if (color.by == UMAP.COLOR.GENOTYPE) {colorer = quote(Genotype_Group); color.label = "Genotype"; select.data['Gene'] <- 0; plot.title.str <- "Strain Genotype"}
  else if (color.by == UMAP.COLOR.COND) {colorer = quote(Condition); color.label = "Condition"; select.data['Gene'] <- 0;  plot.title.str <- "Environmental Growth Condition"}
 
  select.data %>%
    select(Gene, Condition, Genotype_Group, V1, V2, umi) %>%
    filter(Condition %in% input$conditions) %>% 
    filter(Genotype_Group %in% input$genotypes) -> select.data

  pl <-ggplot(select.data, aes(x=V2, y=V1)) +
    geom_point(aes(color=!!colorer), alpha=0.25) +
    theme_bw() +
    labs(title=plot.title.str, x="V1", y="V2", color=color.label) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14, face="bold", hjust=0.5),
          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14, face="bold", hjust=0.5),
          plot.title = element_text(size = 16, face = "bold", hjust=0.5),
          strip.text = element_text(size = 16, face = "bold"),
          legend.text = element_text(size=14))
  
  if (color.by == UMAP.COLOR.GENE || color.by == UMAP.COLOR.UMI) {pl <- pl + scale_color_gradient2(low="blue", mid="yellow", high="red")}
  else if (color.by == UMAP.COLOR.COND) {pl <- pl + scale_color_manual(values=UMAP.COLOR.COND.MAP) + guides(color = guide_legend(override.aes = list(size=8, alpha = 1)))}
  else if (color.by == UMAP.COLOR.GENOTYPE) {pl <- pl + scale_color_manual(values=UMAP.COLOR.GENO.MAP) + guides(color = guide_legend(override.aes = list(size=8, alpha = 1)))}
  
  # Return the plot
  pl
}