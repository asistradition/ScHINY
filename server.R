server <- function(input, output) {
  
  # Load the defaults the first time the server object is created
  if(!exists('shiny.data')) {
    shiny.data <- DEFAULT.SHINY.DATA
    active.data <- DATA.DEFAULT
    allowed.names <- DEFAULT.ALLOWED.NAMES
    gene.map <- DEFAULT.GENE.MAP
  }
  
  # Load new data if needed
  reactive({
    if(input$dataset != active.data) {
      shiny.data <- read.table(file.path(DATA.PATH, META.DATA[META.DATA$Display == input$dataset, 'File']))
      gene.map <- read.table(file.path(DATA.PATH, META.DATA[META.DATA$Display == input$dataset, 'GeneFile'], col.names = c('Systemic', 'Common'), stringsAsFactors = FALSE))
      allowed.names <- as.vector(t(gene.map))
      active.data <- input$dataset
    }
  })
  
  # Validate the gene input
  # Convert it to the shiny.data keys if it's in the gene map
  # Otherwise generate a matching genes error message
  validate.gene <- reactive({
    
    # Get the input string and figure out which genes it could be matching
    select.gene <- toupper(input$gene)
    could.be <- allowed.names[startsWith(allowed.names, select.gene)]
    
    # If the list of matching genes is too long, create an error message with the number of genes
    # Otherwise, number of genes and gene names
    if(length(could.be) > 50){
      could.be <- c(length(could.be), "Matching Genes")
    }
    else{
      could.be <- c(length(could.be), "Matching Genes:", could.be)
    }
    
    # If the input gene name an be mapped to the data, proceed
    # Otherwise print the error message generated above
    validate(
      need(select.gene %in% allowed.names, paste(could.be, collapse=" "))
    )
    
    # Make sure that there's one condition selected
    validate(
      need(!is.null(input$conditions), "At least one condition must be selected")
    )
    
    # Convert the input to the same (systemic) format as the Gene column in the data
    if(toupper(input$gene) %in% gene.map$Common) {
      select.gene <- toString(gene.map[gene.map$Common == select.gene,]['Systemic'])
    }
    
    # TODO: Go fix the prep script to make this unnecessary
    select.gene <- gsub("-", ".", select.gene)
    
    # Return the converted input name
    select.gene
  })
  
  # Change conditions in the UI
  output$conditions <- renderUI({
    checkboxGroupInput(inputId = 'conditions',
                       label = 'Conditions',
                       choices = levels(shiny.data$Condition),
                       selected = levels(shiny.data$Condition))
  })
  
  # Render the plot
  output$plots <- renderPlot({
    # Validate gene name
    select.gene <- validate.gene()
    
    # Filter the data to just what's needed for these plots
    shiny.data %>%
      filter(Gene == select.gene) %>%
      filter(Condition %in% input$conditions) -> select.data
    
    # Set stuff for the data scaling
    y.axis.scale <- input$yaxis
    if(y.axis.scale == YAXIS.UMI.COUNT) {
      y.quote = quote(Count)
      y.text = YAXIS.UMI.COUNT.LABEL
    }
    else if (y.axis.scale == YAXIS.LIB.NORM | y.axis.scale ==  YAXIS.WT.NORM) {
      y.quote = quote(LibraryNormalized)
      y.text = YAXIS.LIB.NORM.LABEL
    }
    
    # Calculate the WT mean
    select.data %>%
      filter(Gene_Group %in% "WT") %>%
      group_by(Condition) %>%
      summarize(wt=mean(!!y.quote)) -> wt_mean
    
    # Normalize data to WT if that option is selected
    if (input$yaxis == YAXIS.WT.NORM) {
      y.quote = quote(FoldChange)
      y.text = YAXIS.WT.NORM.LABEL
      select.data <- merge(select.data, wt_mean, by='Condition')
      select.data$FoldChange <- select.data$LibraryNormalized / select.data$wt
      select.data[is.na(select.data)] = NaN
      wt_mean$wt = 1
    }
    
    # Generate a title line with the common and systemic names
    plot.title.str <- paste("(", gene.map[gene.map$Systemic == gsub("\\.", "-", select.gene), 'Common'], ")", sep="")
    plot.title.str <- paste(gsub("\\.", "-", select.gene), plot.title.str, "Expression\n")
    
    # Draw plots for the data
    pl <- ggplot(select.data, aes_q(quote(factor(Gene_Group)), y.quote)) +
      labs(title=plot.title.str, x="Genotype", y=y.text, color="Genotype") +
      geom_point(aes(color=factor(Gene_Group)), size=3, alpha=0.75) +
      stat_summary(fun.y='mean', size=20, geom='point', shape='-') +
      theme_bw() + 
      theme(axis.text.x = element_text(size = 12, angle=90), axis.title.x = element_text(size = 14),
            axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14),
            legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
            plot.title = element_text(size = 16, face = "bold", hjust=0.5)) +
      
    # Draw a dashed line onto the plots at the WT mean if that option is selected
    if (input$wtline == YAXIS.SHOW){
      pl <- geom_hline(data = wt_mean, aes(yintercept=wt), linetype='dashed', size=0.25)
    }
    
    # Facet_wrap on condition with desired scaling
    if (input$yaxislim == YAXIS.FIXED){
      pl <- pl + facet_wrap(~Condition, ncol=1, scales='fixed')
    }
    else if (input$yaxislim == YAXIS.FREE){
      pl <- pl + facet_wrap(~Condition, ncol=1, scales='free')
    }
    
    # Return the plot
    pl
  })
}