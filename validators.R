## Validators should use validate / need from shiny to parse inputs ##
## They should yield either a systematic, a vector of systematics   ##
## Or NULL - the output from these functions is used as file names  ##
## All validation functions should take the shiny imput object      ##

require(shiny)

# Validate a single gene
validate.gene.input <- function(input.gene){
  # Check for null input
  validate(need(!is.null(input.gene), "No gene input"))
  # Get the input string and figure out which genes it could be matching
  select.gene <- toupper(input.gene)
  could.be <- ALLOWED.NAMES[startsWith(ALLOWED.NAMES, select.gene)]
  
  # If the list of matching genes is too long, create an error message with the number of genes
  # Otherwise, number of genes and gene names
  if(length(could.be) > 50) {could.be <- c(length(could.be), "Matching Genes")}
  else{could.be <- c(length(could.be), "Matching Genes:", could.be)}
  
  # If the input gene name an be mapped to the data, proceed
  # Otherwise print the error message generated above
  select.gene <- process.gene.input(select.gene)
  validate(need(!is.null(select.gene), paste(could.be, collapse=" ")))
  
  # Return the converted input name
  return(select.gene)
}

# Validate multiple genes
validate.multigene.input <- function(input.gene) {
  gene.vec <- strsplit(toupper(input.gene), " ")[[1]]
  new.genes <- NULL
  for (gene.name in gene.vec) {new.genes <- c(new.genes, process.gene.input(gene.name))}
  shiny::validate(shiny::need(!is.null(new.genes), "At least one gene must be entered"))
  return(new.genes)
}

# Don't validate the gene input at all
validate.null <- function(input.gene) {return(NULL)}

# Process a gene name through the map file
process.gene.input <- function(gene.name) {
  gene.name <- toupper(gene.name)
  if (gene.name %in% GENE.MAP$Common) {return(toString(GENE.MAP[GENE.MAP$Common == gene.name, 'Systemic']))}
  else if (gene.name %in% GENE.MAP$Systemic) {return(gene.name)}
  else {return(NULL)}
}