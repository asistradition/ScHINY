## Validators should use validate / need from shiny to parse inputs ##
## They should yield either a systematic, a vector of systematics   ##
## Or NULL - the output from these functions is used as file names  ##
## All validation functions should take the shiny imput object      ##

require(shiny)
require(stringr)

# Validate a single gene
validate.gene.input <- function(input.gene, species, use.common.names=T){
  # Check for null input
  validate(need(!is.null(input.gene), "No gene input"))

  # If the input gene name an be mapped to the data, proceed
  # Otherwise print the error message generated above
  select.gene <- process.gene.input(input.gene, species, use.common.names = use.common.names)

  validate(need(length(select.gene) > 0, "Gene name not found"))
  validate(need(length(select.gene) < 2, "Multiple genes match"))

  # Return the converted input name
  return(select.gene)
}

# Validate multiple genes
validate.multigene.input <- function(input.gene, species, use.common.names=T) {
  gene.vec <- strsplit(toupper(input.gene), " ")[[1]]
  new.genes <- NULL
  for (gene.name in gene.vec) {new.genes <- c(new.genes, process.gene.input(gene.name, species, use.common.names = use.common.names))}
  shiny::validate(shiny::need(!is.null(new.genes), "At least one gene must be entered"))
  return(new.genes)
}

# Don't validate the gene input at all
validate.null <- function(input.gene, use.common.names=T) {return(NULL)}

# Process a gene name through the map file
process.gene.input <- function(gene.name, species, use.common.names=T) {
  
  gene.name <- str_replace(gene.name, coll("*"), "[A-Z0-9]*")
  gene.name <- paste0("\\b", gene.name, "\\b")

  name_matched <- str_detect(GENE.MAP[[species]]$gene_id, regex(gene.name, ignore_case = T))
  
  if (use.common.names) {
    name_matched <- name_matched | str_detect(GENE.MAP[[species]]$gene_name, regex(gene.name, ignore_case = T))
  }

  if (sum(name_matched) > 0) {
    return(GENE.MAP[[species]][name_matched, 'gene_id'])
  }

  else {return(NULL)}
}
