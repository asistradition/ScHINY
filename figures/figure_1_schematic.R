require(ggplot2)

FIGURE.1.FILE.NAME <- "Figure_1.png"

# Don't judge me
plot.figure.1.schematic <- function(shiny.data, gene, input) {
  fig1.img <- grid::rasterGrob(png::readPNG(file.path(IMAGE.PATH, FIGURE.1.FILE.NAME)), interpolate=TRUE)
  fig1.plt.hack <- qplot(1:10, 1:10, geom = "blank") +
    annotation_custom(fig1.img, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    theme_void()
  return(fig1.plt.hack)
}

exp.panel.figure.1.schematic <- function() {list()}

describe.figure.1.schematic <- function(...) {
"Schematic workflows for: (A) Growth of a transcriptionally-barcoded pool of 11 nitrogen metabolism transcription factor (TF) knockout strains and a wild-type control strain (B) Synthesis in microfluidic droplets of single-cell cDNA with a cell-specific index sequence (IDX) attached to the oligo-dT primer, and a common template switch oligo (TSO). cDNA is processed for whole-transcriptome libraries, to quantify gene expression. In parallel, PCR products are amplified containing the genotype-specific transcriptional barcode (BC) encoded on the KanR antibiotic resistance marker mRNA, to identify cell genotype. DNA libraries and PCR products are separately indexed for multiplexed sequencing (C) Processing of single-cell sequencing data into a Unique Molecular Identifier (UMI) count matrix which is used to learn a gene regulatory network with multi-task network inference from several different growth conditions."
}