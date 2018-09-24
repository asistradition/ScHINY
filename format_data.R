library(dplyr)
library(reshape2)

IN.FILES <- c('File1.tsv', 'File2.tsv')
IN.CONDITIONS <- c('C1', 'C2')
GENE.LIST <- read.table('genes.tsv')$V1
MELT.VARS <- c('Condition', 'Gene_Group', 'Num_Cells', 'Total_UMI', 'Mean_UMI', 'Genotype')

melted <- data.frame(Condition = character(), 
                     Gene_Group = factor(), 
                     Num_Cells = integer(), 
                     Total_UMI = integer(), 
                     Mean_UMI = double(),
                     Genotype = factor(), 
                     Gene = factor(), 
                     value = integer())

for (i in 1:length(IN.FILES)){
  data <- read.table(IN.FILES[i], header=TRUE)
  data$Total_UMI <- rowSums(data[GENE.LIST])
  data$Mean_UMI <- data$Total_UMI / data$Num_Cells
  data$Condition <- IN.CONDITIONS[i]
  
  data %>%
    melt(id.vars=MELT.VARS, variable.name='Gene') %>%
    rbind(melted)
}

melted$LibraryNormalized <- melted$value / melted$Total_UMI * 100
melted$Count <- melted$value / melted$Num_Cells
melted$Num_Cells <- NULL
melted$Total_UMI <- NULL
melted$Mean_UMI <- NULL
melted$Genotype <- NULL
melted$value <- NULL
write.table(melted, file="Output.tsv", sep="\t")