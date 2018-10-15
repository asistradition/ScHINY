library(dplyr)
library(reshape2)

IN.FILES <- c('ScUTR_090818_Glutamine_BCDEL1_filtered_bulk.tsv', 'ScUTR_090818_Proline_BCDEL1_filtered_bulk.tsv', 'ScUTR_100118_AS_BCDEL1_filtered_bulk.tsv', 'ScUTR_100118_Urea_BCDEL1_filtered_bulk.tsv')
IN.CONDITIONS <- c('Glutamine', 'Proline', 'AmmoniumSulfate', 'Urea')
GENE.LIST <- read.table('genes.tsv')$V1
MELT.VARS <- c('Condition', 'Gene_Group', 'Num_Cells', 'Total_UMI', 'Mean_UMI', 'Genotype')
GENE.GROUPS <- list(c(rep('WT', 6), rep('dal80', 6), rep('dal81', 6), rep('dal82', 6), rep('gat1', 6), rep('gcn4', 5), rep('gln3', 6), rep('gzf3', 6), rep('rtg1', 6), rep('rtg3', 6), rep('stp1', 6), rep('stp2', 6)),
                    c(rep('WT', 6), rep('dal80', 6), rep('dal81', 6), rep('dal82', 6), rep('gat1', 6), rep('gcn4', 6), rep('gln3', 6), rep('gzf3', 6), rep('rtg1', 6), rep('rtg3', 6), rep('stp1', 6), rep('stp2', 6)),
                    c(rep('WT', 6), rep('dal80', 6), rep('dal81', 6), rep('dal82', 6), rep('gat1', 6), rep('gcn4', 5), rep('gln3', 6), rep('gzf3', 6), rep('rtg1', 6), rep('rtg3', 6), rep('stp1', 6), rep('stp2', 6)),
                    c(rep('WT', 6), rep('dal80', 6), rep('dal81', 6), rep('dal82', 6), rep('gat1', 6), rep('gcn4', 6), rep('gln3', 6), rep('gzf3', 6), rep('rtg1', 6), rep('rtg3', 6), rep('stp1', 6), rep('stp2', 6))
                    )
names(GENE.GROUPS) <- c('Glutamine', 'Proline', 'AmmoniumSulfate', 'Urea')

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
  
  data %>%
    subset(select=-c(Genotype)) %>%
    rowSums -> data$Total_UMI
  
  data$Mean_UMI <- data$Total_UMI / data$Num_Cells
  data$Condition <- as.factor(IN.CONDITIONS[i])
  data$Gene_Group <- as.factor(unlist(GENE.GROUPS[IN.CONDITIONS[i]]))
  
  data %>%
    melt(id.vars=MELT.VARS, variable.name='Gene') -> melt_data
  
  melted <- rbind(melted, melt_data)
}

melted$LibraryNormalized <- melted$value / melted$Total_UMI * 100
melted$Count <- melted$value / melted$Num_Cells
melted$Num_Cells <- NULL
melted$Total_UMI <- NULL
melted$Mean_UMI <- NULL
melted$Genotype <- NULL
melted$value <- NULL
write.table(melted, file="ScHINY_Data.tsv", sep="\t")
