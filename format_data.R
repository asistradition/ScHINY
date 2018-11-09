library(dplyr)
library(reshape2)

# Bulk data pileup

IN.FILES <- c('ScUTR_090818_Glutamine_BCDEL1_filtered_bulk.tsv', 'ScUTR_090818_Proline_BCDEL1_filtered_bulk.tsv', 'ScUTR_100118_AS_BCDEL1_filtered_bulk.tsv', 'ScUTR_100118_Urea_BCDEL1_filtered_bulk.tsv', 'ScUTR_101118_YPD_BCDEL1_filtered_bulk.tsv', 'ScUTR_101118_YPD_Rapa_BCDEL1_filtered_bulk.tsv', 'ScUTR_101118_CStarve_filtered_bulk.tsv', 'ScUTR_102718_MMGlucose_filtered_bulk.tsv', 'ScUTR_102718_MMEtOH_filtered_bulk.tsv', 'ScUTR_102718_YPD_Diauxic_filtered_bulk.tsv', 'ScUTR_102718_YPEtOH_filtered_bulk.tsv')
IN.CONDITIONS <- c('Glutamine', 'Proline', 'AmmoniumSulfate', 'Urea', 'YPD', 'YPDRapa', 'CStarve', 'MinimalGlucose', 'MinimalEtOH', 'YPDDiauxic', 'YPEtOH')
GENE.LIST <- read.table('genes.tsv')$V1
GENE.GROUPS <- list(c(rep('WT', 6), rep('dal80', 6), rep('dal81', 6), rep('dal82', 6), rep('gat1', 6), rep('gcn4', 5), rep('gln3', 6), rep('gzf3', 6), rep('rtg1', 6), rep('rtg3', 6), rep('stp1', 6), rep('stp2', 6)),
                    c(rep('WT', 6), rep('dal80', 6), rep('dal81', 6), rep('dal82', 6), rep('gat1', 6), rep('gcn4', 6), rep('gln3', 6), rep('gzf3', 6), rep('rtg1', 6), rep('rtg3', 6), rep('stp1', 6), rep('stp2', 6)),
                    c(rep('WT', 6), rep('dal80', 6), rep('dal81', 6), rep('dal82', 6), rep('gat1', 6), rep('gcn4', 5), rep('gln3', 6), rep('gzf3', 6), rep('rtg1', 6), rep('rtg3', 6), rep('stp1', 6), rep('stp2', 6)),
                    c(rep('WT', 6), rep('dal80', 6), rep('dal81', 6), rep('dal82', 6), rep('gat1', 6), rep('gcn4', 6), rep('gln3', 6), rep('gzf3', 6), rep('rtg1', 6), rep('rtg3', 6), rep('stp1', 6), rep('stp2', 6)),
                    c(rep('WT', 6), rep('dal80', 6), rep('dal81', 6), rep('dal82', 6), rep('gat1', 6), rep('gcn4', 6), rep('gln3', 6), rep('gzf3', 6), rep('rtg1', 6), rep('rtg3', 6), rep('stp1', 6), rep('stp2', 4)),
                    c(rep('WT', 6), rep('dal80', 6), rep('dal81', 6), rep('dal82', 6), rep('gat1', 6), rep('gcn4', 6), rep('gln3', 6), rep('gzf3', 6), rep('rtg1', 6), rep('rtg3', 6), rep('stp1', 6), rep('stp2', 4)),
                    c(rep('WT', 6), rep('dal80', 6), rep('dal81', 6), rep('dal82', 6), rep('gat1', 6), rep('gcn4', 6), rep('gln3', 6), rep('gzf3', 6), rep('rtg1', 6), rep('rtg3', 6), rep('stp1', 6), rep('stp2', 4)),
                    c(rep('WT', 6), rep('dal80', 6), rep('dal81', 6), rep('dal82', 6), rep('gat1', 6), rep('gcn4', 6), rep('gln3', 6), rep('gzf3', 6), rep('rtg1', 6), rep('rtg3', 6), rep('stp1', 6), rep('stp2', 6)),
                    c(rep('WT', 6), rep('dal80', 6), rep('dal81', 6), rep('dal82', 6), rep('gat1', 6), rep('gcn4', 6), rep('gln3', 6), rep('gzf3', 6), rep('rtg1', 6), rep('rtg3', 6), rep('stp1', 6), rep('stp2', 6)),
                    c(rep('WT', 6), rep('dal80', 6), rep('dal81', 6), rep('dal82', 6), rep('gat1', 6), rep('gcn4', 6), rep('gln3', 6), rep('gzf3', 6), rep('rtg1', 6), rep('rtg3', 6), rep('stp1', 6), rep('stp2', 6)),
                    c(rep('WT', 6), rep('dal80', 6), rep('dal81', 6), rep('dal82', 6), rep('gat1', 6), rep('gcn4', 6), rep('gln3', 6), rep('gzf3', 6), rep('rtg1', 6), rep('rtg3', 6), rep('stp1', 6), rep('stp2', 6)))
names(GENE.GROUPS) <- IN.CONDITIONS

headers <- read.table(IN.FILES[1], nrows = 1, header = FALSE, sep ='\t', stringsAsFactors = FALSE)
headers[,ncol(headers)+1] = "Condition"
headers[,ncol(headers)+1] = "Mean_UMI"
headers[,ncol(headers)+1] = "Total_UMI"
headers[,ncol(headers)+1] = "Gene_Group"
all_data <- data.frame(matrix(ncol = ncol(headers), nrow = 0))
colnames(all_data) <- as.vector(headers)

for (i in 1:length(IN.FILES)){
  data <- read.table(IN.FILES[i], header=TRUE)
  
  data %>%
    subset(select=-c(Genotype, Num_Cells)) %>%
    rowSums -> data$Total_UMI
  
  data$Mean_UMI <- data$Total_UMI / data$Num_Cells
  data$Condition <- as.factor(IN.CONDITIONS[i])
  data$Gene_Group <- as.factor(unlist(GENE.GROUPS[IN.CONDITIONS[i]]))
  
  all_data <- rbind(all_data, data)
}

write.table(all_data, file="ScHINY_Data.tsv", sep="\t")

# Individual cell pileup

IN.FILES <- c('ScUTR_090818_Glutamine_BCDEL1_filtered.tsv', 'ScUTR_090818_Proline_BCDEL1_filtered.tsv', 'ScUTR_100118_AS_BCDEL1_filtered.tsv', 'ScUTR_100118_Urea_BCDEL1_filtered.tsv', 'ScUTR_101118_YPD_BCDEL1_filtered.tsv', 'ScUTR_101118_YPD_Rapa_BCDEL1_filtered.tsv', 'ScUTR_101118_CStarve_filtered.tsv', 'ScUTR_102718_MMGlucose_filtered.tsv', 'ScUTR_102718_MMEtOH_filtered.tsv', 'ScUTR_102718_YPD_Diauxic_filtered.tsv', 'ScUTR_102718_YPEtOH_filtered.tsv')
IN.CONDITIONS <- c('Glutamine', 'Proline', 'AmmoniumSulfate', 'Urea', 'YPD', 'YPDRapa', 'CStarve', 'MinimalGlucose', 'MinimalEtOH', 'YPDDiauxic', 'YPEtOH')
GENE.LIST <- read.table('genes.tsv')$V1

headers <- read.table(IN.FILES[1], nrows = 1, header = FALSE, sep ='\t', stringsAsFactors = FALSE)
headers[,ncol(headers)+1] = "Condition"
headers[,ncol(headers)+1] = "tenXBarcode"
all_data <- data.frame(matrix(ncol = ncol(headers), nrow = 0))
colnames(all_data) <- as.vector(headers)

for (i in 1:length(IN.FILES)){
  data <- read.table(IN.FILES[i], header=TRUE)
  data$Condition = IN.CONDITIONS[i]
  data$tenXBarcode <- rownames(data)
  rownames(data) <- paste(1:nrow(data), "_", IN.CONDITIONS[i], sep="")
  all_data <- rbind(all_data, data)
  data = NULL
}

write.table(all_data, file="103118_SS_Data.tsv", sep="\t")
