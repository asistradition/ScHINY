library(dplyr)
library(reshape2)

# Bulk data pileup
if(!exists('WORKING.PATH')) {WORKING.PATH <- "."}
if(!exists('COUNTS.FILE')) {COUNTS.FILE <- "103118_SS_Data.tsv.gz"}
if(!exists('EXPRESSION.SUMMARY.OUT')) {EXPRESSION.SUMMARY.OUT <- "Expression_Summary"}
if(!exists('FIG2.OUT')) {FIG2.OUT <- "Figure2"}
if(!exists('FIG4.OUT')) {FIG4.OUT <- "Expression_Count"}


setwd(WORKING.PATH)
if (!exists('count_data')) {
  count_data <- read.table(gzfile(COUNTS.FILE), header=TRUE, sep="\t")
  rownames(count_data) <- gsub("-", ".", rownames(count_data))
}

#
cond_order <- c("YPD", "YPDDiauxic", "YPDRapa", "YPEtOH", "MinimalGlucose", "MinimalEtOH", "Glutamine", "Proline", "AmmoniumSulfate", "Urea", "CStarve")

# Prepare data for expression summary
OUT.PATH <- EXPRESSION.SUMMARY
dir.create(OUT.PATH, showWarnings = FALSE)
count_data %>%
  select(-c('Replicate', 'tenXBarcode')) %>%
  group_by(Condition, Genotype_Group, Genotype) %>%
  summarise_all(sum) %>% 
  as.data.frame() -> gene_count

count_data %>%
  select(-c('Replicate', 'tenXBarcode')) %>%
  group_by(Condition, Genotype_Group, Genotype) %>%
  summarize(Num_Cells = n()) %>%
  as.data.frame() -> cell_count

gene_count %>%
  mutate(TotalUMI = rowSums(select(., -c('Condition', 'Genotype_Group', 'Genotype')))) %>%
  left_join(cell_count, by = c("Condition", "Genotype_Group", "Genotype")) %>%
  mutate(MeanUMI = TotalUMI / Num_Cells) %>%
  arrange(match(Condition, cond_order)) -> gene_count

rownames(gene_count) <- paste(gene_count$Condition, gene_count$Genotype, sep="_")

write.table(select(gene_count, c("Condition", "Genotype_Group", "Genotype", "TotalUMI", "MeanUMI", "Num_Cells")), file=file.path(OUT.PATH, "meta_data.tsv"), sep="\t")
gene_count <- select(gene_count, -c("Condition", "Genotype_Group", "Genotype", "TotalUMI", "MeanUMI", "Num_Cells"))

for (i in 1:ncol(gene_count)) {
  write.table(gene_count[,i], file=file.path(OUT.PATH, paste(colnames(gene_count)[i], ".tsv", sep = "")), sep="\t", row.names = FALSE, col.names = FALSE)
}

scaler <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

# Prepare data for figure 2
OUT.PATH <- FIG2.OUT
dir.create(OUT.PATH, showWarnings = FALSE)
count_data %>%
  select(-c('Condition', 'Genotype_Group', 'Genotype', 'Replicate', 'tenXBarcode')) %>%
  mutate_all(funs(scaler(.))) %>% 
  as.data.frame() -> gene_count

for (i in 1:ncol(gene_count)) {
  write.table(gene_count[,i], file=file.path(OUT.PATH, paste(colnames(gene_count)[i], ".tsv", sep = "")), sep="\t", row.names = FALSE, col.names = FALSE)
}

# Prepare data for figure 4
OUT.PATH <- FIG4.OUT
dir.create(OUT.PATH, showWarnings = FALSE)
count_data %>%
  select(-c('Replicate', 'tenXBarcode')) %>%
  mutate(TotalUMI = rowSums(select(., -c('Condition', 'Genotype_Group', 'Genotype')))) %>%
  as.data.frame() -> gene_count

write.table(select(gene_count, c("Condition", "Genotype_Group", "Genotype", "TotalUMI")), file=file.path(OUT.PATH, "meta_data.tsv"), sep="\t")
gene_count <- select(gene_count, -c("Condition", "Genotype_Group", "Genotype", "TotalUMI"))
for (i in 1:ncol(gene_count)) {
  write.table(gene_count[,i], file=file.path(OUT.PATH, paste(colnames(gene_count)[i], ".tsv", sep = "")), sep="\t", row.names = FALSE, col.names = FALSE)
}

