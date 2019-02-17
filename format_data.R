DATA.TYPES <- c('counts', 'logcounts', 'activity')

if (!exists('count_data')) {source('figure_common_functions.R')}
if (!exists('sce.data')) {sce.data <<- dimred_for_umap(load_data_into_sce(count.data = dplyr::select(count_data, -META.DATA.COLUMNS), 
                                                                          meta.data = meta_data), 
                                                       random_state=RANDOM.SEED)}

OUT.DIR <- "shiny_data"
dir.create(OUT.DIR, showWarnings = FALSE)

# Read in the learned network
network.learned <- read.table(file.path(WORKING.PATH, NETWORK.PATH, NETWORK.FILE), sep = "\t", header = TRUE)
rownames(network.learned) <- network.learned$target
network.learned$target <- NULL

# Calculate transcription factor activity from a network and an expression matrix
tfa <- function(prior, exp.mat, noself=TRUE, dup.self=TRUE) {
  tfwt <- apply(prior != 0, 2, sum) > 0
  
  duplicates <- c()
  if (dup.self) {
    duplicates <- duplicated(prior[, tfwt], MARGIN=2) |
      duplicated(prior[, tfwt], MARGIN=2, fromLast = TRUE)
    duplicates <- colnames(prior)[tfwt][duplicates]
  }
  
  tfs <- setdiff(colnames(prior), duplicates)
  tfs <- intersect(tfs, rownames(prior))
  if (noself) {
    diag(prior[tfs, tfs]) <- 0
  }
  
  activities <- matrix(0, ncol(prior), ncol(exp.mat), 
                       dimnames=list(colnames(prior), colnames(exp.mat)))
  
  if (any(tfwt)) {
    activities[tfwt, ] <- corpcor::pseudoinverse(prior[, tfwt, drop=FALSE]) %*% exp.mat
  }
  
  use.exp <- intersect(colnames(prior)[!tfwt], rownames(exp.mat))
  activities[use.exp, ] <- exp.mat[use.exp, ]
  
  return(activities)
}

count.for.activity <- as.matrix(t(dplyr::select(count_data, -META.DATA.COLUMNS)[, rownames(network.learned)]))
activities <- tfa(as.matrix(network.learned), count.for.activity)
SingleCellExperiment::reducedDim(sce.data, "activity") <- t(activities)

sce.obj.fig3 <- load_data_into_sce(count.data = dplyr::select(count_data, -META.DATA.COLUMNS), meta.data = meta_data)
sce.conditions <- list()
for (condition in rev(CONDITION.LEVELS)) {sce.conditions[[condition]] <- sce.obj.fig3[, sce.obj.fig3$Condition == condition]}
rm(sce.obj.fig3)


print("Normalizing Figure 3 Data")
sce.conditions <- lapply(sce.conditions, normalizeSCE)
print("PCA for Figure 3 Data")
sce.conditions <- lapply(sce.conditions, function(x) {runPCA(x, ncomponents=50, ntop=Inf, method='irlba', exprs_values='logcounts')})
print("UMAP for Figure 3 Data")
sce.conditions <- lapply(sce.conditions, function(x) {runUMAP(x, random_state=RANDOM.SEED, min_dist=0.1, use_dimred='PCA')})
print("SNN for Figure 3 Data")
sce.conditions <- lapply(sce.conditions, apply.snn.cluster)

fig3.data <- NULL
for (condition in unique(count_data$Condition)) {
  df.umap <- as.data.frame(reducedDim(sce.conditions[[condition]], 'UMAP'))
  colnames(df.umap) <- c("UMAP.FIG3.1", "UMAP.FIG3.2")
  df.umap$Cluster <- sce.conditions[[condition]]$cluster
  df.umap$Condition <- condition
  fig3.data <- rbind(fig3.data, df.umap)
}

# Write the meta data with UMAP coords
df <- as.data.frame(reducedDim(sce.data, 'UMAP'))
colnames(df) <- c('UMAP1', 'UMAP2')
df <- cbind(meta_data[, c("Condition", "Genotype_Group", "Genotype")], df)
df$UMI <- colSums(SummarizedExperiment::assay(sce.data, 'counts'))

print(all(df$Condition == fig3.data$Condition))
fig3.data$Condition <- NULL
df <- cbind(df, fig3.data)

write.table(df, gzfile(file.path(OUT.DIR, "jackson_2019_meta_data.tsv.gz"), open = "w"), quote = FALSE, sep = "\t", row.names = FALSE)

gene.names <- colnames(count_data)[!colnames(count_data) %in% META.DATA.COLUMNS]

for (gene.name in gene.names) {
  gene.data <- as.data.frame(t(SummarizedExperiment::assay(sce.data[gene.name, ], 'counts')))
  gene.data <- cbind(gene.data, as.data.frame(t(SummarizedExperiment::assay(sce.data[gene.name, ], 'logcounts'))))
  if (!gene.name %in% colnames(SingleCellExperiment::reducedDim(sce.data, 'activity'))) {gene.data$activity <- NA}
  else {gene.data <- cbind(gene.data, as.data.frame(SingleCellExperiment::reducedDim(sce.data, 'activity')[, gene.name]))}
  colnames(gene.data) <- c("counts", "logcounts", "activity")
  file.connection <- gzfile(file.path(OUT.DIR, paste0(gsub("-", "\\.", gene.name), ".tsv.gz")), open = "w")
  write.table(gene.data, file.connection, quote = FALSE, sep = "\t", row.names = FALSE)
  close(file.connection)
}
