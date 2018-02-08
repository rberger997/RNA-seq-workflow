#------------------------------------------------------------
#                RNA seq analysis workflow
#                    6_cluster_heatmap
#                   Ryan Berger 2-5-18
#
#   RNA seq workflow starting with kallisto aligned reads
#
#  Purpose: Make cluster heatmap from differential expression data
#  Input: Differential expression results object (res), rlog transformed (rld)
#  Output: heatmap showing fold change across samples
#
#
# Source: https://www.bioconductor.org/help/workflows/rnaseqGene/
#------------------------------------------------------------

#------------------------------------------------------------
#                    6_cluster_heatmap

# Heatmap of top variance genes
library(genefilter)

# rld <- rlog(dds, blind = FALSE) # rlog transformation of dds 

# Cluster heatmap of all samples in dds
num1 <- 25  # number of genes in heatmap. (e.g. top n variance genes)
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), num1) # Define top variance genes
mat <- assay(rld)[topVarGenes, ]  # Matrix of the top variance genes from the set
mat <- mat - rowMeans(mat)  # expression is relative to mean of all samples. (fold over mean of all samples)
anno <- as.data.frame(colData(rld))  # sample IDs for annotation
topVarGeneIDs <- res[rownames(mat), 7]  # Gene labels as symbols
# Heatmap
pheatmap(mat, annotation_col = anno, labels_row = topVarGeneIDs)

#------------------------------------------------------------
## Heatmap options


# Save heatmap as PDF
hmap <- recordPlot() # Run after plot is finished
# Pick title format
GeneSetTitle <- 'MT genes' # Title of subset for heatmap
GeneSetTitle <- paste('Top',num1,'VarianceGenes', sep = '') # Title of subset for heatmap
pdf(paste(outputdir,'/',sample1,'_over_',sample2,'-',GeneSetTitle,'-heatmap.pdf',sep = ''), width = 5, height = 6, onefile = FALSE)
hmap
dev.off() # End PDF


## Make heatmap from a subset of samples in the dataset
matx <- assay(rld)
colnames(matx) # Look at samples
maty <- matx[, c(1:12)] # Select the samples you want for the heatmap
colnames(maty)
num1 <- 50 # number of genes in heatmap. (e.g. top n variance genes)
topVarGenes <- head(order(rowVars(maty), decreasing = TRUE), num1) # Define top variance genes
mat <- maty[topVarGenes,] # Select the top variance genes from the set
mat <- mat - rowMeans(mat) # expression is relative to mean of all samples
anno <- as.data.frame(colData(rld)) # sample IDs for annotation
topVarGeneIDs <- res[rownames(mat), 7] # Gene labels as symbols
# Heatmap of the 20 genes with the highest variance between samples
pheatmap(mat, annotation_col = anno, labels_row = topVarGeneIDs)




## Make a cluster heatmap of a custom list of genes
gene_list <- c('MT1A', 'MT1E', 'MT1F', 'MT1G', 'MT1H', 'MT1X', 'MT2A')
geneset <- row.names(subset(res, symbol %in% gene_list))

# Pick one of these two:
# mat <- assay(rld)[geneset, ]  # Look at gene set in all samples
# mat <- maty[geneset, ] # Look at gene set in subset of samples

mat <- mat - rowMeans(mat) # expression is relative to mean of all samples
anno <- as.data.frame(colData(rld)) # sample IDs for annotation
topVarGeneIDs <- res[rownames(mat), 7] # Gene labels as symbols
# Heatmap of the 20 genes with the highest variance between samples
pheatmap(mat, annotation_col = anno, labels_row = topVarGeneIDs)



## Make a cluster heatmap of top 25 genes changed 2FC+ and pvalue < 0.05 in SE vs PBS
sig.fc <- filter(res.df, abs(log2FoldChange) > 1 & padj < 0.05)
a <- sig.fc[order(sig.fc$padj), ]
a <- head(a, 25)
gene_list <- a$symbol
# Remove NAs from vector
gene_list <- gene_list[!is.na(gene_list)]
# Remove duplicates from vector
gene_list <- unique(gene_list)
x <- row.names(subset(res, symbol %in% gene_list))
mat1 <- assay(rld)[x, ]
mat1 <- mat1 - rowMeans(mat1)
IDs <- res[rownames(mat1), 7]
pheatmap(mat1, annotation_col = anno, labels_row = IDs)
head(assay(rld))