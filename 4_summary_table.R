#------------------------------------------------------------
#                RNA seq analysis workflow
#                    4_summary_table
#                   Ryan Berger 2-5-18
#
#   RNA seq workflow starting with kallisto aligned reads
#
#  Purpose: Make summary of gene differential expression
#  Input: Differential expression results dataframe (res.df)
#  Output: Dataframe summary of gene changes
#
#
# Source: https://www.bioconductor.org/help/workflows/rnaseqGene/
#------------------------------------------------------------

#------------------------------------------------------------
#                    4_summary_table

# Function to make RNA seq summary table
# Input is: dataframe, p value cutoff, log2fold change cutoff
GeneTable <- function(x, p, fc){
  x <- x[order(x$padj),] # sort by lowest adj p value
  rowIDs <- c('Total genes with counts',
              paste('Adjusted p value <', p),
              paste('Adj. p <', p, '& log2 Fold Change >', fc))
  Total <- c(nrow(x), # Total column
             nrow(filter(x, padj < p)), 
             nrow(filter(x, padj < p & abs(log2FoldChange) > fc)))
  Increasing <- c(nrow(filter(x, log2FoldChange > 0)), # Increasing column
                  nrow(filter(x, padj < p & log2FoldChange > 0)), 
                  nrow(filter(x, padj < p & log2FoldChange > fc)))
  Decreasing <- c(nrow(filter(x, log2FoldChange < 0)), # Decreasing column
                  nrow(filter(x, padj < p & log2FoldChange < 0)), 
                  nrow(filter(x, padj < p & log2FoldChange < -fc)))
  # Make table
  subsetTable <- data_frame(Total, Increasing, Decreasing)
  subsetTable <- as.data.frame(subsetTable)
  rownames(subsetTable) <- rowIDs
  subsetTable
}

GeneSummary <- GeneTable(res.df, 0.05, 1)
GeneSummary
#------------------------------------------------------------