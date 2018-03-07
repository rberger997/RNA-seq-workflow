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
  x <- x[!is.na(x$symbol),] # Remove transcripts with NA gene associated
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

# Save PDF of table
pdf(paste(outputdir,'/',sample1,'_over_',sample2,'_GeneSummary.pdf', sep = ''), width = 6, height = 1.5)
library(gridExtra)
grid.arrange(tableGrob(GeneSummary), 
             top = paste(sample1, ' over ', sample2, ' gene summary', sep = ''))
dev.off()

# Save png of table
png(filename = file.path(outputdir, paste(sample1, '_over_',sample2,'_GeneSummary.png', sep = '')),
    height = 1.5, width = 6, units = 'in', res = 500)
grid.arrange(tableGrob(GeneSummary), 
             top = paste(sample1, ' over ', sample2, ' gene summary', sep = ''))
dev.off()

#------------------------------------------------------------

# Save significant gene list as PDF
SigGeneList <- function(df, p, fc){
  x <- df
  sig.df <- filter(x, padj < p & abs(log2FoldChange) > fc)
  sig.df <- sig.df[,c(1:2,6,7,9)]
  head(sig.df)
  sig.df <- sig.df[order(sig.df$padj),]
  sig.df <- sig.df[!is.na(sig.df$symbol),]
  rownames(sig.df) <- NULL
  sig.df <- sig.df
}


# Get significant gene list - res.df, p cutoff, log 2 FC cutoff
sig.genes <- SigGeneList(res.df, 0.05, 1)

# Save as PDF
pdf(paste(outputdir,'/',sample1,'_over_',sample2,'_significant.pdf', sep = ''), width = 12, height = 12)
library(gridExtra)
grid.arrange(tableGrob(sig.genes),
             top = paste(sample1, ' over ', sample2, ' significant genes', sep = ''))
dev.off()
