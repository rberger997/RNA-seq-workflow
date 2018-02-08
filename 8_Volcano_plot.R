#------------------------------------------------------------
#                RNA seq analysis workflow
#                     8_Volcano_plot
#                   Ryan Berger 2-5-18
#
#   RNA seq workflow starting with kallisto aligned reads
#
#  Purpose: Make volcano plot of fold change and p values  
#  Input: results object (res) from differential expression
#  Output: Volcano plot
#
#
# Source: https://www.bioconductor.org/help/workflows/rnaseqGene/
#------------------------------------------------------------

#------------------------------------------------------------
#                     8_Volcano_plot

head(res)
# Data is in a dataframe with Gene name, log2foldchange, pvalue, padj

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(padj), 
               pch=20, 
               main = paste(sample1, '/', sample2, '\n', 'Volcano Plot'),
               # xlim=c(-5,8),
               col = 'darkgray'
))

# Add colored points: red if padj<0.05, orange if log2FC>1, green if both
with(subset(res, padj<0.05), 
     points(log2FoldChange, -log10(padj), 
            pch=20, 
            col='red'))
with(subset(res, abs(log2FoldChange)>1), 
     points(log2FoldChange, -log10(padj), 
            pch=20, 
            col='green'))
with(subset(res, padj<0.05 & abs(log2FoldChange)>1), 
     points(log2FoldChange, -log10(padj), 
            pch=20, 
            col='orange'))

# Label points with the textxy function from the calibrate plot
# install.packages('calibrate')
library(calibrate)
with(subset(res, -log10(padj) > 20 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(padj), labs=symbol, cex=.5) )

#------------------------------------------------------------

##  Optional: save plot as PDF
volcano <- recordPlot() # Run after plot is finished
pdf(paste(outputdir,'/',sample1,'_over_',sample2,'_volcanoplot.pdf', sep = ''), width = 6, height = 6)
volcano
dev.off()
