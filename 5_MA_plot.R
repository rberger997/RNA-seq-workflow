#------------------------------------------------------------
#                RNA seq analysis workflow
#                        5_MA_plot
#                   Ryan Berger 2-5-18
#
#   RNA seq workflow starting with kallisto aligned reads
#
#  Purpose: Make MA plot from differential expression data
#  Input: Differential expression results object (res)
#  Output: MA plot showing fold change and counts
#
#
# Source: https://www.bioconductor.org/help/workflows/rnaseqGene/
#------------------------------------------------------------

#------------------------------------------------------------
#                        5_MA_plot

plotMA(res, 
       ylim = c(-7,7), 
       main = paste(sample1, '/', sample2,'\n', 'MA plot'))

#------------------------------------------------------------

# Save MA plot as PDF
MAplot <- recordPlot()  # Run after MA plot is done
pdf(paste(outputdir,'/',sample1,'_over_',sample2,'_MAplot.pdf', sep = ''), width = 6, height = 6)
MAplot
dev.off()

# Save MA plot as png
png(filename = file.path(outputdir, paste(sample1, '_over_',sample2,'_MAplot.png', sep = '')),
    height = 1.5, width = 6, units = 'in', res = 500)
MAplot
dev.off()

## Labeling options:


# Label top significant gene on MA plot
topGene <- rownames(res)[which.min(res$padj)]
topGenep <- resOrdered[1,7]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=1, lwd=2, pch = 1)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})


# Add multiple labels to plot
# This adds labels to first 10 rows of ordered results (min p values)
resOrdered <- res[order(res$pvalue),]
num <- 10  # Top n number of genes to label on plot
head(resOrdered, num)
with(resOrdered[1:num,], {
   points(baseMean, log2FoldChange, col="dodgerblue", cex=1, lwd=2, pch = 1) # Circle around points
  text(baseMean, log2FoldChange, resOrdered[1:num, 7], pos=4, col="dodgerblue", cex = .8) # Text labels
}) # pos = 1 (under), 2 (left), 3 (above), 4 (right)


# Label a specific gene(s) of interest
library(dplyr)
gene <- filter(res.df, symbol == 'ST8SIA4') #'MT1G' 'MT2A'
# Filter one value: symbol == 'X', filter a list: symbol %in% c(a,b,c)
with(gene, {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=1, lwd=2, pch = 1)
  text(baseMean, log2FoldChange, gene$symbol, pos=2, col="dodgerblue")
})
