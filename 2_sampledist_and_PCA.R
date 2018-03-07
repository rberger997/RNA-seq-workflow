#------------------------------------------------------------
#                RNA seq analysis workflow
#                    2_sampledist_and_PCA
#                   Ryan Berger 2-5-18
#
#   RNA seq workflow starting with kallisto aligned reads
#
#  Purpose: Assess global sample differences with sampledist/PCA plots 
#  Input: DESeq dataset (dds)
#  Output: sampledistance plot, PCA plot
#
#
# Source: https://www.bioconductor.org/help/workflows/rnaseqGene/
#------------------------------------------------------------


#------------------------------------------------------------
#                    2_sampledist_and_PCA

# Perform variance stabilizing transformations on dds using rlog
# use argument blind = FALSE when multiple replicates are present
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

# Assess overall similarity between samples using Sample Distances
sampleDists <- dist(t(assay(rld)))
sampleDists
# visualize distances
library(pheatmap)
library(RColorBrewer)
# Make sample distance plot
sampleDistMatrix <- as.matrix(sampleDists)
head(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         treeheight_row = 15,
         treeheight_col = 15)
sampDistPlot <- recordPlot()

# PCA plot
plotPCA(rld, intgroup = 'condition')
PCAplot <- recordPlot()

#------------------------------------------------------------
# Optional:

# Save plots
name <- 'PBS_STM'   #File name - state samples in the plot
save_plot <- PCAplot  # set to the plot object to be saved
plot_type <- 'PCAplot'

# Save plot as PDF
pdf(paste(outputdir,'/',name,'_',plot_type,'.pdf', sep = ''), 
    width = 4, height = 4, onefile = FALSE)
save_plot
dev.off()  # End PDF

# Save plot as png
png(filename = file.path(outputdir, paste(name, '_',plot_type,'.png', sep = '')),
    height = 3, width = 4, units = 'in', res = 500)
save_plot
dev.off() # End png
#------------------------------------------------------------
