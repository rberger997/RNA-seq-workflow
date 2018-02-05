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

sampleDistMatrix <- as.matrix(sampleDists)
head(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# PCA plot
plotPCA(rld, intgroup = 'condition')
#------------------------------------------------------------