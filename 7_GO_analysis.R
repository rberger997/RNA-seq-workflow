#------------------------------------------------------------
#                RNA seq analysis workflow
#                     7_GO_analysis
#                   Ryan Berger 2-5-18
#
#   RNA seq workflow starting with kallisto aligned reads
#
#  Purpose: This is the prep step for exploratory analysis. 
#  Input: abundance.tsv files
#  Output: DESeq dataset (dds)
#
#
# Source: https://www.bioconductor.org/help/workflows/rnaseqGene/
#------------------------------------------------------------

#------------------------------------------------------------
#                     7_GO_analysis

# Need a dataframe with log2fold change and gene ID
geneList <- res.df[,c(2,8)]
head(geneList)

# Selection criteria: genes that are at least 2 fold changed
gene <- geneList[abs(geneList$log2FoldChange) > 1, ]
gene <- as.character(row.names(gene))
head(gene)
class(gene)

library(clusterProfiler)
gene.df <- bitr(gene, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gene <- as.character(gene.df$ENTREZID)

## GO over-representation test
# input:
# gene = character vector of selected genes for analysis (use ENTREZID format)
# universe = character vector of all genes in dataset

ont = 'BP'  # Bio process (BP), Mol function (MF), Cell compartment (CC)
ego <- enrichGO(gene = gene,
                universe = geneList$entrez,
                # keytype = 'ENSEMBL',
                OrgDb = org.Hs.eg.db,
                ont = ont,
                pAdjustMethod = 'BH',
                minGSSize = 10,
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable = TRUE)
## Common error:
## Only input the ENTREZID format, not working when you put in ENSEMBL (even if specified in keytype)


# Simplify GO terms in over-representation
ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)

# Visualize over representation 
dotplot(ego2, 
        showCategory = 20, 
        title = paste(sample1, '/', sample2, 'GO enrichment', '-',ont))

#------------------------------------------------------------

# Save GO plot as PDF
GOplot <- recordPlot()  # Run after plot is done
pdf(paste(outputdir,'/',sample1,'_over_',sample2,'_GOplot-',ont,'.pdf', sep = ''), width = 10, height = 6)
GOplot
dev.off()


# Optional plots:

enrichMap(ego2, n = 20) # n = number of top nodes to look at, default is 50 
## categorySize can be scaled by 'pvalue' or 'geneNum'


# Cnetplot
cnetplot(ego2, categorySize = 'pvalue')
cnetplot


# GO graph
plotGOgraph(ego)


# groupGo plot - doesn't seem to be working, plot is the same every time
# Input is a character vector of genes
ggo <- groupGO(gene = gene,
               OrgDb = org.Hs.eg.db,
               #keytype = 'ENSEMBL', # specify the format of gene terms in list
               ont = 'MF',  # MF (molecular function), BP (biol process), CC (cell compart)
               level = 3,
               readable = TRUE)
head(ggo)
?groupGO
# Plot results of GO
barplot(ggo, drop = TRUE, showCategory = 20, decreasing = TRUE)
# Problem: barplot is not in order of most enriched terms, always the same order
# decreasing = TRUE doesn't change it