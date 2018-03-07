#------------------------------------------------------------
#                RNA seq analysis workflow
#                3_differential_expression
#                   Ryan Berger 2-5-18
#
#   RNA seq workflow starting with kallisto aligned reads
#
#  Purpose: compare gene expression between two samples
#  Input: DESeq dataset (dds)
#  Output: results (res) and dataframe (res.df) 
#
#
# Source: https://www.bioconductor.org/help/workflows/rnaseqGene/
#------------------------------------------------------------


#------------------------------------------------------------
#                3_differential_expression

dds <- DESeq(dds)
head(dds)

# Build results table. Define what 2 conditions to compare using contrast vector
sampleTable

# Select what samples are being compared to each other
sample1 <- 'Styphimurium'
sample2 <- 'PBS'
contrast <- c('condition', sample1, sample2) # factor name, numerator condition, denominator condition
res <- results(dds, contrast = contrast)

res
summary(res)

## Add annotation - symbol and entrezID
library(AnnotationDbi)
library(org.Hs.eg.db)
columns(org.Hs.eg.db)

# Add column for gene symbol
res$symbol <- rownames(res)
# Add column for gene Entrez ID
res$entrez <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = 'ENTREZID',
                     keytype = 'SYMBOL',
                     multiVals = 'first')
# Add column for gene name
res$name <- mapIds(org.Hs.eg.db,
                   keys = rownames(res),
                   column = 'GENENAME',
                   keytype = 'SYMBOL',
                   multiVals = 'first')
# Make results dataframe
res.df <- as.data.frame(res)
res.df <- res.df[order(res.df$padj),]
head(res.df)

# Save output file
write.csv(res.df, file = paste(outputdir,'/',sample1,'_over_', sample2,'_diff_expression.csv', sep = ''))
# End module. Proceed to gene summary table, MA plot, heatmaps, etc.
#------------------------------------------------------------