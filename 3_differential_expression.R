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
sample1 <- 'Senteritidis'
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
res$symbol <- mapIds(org.Hs.eg.db,
                     keys = row.names(res),
                     column = 'SYMBOL',
                     keytype = 'ENSEMBL',
                     multiVals = 'first')
# Add column for gene Entrez ID
res$entrez <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = 'ENTREZID',
                     keytype = 'ENSEMBL',
                     multiVals = 'first')
# Add column for gene name
res$name <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = 'GENENAME',
                     keytype = 'ENSEMBL',
                     multiVals = 'first')
# Make results dataframe
res.df <- as.data.frame(res)

# End module. Proceed to gene summary table, MA plot, heatmaps, etc.
#------------------------------------------------------------