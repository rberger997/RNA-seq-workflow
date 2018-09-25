# RNA-seq-workflow

This is a computational workflow for processing and analyzing RNA sequencing data. The input data are sequencing alignment output files (in this case .tsv or .h5 files from Kallisto). Following the workflow in sequential order will generate the following outputs:

Gene-level results:
 - Matrix containing normalized counts for all genes in the genome 
 - Differential expression dataframes (.csv output) comparing gene expression between samples
 - Principal component analysis (PCA) plot of all samples
 - Summary table of gene changes based on selection criteria
 - Clustered heatmap of differential expression results
 - Volcano plots of differential expression results
 - Interactive volcano plots
 
 Pathway-level results:
 - Gene set enrichment analysis (GSEA) using MSigDB Hallmark pathways
 - GSEA using Reactome pathways
 - Clustered heatmaps of global pathway changes
 
