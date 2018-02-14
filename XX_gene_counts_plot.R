#------------------------------------------------------------
#                RNA seq analysis workflow
#                    XX_gene_counts_plot
#                   Ryan Berger 2-5-18
#
#   RNA seq workflow starting with kallisto aligned reads
#
#  Purpose: Make plots of counts for a single gene or set of genes.
#  Input: Differential expression results (res)
#  Output: Plot of normalized counts in all samples
#
#
# Source: https://www.bioconductor.org/help/workflows/rnaseqGene/
#------------------------------------------------------------

#------------------------------------------------------------

# Plot single gene counts
gene_sym <- 'ST8SIA4'
gene_ens <- row.names(subset(res, symbol %in% gene_list))
plotCounts(dds, gene = gene_ens, intgroup = 'condition', main = res[gene_ens, 7], pch = 16, col = 'blue')

#------------------------------------------------------------

# Make multiple gene plots by looping and save plots as PDFs
save_dir <- '~/Desktop/RNAseq stuff/Results/MTgenes'  # Directory for saved plots to go
gene_list <- c('MT1A', 'MT1E', 'MT1F', 'MT1G', 'MT1H', 'MT1X', 'MT2A', 'GAPDH', 'ACTB')
geneset <- row.names(subset(res, symbol %in% gene_list))
# Loop for making plots
xlabel = 'HIOs Injection'
for (x in geneset){
  pdf(paste(save_dir,'/',res[x,7],'_counts.pdf',sep = ''), width = 6, height = 5, onefile = FALSE)
  plotCounts(dds, gene = x, intgroup = 'condition', main = res[x, 7], pch = 16, col = 'blue', xlab = xlabel)
  dev.off() # End PDF
}

#------------------------------------------------------------