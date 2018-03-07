#------------------------------------------------------------
#                RNA seq analysis workflow
#                    1_counts_to_dds
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
#                    Set up directories

# Where are the input (kallisto alignment) files located?
dir <- '~/Desktop/RNAseq stuff/Run1822/kallisto-Run_1822'


# Where should the output files go?
output_folderID <- '2.28.18'  # Name of output folder

dir.create(file.path('~/Desktop/RNAseq stuff/Results', output_folderID))  # Optional: Create new folder for output
outputdir <- file.path('~/Desktop/RNAseq stuff/Results', output_folderID) # Output folder
setwd(outputdir)  # Change working directory to output folder

#------------------------------------------------------------
#                    1_counts_to_dds

# Need to import the .tsv output files using tximport
# First set up a gene reference (tx2gene) for annotation from ENSEMBL transcipt IDs to gene IDs
# source("https://bioconductor.org/biocLite.R")
# biocLite("EnsDb.Hsapiens.v75")
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
edb

# Create a dataframe of transcripts from ensembl human genome
Tx <- transcripts(edb,
                  columns = c(listColumns(edb , "tx"), "gene_name"),
                  return.type = "DataFrame")
head(Tx)
# assign columns 1 (transcript ID) and 9 (gene name) to dataframe tx2gene
tx2gene <- Tx[,c(1,9)]
head(tx2gene)

# File import
library(dplyr)
# Load in sample key.
# This sample key is all the 8h time samples from run 1822
sample_key <- read.csv(file.path(dir,'sample_key1.csv')) # shoule be 32 obs of 21 variables

# Optional: filter out samples from the key that you don't wan to analyze
sample_key <- filter(sample_key, hr == 8) # Only include 8h samples
sample_key <- sample_key[c(1:4, 13:32), ] # Filter out PMN only, PMN+HIOs samples
sample_key <- filter(sample_key, code_name %in% c('PBS', 
                                                  'Styphimurium'
                                                  # 'Senteritidis'
                                                  # 'PBS+PMNs', 
                                                  # 'Styphimurium+PMNs', 
                                                  # 'SEnt+PMNs' 
                                                  # 'PMN', 
                                                  # 'HIOs+PMNs'
))

# Record of samples used in the tximport
txiSamples <- as.vector(unique(sample_key$code_name))
txiSamples

sample_key$file_name
# Set up path to read files into tximport
files <- file.path(dir, sample_key$file_name, 'abundance.tsv')
# Add sample IDs to files
names(files) <- sample_key$short_name
files

library(tximport)
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
tail(txi$counts, 100)

# Prep for DESeq2 object
library(DESeq2)

# Create a table for sample names. Define the factor (for HIOs it's what was injected) and set the sample names (in this case used from the sample_key.csv file)
sampleTable <- data.frame(condition = factor(sample_key$code_name))
rownames(sampleTable) <- colnames(txi$counts)

# Create DESeq dataset
dds <- DESeqDataSetFromTximport(txi, sampleTable, design = ~condition)
head(dds)
dds
# dds is now ready for DESeq() see DESeq2 vignette

# Filter out genes with zero reads
nrow(dds) # 35300 total genes
# Filter out rows with no counts
dds <- dds[rowSums(counts(dds)) > 1, ]
nrow(dds) 

# Account for transcript length
dds <- DESeq2::estimateSizeFactors(dds)
ddscounts <- DESeq2::counts(dds, normalized = TRUE)
# Save output file:

write.csv(ddscounts, file = file.path(outputdir, "complete-dataset_DESeq2-normalized-counts.csv"))

## RNA seq reads are now in a DESeq dataset. End module, proceed to exploratory analysis.
#------------------------------------------------------------