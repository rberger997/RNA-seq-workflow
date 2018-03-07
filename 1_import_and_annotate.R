#------------------------------------------------------------
#                RNA seq analysis workflow
#                  1_import_and_annotate
#                   Ryan Berger 2-5-18
#
#   RNA seq workflow starting with kallisto aligned reads
#
#  Purpose: This is the prep step for exploratory analysis. 
#  Input: abundance.tsv files
#  Output: DESeq dataset (dds), results dataframe (res.df)
#
#
# Source: https://www.bioconductor.org/help/workflows/rnaseqGene/
#------------------------------------------------------------


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
# assign columns 1 (transcript ID) and 7 (gene ID) to dataframe tx2gene
tx2gene <- Tx[,c(1,7)]
head(tx2gene)

# 
# Assign directory to RNA seq folder with alignment files
dir <- '~/Desktop/RNAseq stuff/Run1822/kallisto-Run_1822'

# Load in sample key.
# This sample key is all the 8h time samples from run 1822
setwd('~/Desktop/RNAseq stuff/Run1822/kallisto-Run_1822')
sample_key <- read.csv('sample_key.csv') # shoule be 32 obs of 21 variables
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