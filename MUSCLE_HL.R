#!/usr/bin/env Rscript --vanilla

#Running MUSCLE on sequences
#Installing required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("muscle")
library(BiocManager)
install("genbankr")
install("Biostrings")
install("ggtree")
install("annotate")
install("muscle") 
install("ape")
library(ape)
library(Biostrings)
library(annotate)
library(muscle)
install.packages("bios2mds")
library(bios2mds)

# Load DNAStringSet file from 'Loading_data' script
seq_string <- readDNAStringSet(filepath = "./covseq_DNAStringSet")

# Multiple alignments with MUSCLE - to be done in bash script
seq_align <- muscle::muscle(stringset= seq_string, quiet=T)
seq_align 

# Save alignment output as fasta file
seq_align_as_align <- msaConvert(seq_align, "bios2mds::align") #convert to align object
export.fasta(seq_align_as_align, outfile = "cov_alignment.fasta")