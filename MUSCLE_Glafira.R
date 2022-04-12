#!/usr/bin/env Rscript


library(ape)
library(Biostrings)
library(annotate)
library(muscle)
library(bios2mds)

# Load DNAStringSet file
seq_string <- readDNAStringSet(filepath = "./covseq_DNAStringSet_GE")

# Multiple alignments with MUSCLE
seq_align <- muscle::muscle(stringset= seq_string, quiet=T)
seq_align

# Save alignment output as fasta file
seq_align_as_align <- msaConvert(seq_align, "bios2mds::align") #convert to align object
export.fasta(seq_align_as_align, outfile = "cov_alignment_GE.fasta")

########################

args <- commandArgs(trailingOnly = TRUE)
cat(args, sep = "\n")

seq_align=args[1]

