## Multiple Alignments with MUSCLE ##

# Packages
library(ape)
library(Biostrings)
library(annotate)
library(muscle)
library(bios2mds)
library(msa)

# Load DNAStringSet file
seq_string <- readDNAStringSet(filepath = "./covseq_DNAStringSet")

# Multiple alignments with MUSCLE
seq_align <- muscle::muscle(stringset= seq_string, quiet=T)

# Save alignment output as fasta file
seq_align_as_align <- msaConvert(seq_align, "bios2mds::align") #convert to align object
export.fasta(seq_align_as_align, outfile = "cov_alignment.fasta")