#Loading Data for MUSCLE alignment

#Load packages
library(ape)
install.packages("EnvNJ")
library(EnvNJ)
library(Biostrings)
library(dplyr)

# Read FASTA files
files <- list.files(path="./Sequences", pattern = ".fasta", 
                    full.names=FALSE, recursive=FALSE) #Rename all files to end with '.fasta'

filenames <- regmatches(files, regexpr(".*[^.fasta]", files)) #remove .fasta extension

#concatenate all sequences into a single .fasta file  for MUSCLE alignment
fastaconc(filenames, inputdir = "./Sequences", out.file = "./covseq_concatenated_Heather.fasta") 

#Read concatenated fasta file with read.FASTA
seq <- read.FASTA("covseq_concatenated_Heather.fasta") 
seq

# Convert DNASbin to DNAStringSet 
seq_string <- seq %>% as.character %>% lapply(.,paste0,collapse="") %>%
  unlist %>% DNAStringSet

#save as a FASTA/StringSet file
writeXStringSet(seq_string, filepath = "./covseq_DNAStringSet", 
                format = "fasta") 

# Identify gaps in sequence lengths (if present)
SeqLen <- as.numeric(lapply(seq_string, length))
library(ggplot2)
qplot(SeqLen) + theme_bw() #visualize gaps
