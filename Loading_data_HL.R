#Loading Data for MUSCLE alignment

#Load packages
library(ape)
library(EnvNJ)
library(Biostrings)
library(dplyr)

#Read concatenated fasta file (from python script) with read.FASTA

seq <- read.FASTA("./Sequences/Combined_fasta.fasta") #read file containing combines fasta sequences

# Convert DNAbin to DNAStringSet 
seq_string <- seq %>% as.character %>% lapply(.,paste0,collapse="") %>%
  unlist %>% DNAStringSet

#save as a FASTA/StringSet file
writeXStringSet(seq_string, filepath = "./covseq_DNAStringSet", 
                format = "fasta") 

# Identify gaps in sequence lengths (if present)
SeqLen<- as.numeric(lapply(seq_string, length))
library(ggplot2)
qplot(SeqLen) + theme_bw() #visualize gaps

