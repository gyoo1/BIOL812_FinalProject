## FASTA to DNAStringSet ##

# Packages
library(ape)
library(EnvNJ)
library(Biostrings)
library(dplyr)

# Read FASTA files
files <- list.files(path="./Sequences", pattern = ".fasta", 
                    full.names=FALSE, recursive=FALSE) #get list of file names from directory

filenames <- regmatches(files, regexpr(".*[^.fasta]", files)) #remove .fasta extension

fastaconc(filenames, inputdir = "./Sequences", out.file = "./covseq_concatenated.fasta") 
#concatenate all sequences into a single .fasta file

seq <- read.FASTA("covseq_concatenated.fasta") #read concatenated sequence


# Convert DNASbin to DNAStringSet file
seq_string <- seq %>% as.character %>% lapply(.,paste0,collapse="") %>%
        unlist %>% DNAStringSet

writeXStringSet(seq_string, filepath = "./covseq_DNAStringSet", 
                format = "fasta") #save as fasta file

# Identify gaps in sequence lengths
SeqLen <- as.numeric(lapply(seq_string, length))
library(ggplot2)
qplot(SeqLen) + theme_bw() #visualize gaps
