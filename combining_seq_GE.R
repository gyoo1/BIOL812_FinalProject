#uploading files
#downloading all the packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("muscle")
library(muscle)

BiocManager::install("annotate")
library(annotate)

BiocManager::install("Biostrings")
library(Biostrings)

#installing other packages
install.packages("seqinr")
library(seqinr)

install.packages("ape")
library(ape)

install.packages("dplyr")
library(dplyr)

install.packages("ggplot2")
library(ggplot2)

install.packages(reshape2)
library(reshape2)

install.packages(readr)
library(readr)

install.packages(stringr)
library(stringr)

install.packages(littler)
library(littler)

install.packages(EnvNJ)
library(EnvNJ)
##########################################################
#Loading and naming the files

files <- list.files(path="./Sequences", pattern = ".fasta", 
                    full.names=FALSE, recursive=FALSE) #get list of file names from directory

filenames <- regmatches(files, regexpr(".*[^.fasta]", files)) #remove .fasta extension

fastaconc(filenames, inputdir = "./Sequences", out.file = "./covseq_concatenated.fasta") 

#concatenate all sequences into a single .fasta file

DNAseq <- read.FASTA("covseq_concatenated.fasta")

# Convert DNAseq to DNA_string file
DNA_string = DNAseq %>% as.character %>% lapply(.,paste0,collapse="") %>%
  unlist %>% DNAStringSet

writeXStringSet(DNA_string, filepath = "./covseq_DNAStringSet_GE", 
                format = "fasta") #save it as fasta file