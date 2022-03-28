#downloading all the packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("muscle")
library(muscle)

BiocManager::install("annotate")
library(annotate)

BiocManager::install("Biostrings")
library(Biostrings)

a#installing other packages
install.packages("seqinr")
library(seqinr)

install.packages("ape")
library(ape)

install.packages("dplyr")
library(dplyr)

install.packages("ggplot2")
library(ggplot2)

#Loading and naming the files
HCov=read.FASTA("HCov-HKU1_NC_006577", type = "DNA")
MERS=read.FASTA("MERS-CoV_NC_019843", type = "DNA")
SARS=read.FASTA("SARS-Cov_NC_004718", type = "DNA")

#converting to strings

#HCovstring
HCovstring=HCov %>% 
  as.character %>% 
  lapply(.,paste0,collapse="") %>%
  unlist %>%
  DNAStringSet
head(HCovstring)
  
#MERSstring

MERSstring=MERS %>% 
  as.character %>% 
  lapply(.,paste0,collapse="") %>%
  unlist %>%
  DNAStringSet
head(MERSstring)

#MERSstring

SARSstring=SARS %>% 
  as.character %>% 
  lapply(.,paste0,collapse="") %>%
  unlist %>%
  DNAStringSet
head(SARSstring)

