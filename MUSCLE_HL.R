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
library(genbankr)
library(ggtree)
library(annotate)
library(ape)
?read.FASTA()
library(muscle)
library(dplyr)
library(Biostrings)
?muscle
#Load Fasta Files as variables
HCoV<- read.FASTA("HCoV-HKU1_NC_006577", type="DNA")
MERS<- read.FASTA("MERS-CoV_NC_019843", type="DNA")
SARS<- read.FASTA("SARS-CoV_NC_004718", type="DNA")

#Convert to strings
HCoVstring<-HCoV$Sew %>%
  as.character %>%
  lapply(.,paste0,collapse="") %>%
  unlist%>%
  DNAStringSet

MERSstring<-MERS$Sew %>%
  as.character %>%
  lapply(.,paste0,collapse="") %>%
  unlist%>%
  DNAStringSet

SARSstring<-SARS$Sew %>%
  as.character %>%
  lapply(.,paste0,collapse="") %>%
  unlist%>%
  DNAStringSet
# 1st Commit 

